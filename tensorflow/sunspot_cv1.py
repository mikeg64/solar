#

import os
import numpy as np
import matplotlib.pyplot as plt
import keras
from keras.models import load_model
from keras.callbacks import Callback
import time

import models
import data_process_256

from keras.utils import multi_gpu_model
from keras.callbacks import ModelCheckpoint
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix
import tensorflow as tf
from numpy import save
# from clr_callback import CyclicLR
from keras import backend as K

exp_id = 'cv1'
data_type = 'm'  # c, m, combo - c for continuum, m for magnetogram, combo for combining c & m
model_name = 'inceptionResNetV2'  # cnnV1, inceptionResNetV2, nasNetLarge, lstm, resNet50
batch_size = 64  # 128 with 1 GPU OOM!
epoch = 200  # expected 90hrs
seq_length = 7  # will be ignored if is not lstm model
balance_data = 'max'  # If balance input data to the smallest/largest class, use 'min'/'max'.
desired_train_num = None  # How many times desired train number of the largest class train number
                          # Sample number in each class = original(~1900) + augmentated data
                          # None or a number > 1
num_gpu = 1
seed = 394
fold_num = 5
print("exp_id: ", exp_id)
print("batch_size: ", batch_size)
print("data_type: ", data_type)
print("seed: ", seed)
print("fold_num: ", fold_num)


class CyclicLR(Callback):
    """This callback implements a cyclical learning rate policy (CLR).
    The method cycles the learning rate between two boundaries with
    some constant frequency, as detailed in this paper (https://arxiv.org/abs/1506.01186).
    The amplitude of the cycle can be scaled on a per-iteration or
    per-cycle basis.
    This class has three built-in policies, as put forth in the paper.
    "triangular":
        A basic triangular cycle w/ no amplitude scaling.
    "triangular2":
        A basic triangular cycle that scales initial amplitude by half each cycle.
    "exp_range":
        A cycle that scales initial amplitude by gamma**(cycle iterations) at each
        cycle iteration.
    For more detail, please see paper.

    # Example
        ```python
            clr = CyclicLR(base_lr=0.001, max_lr=0.006,
                                step_size=2000., mode='triangular')
            model.fit(X_train, Y_train, callbacks=[clr])
        ```

    Class also supports custom scaling functions:
        ```python
            clr_fn = lambda x: 0.5*(1+np.sin(x*np.pi/2.))
            clr = CyclicLR(base_lr=0.001, max_lr=0.006,
                                step_size=2000., scale_fn=clr_fn,
                                scale_mode='cycle')
            model.fit(X_train, Y_train, callbacks=[clr])
        ```
    # Arguments
        base_lr: initial learning rate which is the
            lower boundary in the cycle.
        max_lr: upper boundary in the cycle. Functionally,
            it defines the cycle amplitude (max_lr - base_lr).
            The lr at any cycle is the sum of base_lr
            and some scaling of the amplitude; therefore
            max_lr may not actually be reached depending on
            scaling function.
        step_size: number of training iterations per
            half cycle. Authors suggest setting step_size
            2-8 x training iterations in epoch.
        mode: one of {triangular, triangular2, exp_range}.
            Default 'triangular'.
            Values correspond to policies detailed above.
            If scale_fn is not None, this argument is ignored.
        gamma: constant in 'exp_range' scaling function:
            gamma**(cycle iterations)
        scale_fn: Custom scaling policy defined by a single
            argument lambda function, where
            0 <= scale_fn(x) <= 1 for all x >= 0.
            mode paramater is ignored
        scale_mode: {'cycle', 'iterations'}.
            Defines whether scale_fn is evaluated on
            cycle number or cycle iterations (training
            iterations since start of cycle). Default is 'cycle'.
    """

    def __init__(self, base_lr=0.001, max_lr=0.006, step_size=2000., mode='triangular',
                 gamma=1., scale_fn=None, scale_mode='cycle'):
        super(CyclicLR, self).__init__()

        self.base_lr = base_lr
        self.max_lr = max_lr
        self.step_size = step_size
        self.mode = mode
        self.gamma = gamma
        if scale_fn == None:
            if self.mode == 'triangular':
                self.scale_fn = lambda x: 1.
                self.scale_mode = 'cycle'
            elif self.mode == 'triangular2':
                self.scale_fn = lambda x: 1 / (2. ** (x - 1))
                self.scale_mode = 'cycle'
            elif self.mode == 'exp_range':
                self.scale_fn = lambda x: gamma ** (x)
                self.scale_mode = 'iterations'
        else:
            self.scale_fn = scale_fn
            self.scale_mode = scale_mode
        self.clr_iterations = 0.
        self.trn_iterations = 0.
        self.history = {}

        self._reset()

    def _reset(self, new_base_lr=None, new_max_lr=None,
               new_step_size=None):
        """Resets cycle iterations.
        Optional boundary/step size adjustment.
        """
        if new_base_lr != None:
            self.base_lr = new_base_lr
        if new_max_lr != None:
            self.max_lr = new_max_lr
        if new_step_size != None:
            self.step_size = new_step_size
        self.clr_iterations = 0.

    def clr(self):
        cycle = np.floor(1 + self.clr_iterations / (2 * self.step_size))
        x = np.abs(self.clr_iterations / self.step_size - 2 * cycle + 1)
        if self.scale_mode == 'cycle':
            return self.base_lr + (self.max_lr - self.base_lr) * np.maximum(0, (1 - x)) * self.scale_fn(cycle)
        else:
            return self.base_lr + (self.max_lr - self.base_lr) * np.maximum(0, (1 - x)) * self.scale_fn(
                self.clr_iterations)

    def on_train_begin(self, logs={}):
        logs = logs or {}

        if self.clr_iterations == 0:
            K.set_value(self.model.optimizer.lr, self.base_lr)
        else:
            K.set_value(self.model.optimizer.lr, self.clr())

    def on_batch_end(self, epoch, logs=None):

        logs = logs or {}
        self.trn_iterations += 1
        self.clr_iterations += 1

        self.history.setdefault('lr', []).append(K.get_value(self.model.optimizer.lr))
        self.history.setdefault('iterations', []).append(self.trn_iterations)

        for k, v in logs.items():
            self.history.setdefault(k, []).append(v)

        K.set_value(self.model.optimizer.lr, self.clr())


# Get f1_score as ModelCheckPoint metric
class Metrics(tf.keras.callbacks.Callback):
    def __init__(self, valid_data):
        super(Metrics, self).__init__()
        self.validation_data = valid_data

    def on_epoch_end(self, epoch, logs=None):
        logs = logs or {}
        val_predict = np.argmax(self.model.predict(self.validation_data[0]), -1)
        val_targ = self.validation_data[1]
        if len(val_targ.shape) == 2 and val_targ.shape[1] != 1:
            val_targ = np.argmax(val_targ, -1)

        _val_f1 = f1_score(val_targ, val_predict, average=None)
        _val_f1 = _val_f1[1]  # [1] corresponds to the sunspot class beta

        logs['val_f1'] = _val_f1
        print("val_f1: %f" % _val_f1)
        return


print("Start running!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
# Load data
data_c, data_m, data_y = data_process_256.data_preprosess_chunks(model_name, balance_data, desired_train_num,
                                                                 seq_length, seed, fold_num)

# Combine continuum and magnetogram data
data_combo = [[] for i in range(fold_num)]
if model_name in ('inceptionResNetV2', 'nasNetLarge', 'resNet50'):
    for n in range(fold_num):
        data_c_sc = np.expand_dims(data_c[n][:, :, :, 0], axis=3)
        data_m_sc = np.expand_dims(data_m[n][:, :, :, 0], axis=3)
        data_combo[n] = np.concatenate([data_c_sc, data_m_sc, data_m_sc], axis=3)
else:
    for n in range(fold_num):
        data_combo[n] = np.concatenate([data_c[n], data_m[n]], axis=3)
print("data_combo[0] shape: ", data_combo[0].shape)

# Define train/test data to continuum, magnetogram or combo
if data_type == 'c':
    data_x = data_c
elif data_type == 'm':
    data_x = data_m
elif data_type == 'combo':
    data_x = data_combo
print("data_x[0] shape: ", data_x[0].shape)
print("data_y[0] shape: ", data_y[0].shape)

# if num_gpu > 1:
#     if balance_data == 'min':
#         batch_size = int(np.ceil(train_x.shape[0] / 4))  # /7 for Conv2D kernal_size=11
#     elif balance_data == 'max':
#         batch_size = int(np.ceil(train_x.shape[0] / 14))
print("batch_size: ", batch_size)

fold_f1_macro = []
fold_f1_micro = []
fold_f1_weighted = []
fold_f1_single_class = []
fold_val_acc = []
fold_alpha_acc = []
fold_beta_acc = []
fold_betax_acc = []
for n in range(fold_num):
    start_time = time.time()
    print("This is fold {}".format(n))

    multiGPU_model_path = './history/' + str(exp_id) + '_multiGPU_fold' + str(n) + '.h5'
    singleGPU_model_path = './history/' + str(exp_id) + '_singleGPU_fold' + str(n) + '.h5'
    predict_result = './history/' + str(exp_id) + '_pred_fold' + str(n) + '.npy'

    train_ind = np.where(np.arange(fold_num) != n)
    train_x = np.concatenate((data_x[train_ind]))
    train_y = np.concatenate((np.array(data_y)[train_ind]))
    test_x = data_x[n]
    test_y = data_y[n]
    print("train_x.shape: ", train_x.shape)
    print("train_y.shape: ", train_y.shape)
    print("test_x.shape: ", test_x.shape)
    print("test_y.shape: ", test_y.shape)

    # Construct model
    if model_name == 'cnnV1':
        my_model = models.cnn_model_v1(train_x[0].shape)
    elif model_name == 'inceptionResNetV2':
        my_model = models.inceptionResNetV2(train_x[0].shape)
    elif model_name == 'nasNetLarge':
        my_model = models.nasNetLarge(train_x[0].shape)
    elif model_name == 'lstm':
        my_model = models.lstm(train_x[0].shape)
    elif model_name == 'resNet50':
        my_model = models.resNet50(train_x[0].shape)

    # Multiple GPUs
    if num_gpu > 1:  # If num_gpu=1, the function multi_gpu_model gives error!
        my_model_gpu = multi_gpu_model(my_model, gpus=num_gpu)
    elif num_gpu == 1:
        my_model_gpu = my_model
    # my_model_gpu.summary()

    my_model_gpu.compile(loss='categorical_crossentropy', optimizer=keras.optimizers.Adam(), metrics=['accuracy'])
    for k, layer in enumerate(my_model_gpu.layers):
        print(k, layer.name, layer.trainable)

    # initialize the cyclical learning rate callback
    min_lr = 1e-7
    max_lr = 1e-2
    step_size = 8
    clr_method = "triangular"
    NUM_EPOCHS = 96
    clr = CyclicLR(mode=clr_method, base_lr=min_lr, max_lr=max_lr,
                   step_size=step_size * (train_x.shape[0] // batch_size))

    # checkpoint callback
    checkpoint = ModelCheckpoint(multiGPU_model_path, monitor='val_acc', verbose=1, save_best_only=True, mode='max')

    # Callback list
    callbacks_list = [Metrics(valid_data=(test_x, test_y)), checkpoint, clr]

    # fit model
    history = my_model_gpu.fit(train_x, train_y, validation_data=(test_x, test_y), batch_size=batch_size, epochs=epoch,
                               callbacks=callbacks_list, verbose=1)
    fold_val_acc.append(history.history['val_acc'])

    # plot the training loss and accuracy
    N = np.arange(0, epoch)
    plt.style.use("ggplot")
    plt.figure()
    plt.plot(N, history.history["loss"], label="train_loss")
    plt.plot(N, history.history["val_loss"], label="val_loss")
    plt.plot(N, history.history["acc"], label="train_acc")
    plt.plot(N, history.history["val_acc"], label="val_acc")
    plt.title("Training Loss and Accuracy")
    plt.xlabel("Epoch #")
    plt.ylabel("Loss/Accuracy")
    plt.legend(loc="lower left")
    plt.savefig('./history/' + str(exp_id) + '_training_plot_fold' + str(n) + '.png')
    # plot the learning rate history
    N = np.arange(0, len(clr.history["lr"]))
    plt.figure()
    plt.plot(N, clr.history["lr"])
    plt.title("Cyclical Learning Rate (CLR)")
    plt.xlabel("Training Iterations")
    plt.ylabel("Learning Rate")
    plt.savefig('./history/' + str(exp_id) + '_clr_plot_fold' + str(n) + '.png')

    # Evaluate the best model
    del my_model_gpu
    my_model_gpu = load_model(multiGPU_model_path)
    if num_gpu > 1:  # If num_gpu=1, the function multi_gpu_model gives error!
        my_model = my_model_gpu.layers[-2]  # you can use multi_gpus_model.summary() to see the layer of the original model
    elif num_gpu == 1:
        my_model = my_model_gpu
    my_model.save_weights(singleGPU_model_path)

    y_pred = my_model_gpu.predict(test_x)
    print("y_pred shape: ", y_pred.shape)
    print("test_y shape: ", test_y.shape)
    save(predict_result, y_pred)  # Save predict result to .npy for later use

    # Convert one-hot and continuous multi-labels to numerical multi-labels
    y_pred = np.argmax(y_pred, axis=1)
    test_y = np.argmax(test_y, axis=1)

    # Calculate acc
    cm = confusion_matrix(test_y, y_pred)
    # Now the normalize the diagonal entries to get the accuracy in percentage
    cm_percent = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

    # Save fold acc for each class
    fold_alpha_acc.append(cm_percent[0][0])
    fold_beta_acc.append(cm_percent[1][1])
    fold_betax_acc.append(cm_percent[2][2])

    # Calculate fold F1 scores
    f1_macro = f1_score(test_y, y_pred, average='macro')
    f1_micro = f1_score(test_y, y_pred, average='micro')
    f1_weighted = f1_score(test_y, y_pred, average='weighted')
    f1_single_class = f1_score(test_y, y_pred, average=None)

    # Save fold F1 scores
    fold_f1_macro.append(f1_macro)
    fold_f1_micro.append(f1_micro)
    fold_f1_weighted.append(f1_weighted)
    fold_f1_single_class.append(f1_single_class)

    # dat = np.array([cm, cm_percent, f1_macro, f1_micro, f1_weighted, f1_single_class], dtype=object)
    # np.savetxt(os.getcwd() + '/history/' + str(exp_id) + '_report_' + str(model_name) + '_' + str(batch_size) + '_' + balance_data + '_'
    #            + str(seed) + '_' + str(seq_length) + '_fold' + str(n) + '.txt', dat, fmt='%s')

    K.clear_session()

    print("Fold {} training complete!".format(n))
    print("--- %s seconds ---" % (time.time() - start_time))

# Calculated mean single class acc over folders
ave_alpha_acc = np.mean(fold_alpha_acc)
ave_beta_acc = np.mean(fold_beta_acc)
ave_betax_acc = np.mean(fold_betax_acc)

# Save acc and F1.
# Don't receive fold_val_acc here becoz it's as long as epoch number and would be too long.
dat = np.array([ave_alpha_acc, ave_beta_acc, ave_betax_acc, fold_alpha_acc, fold_beta_acc, fold_betax_acc,
                fold_f1_macro, fold_f1_micro, fold_f1_weighted, fold_f1_single_class], dtype=object)
np.savetxt(os.getcwd() + '/history/' + str(exp_id) + '_acc_F1_' + str(model_name) + '_' + str(batch_size) + '_'
           + balance_data + '_' + str(seed) + '_' + str(seq_length) + '.txt', dat, fmt='%s')

