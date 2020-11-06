# train/test ratio: 0.85/0.15

import os
from skimage.color import rgb2gray
from keras.preprocessing.image import img_to_array, load_img
from datetime import datetime

import numpy as np
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import random
from keras import backend as K
from collections import defaultdict
import matplotlib.gridspec as gridspec
# from astropy.io import fits
# from sunpy.io import read_file
from scipy.io.idl import readsav
# import cv2
# from skimage.transform import resize
from keras.utils import to_categorical
from keras.preprocessing.image import ImageDataGenerator


def load_validation_data(model_name):
    data_dir = os.path.join(os.getcwd() + '/data/test_data/square256')

    valid_c = []
    objs_c = []
    data_type_clas = os.path.join(data_dir + '/continuum')
    # Load images
    image_names = [f for f in sorted(os.listdir(data_type_clas)) if f.endswith(".sav")]
    # print("image_names: ", image_names)
    for f in image_names:
        image = readsav(os.path.join(data_type_clas, f))
        image = np.array(image['data'])
        # image = resize(image, (32, 32), anti_aliasing=True)
        image = np.clip(image, 2e4, 6e4)
        # plt.imshow(image, cmap='gray')
        # plt.show()
        valid_c.append(image)
        objs_c.append(int(f.split('.')[3][-8:-4]))

    valid_m = []
    objs_m = []
    data_type_clas = os.path.join(data_dir + '/magnetogram')
    # Load images
    image_names = [f for f in sorted(os.listdir(data_type_clas)) if f.endswith(".sav")]
    for f in image_names:
        image = readsav(os.path.join(data_type_clas, f))
        image = np.array(image['data'])
        # image = resize(image, (32, 32), anti_aliasing=True)
        image = np.clip(image, -800, 800)
        # plt.imshow(image, cmap='gray')
        # plt.show()
        valid_m.append(image)
        objs_m.append(int(f.split('.')[3][-8:-4]))

    if objs_c != objs_m:
        raise Exception('Files have not been read in properly!')

    # Preprocess data
    # Normalise to [0,1]
    valid_c = (np.array(valid_c) - 2e4) / (6e4 - 2e4)
    valid_m = (np.array(valid_m) - (-800)) / (800 - (-800))
    # Expand to have color channel
    valid_c = np.expand_dims(valid_c, 3)
    valid_m = np.expand_dims(valid_m, 3)

    if model_name in ('inceptionResNetV2', 'inceptionResNetV2_2class', 'nasNetLarge', 'resNet50', 'xCeption',
                      'inceptionResNetV2_plus_LSTM', 'mm_tseng'):
        valid_c = np.repeat(valid_c, 3, -1)
        valid_m = np.repeat(valid_m, 3, -1)

    return valid_c, valid_m, objs_c


def load_data():
    data_dir = os.path.join(os.getcwd() + '/data/square256')

    images_c = [[]]
    objs_c = [[]]
    for ind, clas in enumerate(['/alpha', '/beta', '/betax']):
        data_type_clas = os.path.join(data_dir + '/continuum' + clas)
        if images_c[-1]:
            images_c.append([])
            objs_c.append([])

        # Load images
        image_names = [f for f in sorted(os.listdir(data_type_clas)) if f.endswith(".sav")]
        # print("image_names: ", image_names)
        for f in image_names:
            # print("Continuum, %s, %s" % (clas, f))
            # image = read_file(os.path.join(data_type_clas, f))[1][0]
            # image = resize(image, (256, 256), anti_aliasing=True)
            image = readsav(os.path.join(data_type_clas, f))
            image = np.array(image['data'])
            # image = resize(image, (32, 32), anti_aliasing=True)
            # plt.imshow(image, cmap='gray')
            # plt.show()
            image = np.clip(image, 2e4, 6e4)
            # plt.imshow(image, cmap='gray')
            # plt.show()
            # image = np.log(image)
            # plt.imshow(image, cmap='gray')
            # plt.show()

            # image = cv2.imread(os.path.join(label_directory, f))
            # image = cv2.resize(image, (pixels, pixels))
            # image = cv2.medianBlur(image, 5)
            images_c[ind].append(image)
            objs_c[ind].append(int(f.split('.')[2]))

    images_m = [[]]
    objs_m = [[]]
    for ind, clas in enumerate(['/alpha', '/beta', '/betax']):
        data_type_clas = os.path.join(data_dir + '/magnetogram' + clas)
        if images_m[-1]:
            images_m.append([])
            objs_m.append([])

        # Load images
        image_names = [f for f in sorted(os.listdir(data_type_clas)) if f.endswith(".sav")]
        # print("image_names: ", image_names)
        for f in image_names:
            # print("Magnetogram, %s, %s" % (clas, f))
            image = readsav(os.path.join(data_type_clas, f))
            image = np.array(image['data'])
            # image = resize(image, (32, 32), anti_aliasing=True)
            image = np.clip(image, -800, 800)
            # plt.imshow(image, cmap='gray')
            # plt.show()
            images_m[ind].append(image)
            objs_m[ind].append(int(f.split('.')[2]))

    if objs_c != objs_m:
        raise Exception('Files have not been read in properly!')

    # # Create y labels
    # y = []
    # # alpha = np.tile([0, 0, 1], (clas_num[0], 1))
    # # beta = np.tile([0, 1, 0], (clas_num[1], 1))
    # # betax = np.tile([1, 0, 0], (clas_num[2], 1))
    # alpha = np.full(clas_num[0], 0)
    # beta = np.full(clas_num[1], 1)
    # betax = np.full(clas_num[2], 2)
    # y.append(alpha)
    # y.append(beta)
    # y.append(betax)

    return images_c, images_m, objs_c


def sample_data_min(_seed):
    # Load all data from directories
    images_c, images_m, objs_c = load_data()

    # Generate train/test datasets
    train_c = []
    train_m = []
    train_y = []
    test_c = []
    test_m = []
    test_y = []
    # sample_num = min(clas_num)
    sample_num = min(len(objs_c[0]), len(objs_c[1]), len(objs_c[2]))
    train_num = int(np.floor(sample_num * 0.85))
    test_num = int(np.floor(sample_num * 0.15))
    print("train_num: ", train_num)
    print("test_num: ", test_num)

    for i in range(3):  # 3 classes
        train_count = 0  # Number of images in a class
        test_count = 0
        events = list(set(objs_c[i]))
        random.Random(_seed).shuffle(events)
        print("shuffled events: ", events)

        for event in events:
            # Have enough test data, start to load train data
            if test_count >= test_num:
                for ind, obj in enumerate(objs_c[i]):
                    if obj == event:
                        train_c.append(images_c[i][ind])
                        train_m.append(images_m[i][ind])
                        train_y.append(i)
                        train_count += 1
                if train_count >= train_num:
                    break
                else:
                    continue

            # Load test data
            for ind, obj in enumerate(objs_c[i]):
                if obj == event:
                    test_c.append(images_c[i][ind])
                    test_m.append(images_m[i][ind])
                    test_y.append(i)
                    test_count += 1

    clas_num = [len(objs_c[0]), len(objs_c[1]), len(objs_c[2])]
    print("clas_num: ", clas_num)
    print("train data class 0 num: ", train_y.count(0))
    print("train data class 1 num: ", train_y.count(1))
    print("train data class 2 num: ", train_y.count(2))
    print("test data class 0 num: ", test_y.count(0))
    print("test data class 1 num: ", test_y.count(1))
    print("test data class 2 num: ", test_y.count(2))

    return np.array(train_c), np.array(train_m), train_y, np.array(test_c), np.array(test_m), test_y


def sample_data_max(desired_train_num, _seed):
    images_c, images_m, objs_c = load_data()

    # Calculate test_num
    clas_num = [len(objs_c[0]), len(objs_c[1]), len(objs_c[2])]
    smallest_clas_num = min(clas_num)
    test_num = int(np.floor(smallest_clas_num * 0.15))

    # Generate train/test datasets
    train_c = [[]]
    train_m = [[]]
    train_y = [[]]
    test_c = []
    test_m = []
    test_y = []
    for i in range(3):  # 3 classes
        if train_c[-1]:
            train_c.append([])
            train_m.append([])
            train_y.append([])

        test_count = 0
        events = list(set(objs_c[i]))
        random.Random(_seed).shuffle(events)
        print("shuffled events: ", events)

        for event in events:
            # Have enough test data, start to load all the remaining data to training data
            if test_count >= test_num:
                for ind, obj in enumerate(objs_c[i]):
                    if obj == event:
                        train_c[i].append(images_c[i][ind])
                        train_m[i].append(images_m[i][ind])
                        train_y[i].append(i)
                continue

            # Load test data
            for ind, obj in enumerate(objs_c[i]):
                if obj == event:
                    test_c.append(images_c[i][ind])
                    test_m.append(images_m[i][ind])
                    test_y.append(i)
                    test_count += 1

    # Calculate actual train_num
    train_clas_num = [len(train_c[0]), len(train_c[1]), len(train_c[2])]
    if desired_train_num is None:
        train_num = max(train_clas_num)
    else:
        train_num = int(desired_train_num * max(train_clas_num))

    # Augment training data
    datagen = ImageDataGenerator(
        featurewise_center=True,
        featurewise_std_normalization=True,
        zca_whitening=True,
        rotation_range=20,
        width_shift_range=0.2,
        height_shift_range=0.2,
        shear_range=20,
        zoom_range=0.2,
        horizontal_flip=True,
        vertical_flip=True,
        fill_mode='nearest')

    for i in range(3):
        train_c[i] = np.expand_dims(train_c[i], 3)
        train_m[i] = np.expand_dims(train_m[i], 3)
        train_y[i] = np.array(train_y[i])
        if len(train_c[i]) < train_num:
            for batch_c, batch_y in datagen.flow(train_c[i], train_y[i], batch_size=32):
                train_c[i] = np.concatenate((train_c[i], batch_c))
                if len(train_c[i]) >= train_num:
                    break
            for batch_m, batch_y in datagen.flow(train_m[i], train_y[i], batch_size=32):
                train_m[i] = np.concatenate((train_m[i], batch_m))
                train_y[i] = np.concatenate((train_y[i], batch_y))
                if len(train_m[i]) >= train_num:
                    break

    # # Visualise training dataset after augmentation
    # for i in range(68):
    #     fig, axs = plt.subplots(1,2)
    #     axs[0].imshow(np.squeeze(train_c[1][i]), cmap='gray')
    #     axs[0].set_title('continuum')
    #     axs[1].imshow(np.squeeze(train_m[1][i]), cmap='gray')
    #     axs[1].set_title('magnetogram')
    #     fig.suptitle('beta')
    #     plt.show()
    return np.concatenate(train_c), np.concatenate(train_m), np.concatenate(train_y), np.array(test_c), np.array(test_m), test_y


def sample_data_all_available(_seed):
    images_c, images_m, objs_c = load_data()

    # Calculate test_num
    clas_num = [len(objs_c[0]), len(objs_c[1]), len(objs_c[2])]
    smallest_clas_num = min(clas_num)
    test_num = int(np.floor(smallest_clas_num * 0.15))
    print("test_num: ", test_num)

    # Generate train/test datasets
    train_c = []
    train_m = []
    train_y = []
    test_c = []
    test_m = []
    test_y = []
    for i in range(3):  # 3 classes
        test_count = 0
        events = list(set(objs_c[i]))
        random.Random(_seed).shuffle(events)
        print("shuffled events: ", events)

        for event in events:
            # Have enough test data, start to load train data
            if test_count >= test_num:
                for ind, obj in enumerate(objs_c[i]):
                    if obj == event:
                        train_c.append(images_c[i][ind])
                        train_m.append(images_m[i][ind])
                        train_y.append(i)
                continue

            # Load test data
            for ind, obj in enumerate(objs_c[i]):
                if obj == event:
                    test_c.append(images_c[i][ind])
                    test_m.append(images_m[i][ind])
                    test_y.append(i)
                    test_count += 1

    # # Visualise training dataset after augmentation
    # for i in range(68):
    #     fig, axs = plt.subplots(1,2)
    #     axs[0].imshow(train_c[i], cmap='gray')
    #     axs[0].set_title('continuum')
    #     axs[1].imshow(train_m[i], cmap='gray')
    #     axs[1].set_title('magnetogram')
    #     fig.suptitle('beta')
    #     plt.show()

    return np.array(train_c), np.array(train_m), train_y, np.array(test_c), np.array(test_m), test_y


def load_data_lstm(seq_length, _seed):
    images_c, images_m, objs_c = load_data()

    # Stitch images into time sequence data format for lstm
    data_c, data_m, data_y, seq_events = construct_timeSeq_data(images_c, images_m, objs_c, seq_length)

    # Split sequence data into training/testing data - data_c, data_m and data_y using seq_events
    train_c = [[], [], []]
    train_m = [[], [], []]
    train_y = [[], [], []]
    test_c = [[], [], []]
    test_m = [[], [], []]
    test_y = [[], [], []]
    for i in range(3):
        seq_num = len(data_c[i])
        test_num = int(np.floor(seq_num * 0.15))
        print("class: ", i)
        print("seq_num: ", seq_num)
        print("test_num: ", test_num)
        test_count = 0

        events = list(set(seq_events[i]))
        random.Random(_seed).shuffle(events)
        print("shuffled events: ", events)

        for event in events:
            # Have enough test data, start to load all the remaining data to training data
            if test_count >= test_num:
                for ind, obj in enumerate(seq_events[i]):
                    if obj == event:
                        train_c[i].append(data_c[i][ind])
                        train_m[i].append(data_m[i][ind])
                        train_y[i].append(data_y[i][ind])
                continue

            # Load test data
            for ind, obj in enumerate(seq_events[i]):
                if obj == event:
                    test_c[i].append(data_c[i][ind])
                    test_m[i].append(data_m[i][ind])
                    test_y[i].append(data_y[i][ind])
                    test_count += 1

    train_c = np.concatenate(train_c)
    train_m = np.concatenate(train_m)
    train_y = np.concatenate(train_y).tolist()
    test_c = np.concatenate(test_c)
    test_m = np.concatenate(test_m)
    test_y = np.concatenate(test_y).tolist()

    print("train data class 0 num: ", train_y.count(0))
    print("train data class 1 num: ", train_y.count(1))
    print("train data class 2 num: ", train_y.count(2))
    print("test data class 0 num: ", test_y.count(0))
    print("test data class 1 num: ", test_y.count(1))
    print("test data class 2 num: ", test_y.count(2))

    #########################  Visualisation of data   ######################
    # fig = plt.figure(figsize=(70, 10))
    # for j in range(seq_length):
    #     ax = fig.add_subplot(1, 7, j + 1)
    #     ax.imshow(train_m[0][j], cmap='gray')
    #     ax.axis('off')
    # plt.show()
    #########################################################################

    return np.array(train_c), np.array(train_m), train_y, np.array(test_c), np.array(test_m), test_y


def sample_data_nfolds_nonTimeSeq(desired_train_num, _seed, fold_num, model_name):
    images_c, images_m, objs_c = load_data()

    # Calculate test_num
    clas1_num = int(np.floor(len(objs_c[0]) / fold_num))
    clas2_num = int(np.floor(len(objs_c[1]) / fold_num))
    clas3_num = int(np.floor(len(objs_c[2]) / fold_num))
    data_num = [clas1_num, clas2_num, clas3_num]

    # Generate n-fold dataset
    fold_c = [[] for j in range(fold_num)]
    fold_m = [[] for j in range(fold_num)]
    fold_y = [[] for j in range(fold_num)]
    for i in range(3):  # 3 classes
        n = 0
        data_count = 0
        events = list(set(objs_c[i]))
        random.Random(_seed).shuffle(events)
        print("shuffled events: ", events)

        for event in events:
            for ind, obj in enumerate(objs_c[i]):
                if obj == event:
                    fold_c[n].append(images_c[i][ind])
                    fold_m[n].append(images_m[i][ind])
                    fold_y[n].append(i)
                    data_count += 1
            if data_count >= data_num[i]:
                n += 1
                data_count = 0
            if n == fold_num:
                break

    for n in range(fold_num):
        print("data augmentation for fold: ", n)
        # Calculate desired class num in each fold
        actual_clas_num = [fold_y[n].count(0), fold_y[n].count(1), fold_y[n].count(2)]
        if desired_train_num is None:
            clas_num = max(actual_clas_num)
        else:
            clas_num = int(desired_train_num * max(actual_clas_num))
        print("clas_num: ", clas_num)

        # Augment training data
        datagen = ImageDataGenerator(
            featurewise_center=True,
            featurewise_std_normalization=True,
            zca_whitening=True,
            rotation_range=20,
            width_shift_range=0.2,
            height_shift_range=0.2,
            shear_range=20,
            zoom_range=0.2,
            horizontal_flip=True,
            vertical_flip=True,
            fill_mode='nearest')

        fold_c[n] = np.expand_dims(fold_c[n], 3)
        fold_m[n] = np.expand_dims(fold_m[n], 3)
        fold_y[n] = np.array(fold_y[n])
        fold_y_original = fold_y[n]
        for i in range(3):
            print("Class: ", i)
            print("np.count_nonzero(data_y[n]==i): ", np.count_nonzero(fold_y[n]==i))
            if np.count_nonzero(fold_y[n]==i) < clas_num:
                ind = np.where(fold_y[n] == i)
                for batch_c, batch_y in datagen.flow(fold_c[n][ind], fold_y[n][ind], batch_size=32):
                    fold_c[n] = np.concatenate((fold_c[n], batch_c))
                    fold_y[n] = np.concatenate((fold_y[n], batch_y))
                    if np.count_nonzero(fold_y[n]==i) >= clas_num:
                        break
                fold_y[n] = fold_y_original
                for batch_m, batch_y in datagen.flow(fold_m[n][ind], fold_y[n][ind], batch_size=32):
                    fold_m[n] = np.concatenate((fold_m[n], batch_m))
                    fold_y[n] = np.concatenate((fold_y[n], batch_y))
                    if np.count_nonzero(fold_y[n]==i) >= clas_num:
                        break
                fold_y_original = fold_y[n]
    print("Finish data augmentation.")

    # Normalise to [0,1]
    fold_c = (np.array(fold_c) - 2e4) / (6e4 - 2e4)
    fold_m = (np.array(fold_m) - (-800)) / (800 - (-800))
    for n in range(fold_num):
        print("data process for fold: ", n)
        fold_y[n] = to_categorical(fold_y[n])

        if model_name in ('inceptionResNetV2', 'inceptionResNetV2_2class', 'nasNetLarge', 'resNet50', 'xCeption',
                          'inceptionResNetV2_plus_LSTM', 'mm_tseng'):
            fold_c[n] = np.repeat(fold_c[n], 3, -1)
            # print("train_c shape: ", train_c.shape)  # (?, 256, 256, 3)
            fold_m[n] = np.repeat(fold_m[n], 3, -1)

        print("fold_c[" + str(n) + "].shape: ", fold_c[n].shape)
    print("len(fold_c): ", len(fold_c))

    return fold_c, fold_m, fold_y


def construct_timeSeq_data(images_c, images_m, objs_c, seq_length):
    # Generate all sequence data
    data_c = [[[]], [[]], [[]]]
    data_m = [[[]], [[]], [[]]]
    data_y = [[], [], []]
    seq_events = [[], [], []]
    pad_zeros = np.zeros((images_c[0][0].shape[0], images_c[0][0].shape[1]))  # init 2D array with 0s
    for i in range(3):
        ind_seq = 0
        events = list(set(objs_c[i]))

        for event in events:
            event_occurrence = objs_c[i].count(event)
            # If the obj doesn't have seq_length images
            if event_occurrence < seq_length:
                if data_c[i][-1]:
                    data_c[i].append([])
                    data_m[i].append([])
                for ind, obj in enumerate(objs_c[i]):
                    if obj == event:
                        data_c[i][ind_seq].append(images_c[i][ind])
                        data_m[i][ind_seq].append(images_m[i][ind])

                while len(data_c[i][ind_seq]) < seq_length:
                    data_c[i][ind_seq].append(pad_zeros)
                    data_m[i][ind_seq].append(pad_zeros)

                data_y[i].append(i)
                seq_events[i].append(event)
                ind_seq += 1
                continue
            # If obj has seq_length images
            for ind, obj in enumerate(objs_c[i]):
                if obj == event and ind + seq_length <= len(objs_c[i]):
                    if objs_c[i][ind + seq_length - 1] == event:
                        if data_c[i][-1]:
                            data_c[i].append([])
                            data_m[i].append([])
                        while len(data_c[i][ind_seq]) < seq_length:
                            data_c[i][ind_seq].append(images_c[i][ind])
                            data_m[i][ind_seq].append(images_m[i][ind])
                            ind += 1
                        data_y[i].append(i)
                        ind_seq += 1
                        seq_events[i].append(event)

        # Preprocess data
        data_c[i] = (np.array(data_c[i]) - 2e4) / (6e4 - 2e4)
        data_m[i] = (np.array(data_m[i]) - (-800)) / (800 - (-800))
        data_c[i] = np.expand_dims(data_c[i], axis=-1)
        data_m[i] = np.expand_dims(data_m[i], axis=-1)

    return np.array(data_c), np.array(data_m), np.array(data_y), seq_events


def sample_data_nfolds_timeSeq(desired_train_num, _seed, fold_num, model_name, seq_length):
    images_c, images_m, objs_c = load_data()

    # Stitch images into time sequence data format for lstm
    data_c, data_m, data_y, seq_events = construct_timeSeq_data(images_c, images_m, objs_c, seq_length)

    # Calculate fold num
    clas1_num = int(np.floor(len(seq_events[0]) / fold_num))
    clas2_num = int(np.floor(len(seq_events[1]) / fold_num))
    clas3_num = int(np.floor(len(seq_events[2]) / fold_num))
    data_num = [clas1_num, clas2_num, clas3_num]

    # Generate n-fold dataset
    fold_c = [[] for j in range(fold_num)]
    fold_m = [[] for j in range(fold_num)]
    fold_y = [[] for j in range(fold_num)]
    for i in range(3):  # 3 classes
        n = 0
        data_count = 0
        events = list(set(seq_events[i]))
        random.Random(_seed).shuffle(events)
        print("shuffled events: ", events)

        for event in events:
            for ind, obj in enumerate(seq_events[i]):
                if obj == event:
                    fold_c[n].append(data_c[i][ind])
                    fold_m[n].append(data_m[i][ind])
                    fold_y[n].append(i)
                    data_count += 1
            if data_count >= data_num[i]:
                n += 1
                data_count = 0
            if n == fold_num:
                break

    # Preprocess data
    for n in range(fold_num):
        print("data pre-process for fold: ", n)
        fold_y[n] = to_categorical(fold_y[n])

        if model_name in ('inceptionResNetV2', 'inceptionResNetV2_2class', 'nasNetLarge', 'resNet50', 'xCeption',
                          'inceptionResNetV2_plus_LSTM', 'mm_tseng'):
            fold_c[n] = np.repeat(fold_c[n], 3, -1)
            # print("train_c shape: ", train_c.shape)  # (?, 256, 256, 3)
            fold_m[n] = np.repeat(fold_m[n], 3, -1)

        print("fold_c[" + str(n) + "].shape: ", fold_c[n].shape)
    print("len(fold_c): ", len(fold_c))

    return np.array(fold_c), np.array(fold_m), np.array(fold_y)


def data_preprosess_chunks(model_name, balance_data, desired_train_num, seq_length, _seed, fold_num):
    print("model_name: ", model_name)

    if fold_num is None:  # train/test split
        if model_name in ('lstm', 'inceptionResNetV2_plus_LSTM', 'mm_tseng'):
            train_c, train_m, train_y, test_c, test_m, test_y = load_data_lstm(seq_length, _seed)
            train_c = np.expand_dims(train_c, 4)
            train_m = np.expand_dims(train_m, 4)
        else:
            if balance_data == 'min':
                train_c, train_m, train_y, test_c, test_m, test_y = sample_data_min(_seed)
                # Expand color channel to be 1
                train_c = np.expand_dims(train_c, 3)
                train_m = np.expand_dims(train_m, 3)
            elif balance_data == 'max':
                train_c, train_m, train_y, test_c, test_m, test_y = sample_data_max(desired_train_num, _seed)
            elif balance_data is None:
                train_c, train_m, train_y, test_c, test_m, test_y = sample_data_all_available(_seed)
                train_c = np.expand_dims(train_c, 3)
                train_m = np.expand_dims(train_m, 3)
    else:  # cross validation
        if model_name in ('lstm', 'inceptionResNetV2_plus_LSTM', 'mm_tseng'):
            fold_c, fold_m, fold_y = sample_data_nfolds_timeSeq(desired_train_num, _seed, fold_num, model_name, seq_length)
        else:
            fold_c, fold_m, fold_y = sample_data_nfolds_nonTimeSeq(desired_train_num, _seed, fold_num, model_name)

    # Pre-process data
    if fold_num is None:
        if model_name in ('cnnV1', 'cnnV2', 'lstm'):
            # Normalise to [0,1]
            train_c = (train_c - 2e4) / (6e4-2e4)
            # train_c = (train_c - np.log(2e4)) / (np.log(6e4) - np.log(2e4))
            train_m = (train_m - (-800)) / (800-(-800))
            test_c = (test_c - 2e4) / (6e4-2e4)
            # test_c = (test_c - np.log(2e4)) / (np.log(6e4) - np.log(2e4))
            test_m = (test_m - (-800)) / (800 - (-800))

            # Expand color channel to be 1
            train_y = to_categorical(train_y)
            test_c = np.expand_dims(test_c, 4)
            test_m = np.expand_dims(test_m, 4)
            test_y = to_categorical(test_y)
        elif model_name in ('inceptionResNetV2', 'inceptionResNetV2_2class', 'nasNetLarge', 'resNet50', 'xCeption',
                            'inceptionResNetV2_plus_LSTM', 'mm_tseng'):
            # Normalise to [0,1]
            train_c = (train_c - 2e4) / (6e4 - 2e4)
            # train_c = (train_c - np.log(2e4)) / (np.log(6e4) - np.log(2e4))
            train_m = (train_m - (-800)) / (800 - (-800))
            test_c = (test_c - 2e4) / (6e4 - 2e4)
            # test_c = (test_c - np.log(2e4)) / (np.log(6e4) - np.log(2e4))
            test_m = (test_m - (-800)) / (800 - (-800))

            # Expand color channel to be 3 for fit in incpetionResNetV2
            # print("train_c shape: ", train_c.shape)  # (?, 256, 256, 1)
            train_c = np.repeat(train_c, 3, -1)
            # print("train_c shape: ", train_c.shape)  # (?, 256, 256, 3)
            train_m = np.repeat(train_m, 3, -1)
            train_y = to_categorical(train_y)

            # print("test_c shape: ", test_c.shape)  # (?, 256, 256)
            test_c = np.repeat(test_c[..., np.newaxis], 3, -1)
            # print("test_c shape: ", test_c.shape)  # (?, 256, 256, 3)
            test_m = np.repeat(test_m[..., np.newaxis], 3, -1)
            test_y = to_categorical(test_y)
        return train_c, train_m, train_y, test_c, test_m, test_y
    else:
        return fold_c, fold_m, fold_y


# Visualize only the first 15 feature maps from each layer
def model_to_visualize_15(model, model_name, img_to_visualize, location):

    num_visible_layers = len(model.layers)-4

    for i in range(num_visible_layers):
        layer = model.get_layer(index=i)
        print("layer: ", layer.name)
        inputs = [K.learning_phase()] + model.inputs

        _convout1_f = K.function(inputs, [layer.output])

        def convout1_f(X):
            # The [0] is to disable the training phase flag
            return _convout1_f([0] + [X])

        convolutions = convout1_f(img_to_visualize)
        convolutions = np.squeeze(convolutions)

        if model_name == 'lstm_consecutive_images':
            convolutions = convolutions[0]  # Choose 1 image from the image sequence
        print('Shape of convolutions:', convolutions.shape)

        n_filter = convolutions.shape[2]
        n = int(np.ceil(np.sqrt(n_filter)))
        print("n: ", n)

        # # Visualization of each filter of the layer
        # fig = plt.figure(1)
        # for j in range(n_filter):
        #     # ax = fig.add_subplot(n, n, j + 1)                     # Plot to several columns
        #     ax = fig.add_subplot(n_filter, 1, j + 1)   # Plot to 1 column
        #     ax.imshow(convolutions[:,:,j], cmap='gray')
        #     ax.axis('off')
        # fig.suptitle("Activation: " + layer.name)
        # if location == 'local':
        #     plt.savefig(os.getcwd() + '/figs_cnnV2/activation_layer_' + str(i) + '_' + layer.name + '.png')
        # elif location == 'sharc':
        #     plt.savefig('/home/ms1yw/figs_cnnV2/activation_layer_' + str(i) + '_' + layer.name + '.png')

        fig = plt.figure(1)
        fig.set_size_inches(30, 30)
        # set up subplot grid
        gridspec.GridSpec(15, num_visible_layers)
        for j in range(15):
            plt.subplot2grid((15, num_visible_layers), (j, i))
            plt.locator_params(axis='x', nbins=5)
            plt.locator_params(axis='y', nbins=5)
            plt.imshow(convolutions[:, :, j], cmap='gray')
            plt.axis('off')
        plt.subplots_adjust(wspace=0.0, hspace=0.1)
    # if location == 'local':
    plt.savefig(os.getcwd() + '/figs_cnnV2/visualise_layers.png')
    # elif location == 'sharc':
    #     plt.savefig('/home/ms1yw/figs_cnnV2/visualise_layers.png')


# Visualize all feature maps from each layer
def model_to_visualize_all(model, model_name, img_to_visualize, which_image_to_print):
    os.makedirs(os.getcwd() + '/figs_cnnV2_' + str(which_image_to_print), exist_ok=True)

    num_visible_layers = len(model.layers)-4

    for i in range(num_visible_layers):
        layer = model.get_layer(index=i)
        print("layer: ", layer.name)
        inputs = [K.learning_phase()] + model.inputs

        _convout1_f = K.function(inputs, [layer.output])

        def convout1_f(X):
            # The [0] is to disable the training phase flag
            return _convout1_f([0] + [X])

        convolutions = convout1_f(img_to_visualize)
        convolutions = np.squeeze(convolutions)

        if model_name == 'lstm_consecutive_images':
            convolutions = convolutions[0]  # Choose 1 image from the image sequence
        print('Shape of convolutions:', convolutions.shape)

        n_filter = convolutions.shape[2]
        if n_filter == 32:
            row = 4
            col = 8
        elif n_filter == 64:
            row = 8
            col = 8
        elif n_filter == 128:
            row = 8
            col = 16
        elif n_filter == 256:
            row = 16
            col = 16

        # n = int(np.ceil(np.sqrt(n_filter)))
        # print("n: ", n)

        # Visualization of each filter of the layer
        fig = plt.figure(figsize=(col, row))
        for j in range(n_filter):
            ax = fig.add_subplot(row, col, j + 1)
            ax.imshow(convolutions[:, :, j], cmap='gray')
            # ax.imshow(convolutions[:, :, j])
            ax.axis('off')
        # fig.suptitle("Activation: " + layer.name)
        fig.subplots_adjust(wspace=0.05, hspace=0.05)

        # if location == 'local':
        plt.savefig(os.getcwd() + '/figs_cnnV2_' + str(which_image_to_print) + '/activation_layer_' + str(i) + '_' + layer.name + '.png')

        # This is another way of visualization, the effect is the same as the above way.
        # from keras.models import Model
        #
        # layer_outputs = [layer.output for layer in model.layers]
        # activation_model = Model(inputs=model.input, outputs=layer_outputs)
        # activations = activation_model.predict(img_to_visualize)
        #
        # def display_activation(activations, col_size, row_size, act_index):
        #     activation = activations[act_index]
        #     activation_index = 0
        #     fig, ax = plt.subplots(row_size, col_size, figsize=(row_size * 2.5, col_size * 1.5))
        #     for row in range(0, row_size):
        #         for col in range(0, col_size):
        #             ax[row][col].imshow(activation[0, :, :, activation_index], cmap='gray')
        #             activation_index += 1
        #     plt.show()
        #
        # display_activation(activations, 8, 4, 0)

def weights_to_visualize(layer, which_image_to_print):
    print("layer: ", layer.name)
    layer_weights = layer.get_weights()
    layer_weights = np.array(layer_weights)
    print("layer_weights shape: {}, data type: {}".format(layer_weights.shape, layer_weights[0].dtype))
    layer_weights = layer_weights[0]

    layer_shape = layer_weights.shape
    if layer_shape[2] == 1:  # if it is the first conv layer whose format is (x,x,1,x)

        if layer_shape[3] == 32:
            row = 4
            col = 8
        elif layer_shape[3] == 64:
            row = 8
            col = 8
        elif layer_shape[3] == 128:
            row = 8
            col = 16
        elif layer_shape[3] == 256:
            row = 16
            col = 16

        # n = int(np.ceil(np.sqrt(layer_shape[3])))
        # print("n: ", n)

        layer_weights = np.squeeze(layer_weights)
        print("layer_weights shape: ", layer_weights.shape)

        # Visualization of each filter of the layer
        fig = plt.figure(figsize=(col, row))
        for i in range(layer_shape[3]):
            ax = fig.add_subplot(row, col, i + 1)
            ax.imshow(layer_weights[:, :, i], cmap='gray', interpolation='none')
            ax.axis('off')
        # fig.suptitle("Weights: " + layer.name)
        fig.subplots_adjust(wspace=0.05, hspace=0.05)
        plt.savefig(os.getcwd() + '/figs_cnnV2_' + str(which_image_to_print) + '/weights_' + layer.name + '.png')
    else:
        n = layer_shape[2] * layer_shape[3]
        n = int(np.ceil(np.sqrt(n)))
        print("n: ", n)

        fig = plt.figure(figsize=(24, 15))
        for i in range(layer_shape[3]):
            for j in range(layer_shape[2]):
                ax = fig.add_subplot(n, n, i*32+j+1)
                single_layer_weights = layer_weights[:, :, j, i]
                # print("single_layer_weights shape1: ", single_layer_weights.shape)
                ax.imshow(single_layer_weights, cmap='gray', interpolation='none')
                ax.axis('off')
        fig.suptitle("Weights: " + layer.name)
        plt.savefig(os.getcwd() + '/figs_cnnV2_' + str(which_image_to_print) + '/weights_' + layer.name + '.png')

# def weights_to_visualize(layer, which_image_to_print):
#     print("layer: ", layer.name)
#     layer_weights = layer.get_weights()
#     layer_weights = np.array(layer_weights)
#     print("layer_weights shape: {}, data type: {}".format(layer_weights.shape, layer_weights[0].dtype))
#     layer_weights = layer_weights[0]
#
#     layer_shape = layer_weights.shape
#     if layer_shape[2] == 1:  # if it is the first conv layer whose format is (x,x,1,x)
#         n = int(np.ceil(np.sqrt(layer_shape[3])))
#         print("n: ", n)
#
#         layer_weights = np.squeeze(layer_weights)
#         print("layer_weights shape: ", layer_weights.shape)
#
#         # Visualization of each filter of the layer
#         fig = plt.figure(figsize=(24, 15))
#         for i in range(layer_shape[3]):
#             ax = fig.add_subplot(n, n, i + 1)
#             ax.imshow(layer_weights[:, :, i], cmap='gray', interpolation='none')
#             ax.axis('off')
#         fig.suptitle("Weights: " + layer.name)
#         plt.savefig(os.getcwd() + '/figs_cnnV2_' + str(which_image_to_print) + '/weights_' + layer.name + '.png')
#     else:
#         n = layer_shape[2] * layer_shape[3]
#         n = int(np.ceil(np.sqrt(n)))
#         print("n: ", n)
#
#         fig = plt.figure(figsize=(24, 15))
#         for i in range(layer_shape[3]):
#             for j in range(layer_shape[2]):
#                 ax = fig.add_subplot(n, n, i*32+j+1)
#                 single_layer_weights = layer_weights[:, :, j, i]
#                 # print("single_layer_weights shape1: ", single_layer_weights.shape)
#                 ax.imshow(single_layer_weights, cmap='gray', interpolation='none')
#                 ax.axis('off')
#         fig.suptitle("Weights: " + layer.name)
#         plt.savefig(os.getcwd() + '/figs_cnnV2_' + str(which_image_to_print) + '/weights_' + layer.name + '.png')