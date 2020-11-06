
from keras.models import Sequential
from keras.layers import Activation, Dropout, Flatten, Conv2D, MaxPooling2D, Dense, TimeDistributed
from keras.layers.recurrent import LSTM
from keras.layers.normalization import BatchNormalization
from keras.layers.advanced_activations import LeakyReLU

from keras.applications.inception_resnet_v2 import InceptionResNetV2
from keras.applications.nasnet import NASNetLarge
from keras.applications.resnet50 import ResNet50
from keras.applications.xception import Xception
from keras.layers import Input, GlobalAveragePooling2D, Masking, Concatenate
from keras.models import Model


def mm_tseng(input_shape):  # multi-modality Tseng 2017
    # Build Xception over a custom input tensor
    input_tensor_1 = Input(shape=input_shape)
    input_tensor_2 = Input(shape=input_shape)

    base_model1 = InceptionResNetV2(input_tensor=input_tensor_1, weights='imagenet', include_top=False)
    for layer in base_model1.layers:
        layer.name = layer.name + str('_m')
    x1 = base_model1.output

    base_model2 = InceptionResNetV2(input_tensor=input_tensor_2, weights='imagenet', include_top=False)
    for layer in base_model2.layers:
        layer.name = layer.name + str('_c')
    x2 = base_model2.output

    x = Concatenate()([x1, x2])



    seq_len = input_shape[0]
    rows = input_shape[1]
    cols = input_shape[2]
    channels = input_shape[3]

    seq_input = Input(shape=input_shape)
    cnn_base = InceptionResNetV2(input_shape=(rows, cols, channels), weights='imagenet', include_top=False)
    # cnn_base.trainable = False

    # Add 3D convolution
    cnn_base

    # Add a global spatial average pooling layer
    # cnn_out = GlobalAveragePooling2D()(cnn_base.output)
    cnn_out = Flatten()(cnn_base.output)
    cnn = Model(inputs=cnn_base.input, outputs=cnn_out)

    x = TimeDistributed(cnn)(seq_input)
    x = LSTM(1024)(x)
    # x = LSTM(512)(x) ???
    x = Dense(1024, activation='relu')(x)
    x = Dense(512, activation='relu')(x)
    predictions = Dense(3, activation='softmax')(x)
    model = Model([seq_input], predictions)
    # model = Model(inputs=cnn_base.input, outputs=predictions)

    return model


def inceptionResNetV2_plus_LSTM(input_shape):
    seq_len = input_shape[0]
    rows = input_shape[1]
    cols = input_shape[2]
    channels = input_shape[3]

    seq_input = Input(shape=input_shape)
    cnn_base = InceptionResNetV2(input_shape=(rows, cols, channels), weights='imagenet', include_top=False)
    # cnn_base.trainable = False

    # add a global spatial average pooling layer
    # cnn_out = GlobalAveragePooling2D()(cnn_base.output)
    cnn_out = Flatten()(cnn_base.output)
    cnn = Model(inputs=cnn_base.input, outputs=cnn_out)

    x = TimeDistributed(cnn)(seq_input)
    x = LSTM(1024)(x)
    # x = LSTM(512)(x) ???
    x = Dense(1024, activation='relu')(x)
    x = Dense(512, activation='relu')(x)
    predictions = Dense(3, activation='softmax')(x)
    model = Model([seq_input], predictions)
    # model = Model(inputs=cnn_base.input, outputs=predictions)

    return model


def lstm(input_shape):
    # exp5
    ks = 3
    model = Sequential()  # Change receptive field the first Conv2D from 7 to 3
    model.add(TimeDistributed(Conv2D(32, (ks, ks), strides=(2, 2),
                                     activation='relu', padding='same'), input_shape=input_shape))
    model.add(TimeDistributed(Conv2D(32, (ks, ks),
                                     kernel_initializer="he_normal", activation='relu')))
    model.add(TimeDistributed(MaxPooling2D((2, 2), strides=(2, 2))))
    model.add(TimeDistributed(Conv2D(64, (ks, ks),
                                     padding='same', activation='relu')))
    model.add(TimeDistributed(Conv2D(64, (ks, ks),
                                     padding='same', activation='relu')))
    model.add(TimeDistributed(MaxPooling2D((2, 2), strides=(2, 2))))
    model.add(TimeDistributed(Conv2D(128, (ks, ks),
                                     padding='same', activation='relu')))
    model.add(TimeDistributed(Conv2D(128, (ks, ks),
                                     padding='same', activation='relu')))
    model.add(TimeDistributed(MaxPooling2D((2, 2), strides=(2, 2))))
    model.add(TimeDistributed(Conv2D(256, (ks, ks),
                                     padding='same', activation='relu')))
    model.add(TimeDistributed(Conv2D(256, (ks, ks),
                                     padding='same', activation='relu')))
    model.add(TimeDistributed(MaxPooling2D((2, 2), strides=(2, 2))))
    model.add(TimeDistributed(Conv2D(512, (ks, ks),
                                     padding='same', activation='relu')))
    model.add(TimeDistributed(Conv2D(512, (ks, ks),
                                     padding='same', activation='relu')))
    model.add(TimeDistributed(MaxPooling2D((2, 2), strides=(2, 2))))
    model.add(TimeDistributed(Flatten()))
    # model.add(LSTM(512, activation='relu', return_sequences=False))
    model.add(LSTM(512))
    model.add(Dense(512, activation='relu'))
    model.add(Dense(256, activation='relu'))
    model.add(Dense(3, activation='softmax'))

    return model


def inceptionResNetV2_doubleV2(input_shape):
    # Build Xception over a custom input tensor
    input_tensor_1 = Input(shape=input_shape)
    input_tensor_2 = Input(shape=input_shape)

    base_model1 = InceptionResNetV2(input_tensor=input_tensor_1, weights='imagenet', include_top=False)
    for layer in base_model1.layers:
        layer.name = layer.name + str('_m')
    x1 = base_model1.output

    base_model2 = InceptionResNetV2(input_tensor=input_tensor_2, weights='imagenet', include_top=False)
    for layer in base_model2.layers:
        layer.name = layer.name + str('_c')
    x2 = base_model2.output

    x = Concatenate()([x1, x2])

    # add a global spatial average pooling layer
    x = Dense(1024, activation='relu')(x)
    x = GlobalAveragePooling2D()(x)
    x = Dropout(0.2)(x)
    x = Dense(512, activation='relu')(x)
    x = Dense(128, activation='relu')(x)
    predictions = Dense(3, activation='softmax')(x)

    # this is the model we will train
    # model = Model(inputs=base_model.input, outputs=predictions)
    model = Model(inputs=[input_tensor_1, input_tensor_2], outputs=predictions)

    return model


def inceptionResNetV2_doubleV1(input_shape):
    # Build Xception over a custom input tensor
    input_tensor_1 = Input(shape=input_shape)
    input_tensor_2 = Input(shape=input_shape)

    base_model1 = InceptionResNetV2(input_tensor=input_tensor_1, weights='imagenet', include_top=False)
    for layer in base_model1.layers:
        layer.name = layer.name + str('_m')
    # add a global spatial average pooling layer
    x1 = base_model1.output
    x1 = Dense(1024, activation='relu')(x1)
    x1 = GlobalAveragePooling2D()(x1)
    x1 = Dropout(0.2)(x1)
    x1 = Dense(512, activation='relu')(x1)
    x1 = Dense(128, activation='relu')(x1)

    base_model2 = InceptionResNetV2(input_tensor=input_tensor_2, weights='imagenet', include_top=False)
    for layer in base_model2.layers:
        layer.name = layer.name + str('_c')
    # add a global spatial average pooling layer
    x2 = base_model2.output
    x2 = Dense(1024, activation='relu')(x2)
    x2 = GlobalAveragePooling2D()(x2)
    x2 = Dropout(0.2)(x2)
    x2 = Dense(512, activation='relu')(x2)
    x2 = Dense(128, activation='relu')(x2)

    x = Concatenate()([x1, x2])
    predictions = Dense(3, activation='softmax')(x)

    # this is the model we will train
    # model = Model(inputs=base_model.input, outputs=predictions)
    model = Model(inputs=[input_tensor_1, input_tensor_2], outputs=predictions)

    return model

def xCeption(input_shape):
    # exp4 & exp6 & exp8
    # Build Xception over a custom input tensor
    input_tensor = Input(shape=input_shape)
    base_model = Xception(input_tensor=input_tensor, weights='imagenet', include_top=False)
    base_model.pre
    # base_model.save(os.getcwd() + '/models/InceptionResNetV2_baseModel.h5')

    # # let's visualize layer names and layer indices
    # for i, layer in enumerate(base_model.layers):
    #     print(i, layer.name)

    # # inter_model = Model(inputs=base_model.input, outputs=base_model.get_layer('block1_conv1').output)
    # inter_model = Model(inputs=base_model.input, outputs=base_model.get_layer('conv2d_1').output)
    # feature = inter_model.predict(test_X)
    # print("feature shape from keras.models import load_model11: ", feature.shape)
    # feature_to_visualize(feature, 11)
    # feature = base_model.predict(test_X)
    # print("feature shape12: ", feature.shape)
    # feature_to_visualize(feature, 12)

    # add a global spatial average pooling layer
    x = base_model.output
    x = Dense(1024, activation='relu')(x)
    x = GlobalAveragePooling2D()(x)
    x = Dropout(0.2)(x)
    x = Dense(512, activation='relu')(x)
    x = Dense(128, activation='relu')(x)
    predictions = Dense(3, activation='softmax')(x)

    # this is the model we will train
    model = Model(inputs=base_model.input, outputs=predictions)

    return model


def resNet50(input_shape):
    # exp4 & exp6 & exp8
    # Build Xception over a custom input tensor
    input_tensor = Input(shape=input_shape)
    base_model = ResNet50(input_tensor=input_tensor, weights='imagenet', include_top=False)
    # base_model.save(os.getcwd() + '/models/InceptionResNetV2_baseModel.h5')

    # # let's visualize layer names and layer indices
    # for i, layer in enumerate(base_model.layers):
    #     print(i, layer.name)
    # base_model.summary()

    # # inter_model = Model(inputs=base_model.input, outputs=base_model.get_layer('block1_conv1').output)
    # inter_model = Model(inputs=base_model.input, outputs=base_model.get_layer('conv2d_1').output)
    # feature = inter_model.predict(test_X)
    # print("feature shape from keras.models import load_model11: ", feature.shape)
    # feature_to_visualize(feature, 11)
    # feature = base_model.predict(test_X)
    # print("feature shape12: ", feature.shape)
    # feature_to_visualize(feature, 12)

    # add a global spatial average pooling layer
    x = base_model.output
    # x = Dense(1024, activation='relu')(x)
    x = GlobalAveragePooling2D()(x)
    # x = Dropout(0.2)(x)
    x = Dense(1024, activation='relu')(x)
    x = Dense(512, activation='relu')(x)
    predictions = Dense(3, activation='softmax')(x)

    # this is the model we will train
    model = Model(inputs=base_model.input, outputs=predictions)

    return model


def inceptionResNetV2_2class(input_shape):
    # Build Xception over a custom input tensor
    input_tensor = Input(shape=input_shape)
    base_model = InceptionResNetV2(input_tensor=input_tensor, weights='imagenet', include_top=False)
    # base_model.save(os.getcwd() + '/models/InceptionResNetV2_baseModel.h5')

    # # let's visualize layer names and layer indices
    # for i, layer in enumerate(base_model.layers):
    #     print(i, layer.name)

    # # inter_model = Model(inputs=base_model.input, outputs=base_model.get_layer('block1_conv1').output)
    # inter_model = Model(inputs=base_model.input, outputs=base_model.get_layer('conv2d_1').output)
    # feature = inter_model.predict(test_X)
    # print("feature shape from keras.models import load_model11: ", feature.shape)
    # feature_to_visualize(feature, 11)
    # feature = base_model.predict(test_X)
    # print("feature shape12: ", feature.shape)
    # feature_to_visualize(feature, 12)

    # add a global spatial average pooling layer
    x = base_model.output
    x = Dense(1024, activation='relu')(x)
    x = GlobalAveragePooling2D()(x)
    x = Dropout(0.2)(x)
    x = Dense(512, activation='relu')(x)
    x = Dense(128, activation='relu')(x)
    predictions = Dense(2, activation='softmax')(x)

    # this is the model we will train
    model = Model(inputs=base_model.input, outputs=predictions)

    return model


def inceptionResNetV2(input_shape):
    # exp4 & exp6 & exp8
    # Build Xception over a custom input tensor
    input_tensor = Input(shape=input_shape)
    base_model = InceptionResNetV2(input_tensor=input_tensor, weights='imagenet', include_top=False)
    # base_model.save(os.getcwd() + '/models/InceptionResNetV2_baseModel.h5')

    # # let's visualize layer names and layer indices
    # for i, layer in enumerate(base_model.layers):
    #     print(i, layer.name)

    # # inter_model = Model(inputs=base_model.input, outputs=base_model.get_layer('block1_conv1').output)
    # inter_model = Model(inputs=base_model.input, outputs=base_model.get_layer('conv2d_1').output)
    # feature = inter_model.predict(test_X)
    # print("feature shape from keras.models import load_model11: ", feature.shape)
    # feature_to_visualize(feature, 11)
    # feature = base_model.predict(test_X)
    # print("feature shape12: ", feature.shape)
    # feature_to_visualize(feature, 12)

    # add a global spatial average pooling layer
    x = base_model.output
    x = Dense(1024, activation='relu')(x)
    x = GlobalAveragePooling2D()(x)
    x = Dropout(0.2)(x)
    x = Dense(512, activation='relu')(x)
    x = Dense(128, activation='relu')(x)
    predictions = Dense(3, activation='softmax')(x)

    # this is the model we will train
    model = Model(inputs=base_model.input, outputs=predictions)

    return model


def nasNetLarge(input_shape):
    # Build Xception over a custom input tensor
    input_tensor = Input(shape=input_shape)
    base_model = NASNetLarge(input_tensor=input_tensor, include_top=False, weights='imagenet')

    # add a global spatial average pooling layer
    x = base_model.output
    x = Dense(1024, activation='relu')(x)
    x = GlobalAveragePooling2D()(x)
    x = Dropout(0.2)(x)
    x = Dense(512, activation='relu')(x)
    x = Dense(128, activation='relu')(x)
    predictions = Dense(3, activation='softmax')(x)

    # this is the model we will train
    model = Model(inputs=base_model.input, outputs=predictions)

    return model


def cnn_model_v1(input_shape):
    ######################## BN makes bad performance: increase loss & decrease accuracy from training begins!!!
    # exp1 & exp2
    model = Sequential()
    model.add(Conv2D(32, kernel_size=(5, 5), input_shape=input_shape, padding='same'))
    # model.add(BatchNormalization(momentum=0.7))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2), padding='same'))
    model.add(Conv2D(64, kernel_size=(5, 5), padding='same'))
    # model.add(BatchNormalization(momentum=0.7))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2), padding='same'))
    model.add(Conv2D(128, kernel_size=(5, 5), padding='same'))
    # model.add(BatchNormalization(momentum=0.7))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2), padding='same'))
    model.add(Flatten())
    model.add(Dense(512))
    model.add(Activation('relu'))
    model.add(Dense(3))
    model.add(Activation('softmax'))

    # # exp3
    # model = Sequential()
    # model.add(Conv2D(64, kernel_size=(5, 5), input_shape=input_shape, padding='same'))
    # # model.add(BatchNormalization(momentum=0.7))
    # model.add(Activation('relu'))
    # model.add(MaxPooling2D(pool_size=(2, 2), padding='same'))
    # model.add(Conv2D(128, kernel_size=(5, 5), padding='same'))
    # # model.add(BatchNormalization(momentum=0.7))
    # model.add(Activation('relu'))
    # model.add(MaxPooling2D(pool_size=(2, 2), padding='same'))
    # model.add(Conv2D(256, kernel_size=(5, 5), padding='same'))
    # # model.add(BatchNormalization(momentum=0.7))
    # model.add(Activation('relu'))
    # model.add(MaxPooling2D(pool_size=(2, 2), padding='same'))
    # model.add(Flatten())
    # model.add(Dense(256))
    # model.add(Activation('relu'))
    # model.add(Dense(3))
    # model.add(Activation('softmax'))

    # # exp7
    # model = Sequential()
    # model.add(Conv2D(32, kernel_size=(3, 3), input_shape=input_shape, padding='same'))
    # # model.add(BatchNormalization(momentum=0.7))
    # model.add(Activation('relu'))
    # model.add(MaxPooling2D(pool_size=(2, 2), padding='same'))
    # model.add(Conv2D(64, kernel_size=(3, 3), padding='same'))
    # # model.add(BatchNormalization(momentum=0.7))
    # model.add(Activation('relu'))
    # model.add(MaxPooling2D(pool_size=(2, 2), padding='same'))
    # model.add(Conv2D(128, kernel_size=(3, 3), padding='same'))
    # # model.add(BatchNormalization(momentum=0.7))
    # model.add(Activation('relu'))
    # model.add(MaxPooling2D(pool_size=(2, 2), padding='same'))
    # model.add(Flatten())
    # model.add(Dense(512))
    # model.add(Activation('relu'))
    # model.add(Dense(3))
    # model.add(Activation('softmax'))

    return model
