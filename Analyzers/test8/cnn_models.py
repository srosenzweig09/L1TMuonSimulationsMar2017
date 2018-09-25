import os
os.environ['KERAS_BACKEND'] = 'tensorflow'

import numpy as np

import tensorflow as tf

from keras import backend as K
from keras.models import Sequential, Model, clone_model, load_model, model_from_json
from keras.layers import Dense, Activation, Dropout, Input, Concatenate, Lambda, BatchNormalization
from keras.layers import Conv2D, MaxPooling2D, Flatten
from keras.callbacks import LearningRateScheduler, TerminateOnNaN, ModelCheckpoint
from keras.regularizers import Regularizer
from keras.constraints import Constraint
from keras import initializers, regularizers, optimizers, losses, metrics

from cnn_globals import (superstrip_size, n_zones, rows_per_zone, n_rows, n_columns, n_channels, n_classes)


def create_model(nvariables=None, lr=0.001, clipnorm=10., dropout=0.2, use_bn=True, use_dropout=False):
  input_shape = (n_rows, n_columns, n_channels)
  inputs = Input(shape=input_shape, dtype='float32')

  inputs_zone0 = Lambda(lambda x: x[:, 0*rows_per_zone:(0+1)*rows_per_zone])(inputs)
  inputs_zone1 = Lambda(lambda x: x[:, 1*rows_per_zone:(1+1)*rows_per_zone])(inputs)
  inputs_zone2 = Lambda(lambda x: x[:, 2*rows_per_zone:(2+1)*rows_per_zone])(inputs)
  inputs_zone3 = Lambda(lambda x: x[:, 3*rows_per_zone:(3+1)*rows_per_zone])(inputs)
  inputs_zone4 = Lambda(lambda x: x[:, 4*rows_per_zone:(4+1)*rows_per_zone])(inputs)
  inputs_zone5 = Lambda(lambda x: x[:, 5*rows_per_zone:(5+1)*rows_per_zone])(inputs)
  inputs_zone6 = Lambda(lambda x: x[:, 6*rows_per_zone:(6+1)*rows_per_zone])(inputs)

  def first_conv2d_layer_fn(inputs):
    x = Conv2D(
            12,
            (rows_per_zone,63),
            strides=(rows_per_zone,3),
            padding='same',
            kernel_initializer='he_normal',
            activation='relu')(inputs)
    return x

  x_zone0 = first_conv2d_layer_fn(inputs_zone0)
  x_zone1 = first_conv2d_layer_fn(inputs_zone1)
  x_zone2 = first_conv2d_layer_fn(inputs_zone2)
  x_zone3 = first_conv2d_layer_fn(inputs_zone3)
  x_zone4 = first_conv2d_layer_fn(inputs_zone4)
  x_zone5 = first_conv2d_layer_fn(inputs_zone5)
  x_zone6 = first_conv2d_layer_fn(inputs_zone6)

  x_all_zones = [x_zone0, x_zone1, x_zone2, x_zone3, x_zone4, x_zone5, x_zone6]
  x = Concatenate()(x_all_zones)

  x = Conv2D(
          48,
          (5,5),
          strides=(1,1),
          padding='same',
          kernel_initializer='he_normal',
          activation='relu')(x)
  x = MaxPooling2D(
          (1,5),
          strides=(1,3),
          padding='valid')(x)
  x = MaxPooling2D(
          (1,3),
          strides=(1,2),
          padding='valid')(x)
  x = Flatten()(x)

  if use_bn:
    x = Dense(40, use_bias=False, activation=None)(x)
    x = BatchNormalization(momentum=0.9, epsilon=1e-4)(x)
    x = Activation('relu')(x)
    if use_dropout:
      x = Dropout(dropout)(x)
    x = Dense(20, use_bias=False, activation=None)(x)
    x = BatchNormalization(momentum=0.9, epsilon=1e-4)(x)
    x = Activation('relu')(x)
    if use_dropout:
      x = Dropout(dropout)(x)
  else:
    x = Dense(40, activation='relu')(x)
    if use_dropout:
      x = Dropout(dropout)(x)
    x = Dense(20, activation='relu')(x)
    if use_dropout:
      x = Dropout(dropout)(x)

  outputs = Dense(n_classes, activation='softmax')(x)

  model = Model(inputs=inputs, outputs=outputs)

  # Set loss and optimizers
  adam = optimizers.Adam(lr=lr, clipnorm=clipnorm)

  def categorical_accuracy(y_true, y_pred):
    return K.cast(K.equal(K.argmax(y_true, axis=-1),
                          K.argmax(y_pred, axis=-1)),
                  K.floatx())

  def top_k_categorical_accuracy(y_true, y_pred, k=2):
    return K.mean(K.in_top_k(y_pred, K.argmax(y_true, axis=-1), k), axis=-1)

  model.compile(optimizer=adam, loss='categorical_crossentropy',
                metrics=['accuracy', categorical_accuracy, top_k_categorical_accuracy])
  model.summary()
  return model


# ______________________________________________________________________________
# Save/Load models

import h5py
import json

def save_my_model(model, name='model'):
  # Store model to file
  #model.summary()
  model.save(name + '.h5')
  model.save_weights(name + '_weights.h5')
  # Store model to json
  with open(name + '.json', 'w') as outfile:
    outfile.write(model.to_json())
  return

def load_my_model(name='model', weights_name='model_weights'):
  with open(name + '.json', 'r') as f:
    json_string = json.dumps(json.load(f))
    model = model_from_json(json_string)
  #model = load_model(name + '.h5')
  model.load_weights(weights_name + '.h5')
  return model
