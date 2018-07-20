import os
os.environ['KERAS_BACKEND'] = 'tensorflow'

import tensorflow as tf

from keras import backend as K
from keras.models import Sequential, Model, load_model, model_from_json
from keras.layers import Dense, Activation, Dropout, Input, BatchNormalization
from keras import initializers, regularizers, optimizers, losses

import h5py
import json

from nn_logging import getLogger
logger = getLogger()


# ______________________________________________________________________________
# New leaky relu
def NewLeakyReLU(x, alpha=0., max_value=None):
  return K.relu(x, alpha=alpha, max_value=max_value)

# ______________________________________________________________________________
# New tanh
def NewTanh(x):
  return K.tanh(x)
  #return 1.7159 * K.tanh(x * 2./3.)
  #return K.clip(x, -1., 1.)

# ______________________________________________________________________________
# Huber loss

def huber_loss(y_true, y_pred, delta=1.345):
  x = K.abs(y_true - y_pred)
  squared_loss = 0.5*K.square(x)
  absolute_loss = delta * (x - 0.5*delta)
  #xx = K.switch(x < delta, squared_loss, absolute_loss)
  xx = tf.where(x < delta, squared_loss, absolute_loss)  # needed for tensorflow
  return K.mean(xx, axis=-1)

def masked_huber_loss(y_true, y_pred, delta=1.345, mask_value=100.):
  x = K.abs(y_true - y_pred)
  squared_loss = 0.5*K.square(x)
  absolute_loss = delta * (x - 0.5*delta)
  #xx = K.switch(x < delta, squared_loss, absolute_loss)
  xx = tf.where(x < delta, squared_loss, absolute_loss)  # needed for tensorflow

  mask = K.not_equal(y_true, mask_value)
  mask = K.cast(mask, K.floatx())
  xx *= mask
  xx /= K.mean(mask)
  return K.mean(xx, axis=-1)

# ______________________________________________________________________________
# Binary crossentropy

# See: https://github.com/tensorflow/tensorflow/blob/master/tensorflow/python/keras/losses.py
#      https://github.com/tensorflow/tensorflow/blob/master/tensorflow/python/keras/backend.py
#def binary_crossentropy(y_true, y_pred):
#  return K.mean(K.binary_crossentropy(y_true, y_pred), axis=-1)

def masked_binary_crossentropy(y_true, y_pred, mask_value=100.):
  xx = K.binary_crossentropy(y_true, y_pred)

  mask = K.not_equal(y_true, mask_value)
  mask = K.cast(mask, K.floatx())
  xx *= mask
  xx /= K.mean(mask)
  return K.mean(xx, axis=-1)

# ______________________________________________________________________________
# Learning rate decay by epoch number

def lr_schedule(epoch):
  if (epoch % 10) == 0:
    lr = K.get_value(model.optimizer.lr)
    K.set_value(model.optimizer.lr, lr*0.95)
    print("lr changed to {}".format(lr*0.95))
  return K.get_value(model.optimizer.lr)

from keras.callbacks import LearningRateScheduler
lr_decay = LearningRateScheduler(lr_schedule)

# ______________________________________________________________________________
# Terminate training on NaN loss

from keras.callbacks import TerminateOnNaN
terminate_on_nan = TerminateOnNaN()


# ______________________________________________________________________________
def create_model(nvariables, lr=0.001, discr_loss_weight=1.0):
  inputs = Input(shape=(nvariables,), dtype='float32')

  x = Dense(64, activation='tanh', kernel_initializer='glorot_uniform', kernel_regularizer=regularizers.l2(0.0000))(inputs)
  #x = Dropout(0.2)(x)
  x = Dense(32, activation='tanh', kernel_initializer='glorot_uniform', kernel_regularizer=regularizers.l2(0.0000))(x)
  #x = Dropout(0.2)(x)
  x = Dense(16, activation='tanh', kernel_initializer='glorot_uniform', kernel_regularizer=regularizers.l2(0.0000))(x)
  #x = Dropout(0.2)(x)

  regr = Dense(1, activation='linear', kernel_initializer='glorot_uniform', name='regr')(x)
  discr = Dense(1, activation='sigmoid', kernel_initializer='glorot_uniform', name='discr')(x)

  # This creates a model that includes
  # the Input layer, three Dense layers and the Output layer
  model = Model(inputs=inputs, outputs=[regr, discr])

  # Set loss and optimizers
  #binary_crossentropy = losses.binary_crossentropy
  #mean_squared_error = losses.mean_squared_error

  adam = optimizers.Adam(lr=lr)
  model.compile(optimizer=adam,
    loss={'regr': masked_huber_loss, 'discr': masked_binary_crossentropy},
    loss_weights={'regr': 1.0, 'discr': discr_loss_weight},
    #metrics={'regr': ['acc', 'mse', 'mae'], 'discr': ['acc',]}
    )
  return model


# ______________________________________________________________________________
def create_model_sequential(nvariables, lr=0.001):
  model = Sequential()
  model.add(Dense(64, input_dim=nvariables, activation='tanh', kernel_initializer='glorot_uniform', kernel_regularizer=regularizers.l2(0.0000)))
  #model.add(Dropout(0.2))
  model.add(Dense(32, activation='tanh', kernel_initializer='glorot_uniform', kernel_regularizer=regularizers.l2(0.0000)))
  #model.add(Dropout(0.2))
  model.add(Dense(16, activation='tanh', kernel_initializer='glorot_uniform', kernel_regularizer=regularizers.l2(0.0000)))
  #model.add(Dropout(0.2))
  model.add(Dense(1, activation='linear', kernel_initializer='glorot_uniform'))

  adam = optimizers.Adam(lr=lr)
  model.compile(loss=huber_loss, optimizer=adam, metrics=['acc'])
  return model


# ______________________________________________________________________________
def save_my_model(model, name='model'):
  # Store model to file
  model.summary()
  model.save(name + '.h5')
  model.save_weights(name + '_weights.h5')
  # Store model to json
  with open(name + '.json', 'w') as outfile:
    outfile.write(model.to_json())
  logger.info('Saved model as {0}.h5, {0}.json and {0}_weights.h5'.format(name))
  return

def load_my_model(name='model'):
  with open(name + '.json', 'r') as f:
    json_string = json.dumps(json.load(f))
    model = model_from_json(json_string)
  #model = load_model(name + '.h5')
  model.load_weights(name + '_weights.h5')
  logger.info('Loaded model from {0} and weights from {1}'.format(name + '.json', name + '_weights.h5'))
  return model
