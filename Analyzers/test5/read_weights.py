model_file = 'model.h5'
model_weights_file = 'model_weights.h5'

import numpy as np
import h5py

def read_h5py():
  f = h5py.File(model_weights_file)
  print f.keys()
  
  keys = [u'dense_1', u'dense_2', u'dense_3', u'regr', u'discr']
  for k in keys:
    try:
      w = f[k][k]['kernel:0'].value
      b = f[k][k]['bias:0'].value
      print k, w, b, np.min(np.abs(w)), np.max(np.abs(w))
    except:
      pass


import os
import sys
import time
os.environ['KERAS_BACKEND'] = 'tensorflow'
from keras import backend as K
from keras.models import load_model
import tensorflow as tf

def NewTanh(x):
  return K.tanh(x)

def masked_huber_loss(y_true, y_pred, delta=1.345):
  mask_value = 100.
  mask_alpha = 0.02
  mask = K.equal(y_true, mask_value)
  
  #x = K.abs(y_true - y_pred)
  x = tf.where(mask, mask_alpha * K.abs(0.5 * 5 - K.abs(y_pred)), K.abs(y_true - y_pred))  #FIXME
  squared_loss = 0.5*K.square(x)
  absolute_loss = delta * (x - 0.5*delta)
  #xx = K.switch(x < delta, squared_loss, absolute_loss)
  xx = tf.where(x < delta, squared_loss, absolute_loss)  # needed for tensorflow
  return K.mean(xx, axis=-1)

def masked_binary_crossentropy(y_true, y_pred, from_logits=False):
  target, output = y_true, y_pred

  # transform back to logits
  if not from_logits:
    output = K.clip(output, K.epsilon(), 1 - K.epsilon())
    output = K.log(output / (1 - output))
  
  xx =  tf.nn.sigmoid_cross_entropy_with_logits(labels=target, logits=output)
  #xx =  tf.nn.weighted_cross_entropy_with_logits(targets=target, logits=output, pos_weight=0.5)  # pos_weight < 1 decreases the false positive count

  mask_value = 100.
  mask = K.equal(y_true, mask_value)
  mask = 1 - K.cast(mask, K.floatx())
  return K.sum(xx * mask, axis=-1)

def read_keras():
  custom_objects = {'masked_huber_loss': masked_huber_loss, 'masked_binary_crossentropy': masked_binary_crossentropy, 'NewTanh': NewTanh}
  loaded_model = load_model(model_file, custom_objects=custom_objects)
  loaded_model.load_weights(model_weights_file)
  for ilayer, layer in enumerate(loaded_model.layers):
    try:
      k = layer.name
      w =  layer.get_weights()
      w, b = w[0], w[1]
      print k, w, b, np.min(np.abs(w)), np.max(np.abs(w))
    except:
      pass
  

if __name__ == '__main__':

  #read_h5py()

  read_keras()
