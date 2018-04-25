import numpy as np
np.random.seed(2023)
import random
random.seed(2023)
print('[INFO] Using numpy {0}'.format(np.__version__))

import os
import sys
import time
os.environ['KERAS_BACKEND'] = 'tensorflow'
old_stdout = sys.stdout

import keras
import keras.backend as K
from keras.models import Sequential, Model
from keras.layers import Dense, Activation, Dropout, Input, BatchNormalization
from keras import initializers, regularizers, optimizers, losses
K.set_epsilon(1e-08)
print('[INFO] Using keras {0}'.format(keras.__version__))

import tensorflow as tf
print('[INFO] Using tensorflow {0}'.format(tf.__version__))

import sklearn
#from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
print('[INFO] Using sklearn {0}'.format(sklearn.__version__))

#import matplotlib.pyplot as plt
#from matplotlib import colors
#%matplotlib inline

# Misc
from encoder import Encoder, NewLeakyReLU, huber_loss, masked_huber_loss, lr_decay
from utils import slice_arrays, merge_arrays, make_batches


# ______________________________________________________________________________
# Globals
nlayers = 12  # 5 (CSC) + 4 (RPC) + 3 (GEM)

nvariables = 68

do_skim = True

add_noise = True


# ______________________________________________________________________________
# Functions

def muon_data():
  try:
    print('[INFO] Loading muon data ...')
    infile = '../test2/histos_tba.12.npz'
    loaded = np.load(infile)
    the_variables = loaded['variables']
    the_parameters = loaded['parameters']
    print('[INFO] Loaded the variables with shape {0}'.format(the_variables.shape))
    print('[INFO] Loaded the parameters with shape {0}'.format(the_parameters.shape))
  except:
    print('[ERROR] Failed to load data from file: {0}'.format(infile))

  if do_skim:
    the_variables = the_variables[:1000]
    the_parameters = the_parameters[:1000]

  encoder = Encoder(the_variables, the_parameters, adjust_scale=2)
  x, y, w, x_mask = encoder.get_x(), encoder.get_y(), encoder.get_w(), encoder.get_x_mask()
  assert np.isfinite(x).all()

  # Split dataset in training and testing
  x_train, x_test, y_train, y_test, w_train, w_test, x_mask_train, x_mask_test = train_test_split(x, y, w, x_mask, test_size=0.4)
  return x_train, x_test, y_train, y_test, w_train, w_test, x_mask_train, x_mask_test

def pileup_data():
  try:
    print('[INFO] Loading pileup data ...')
    infile = '../test2/histos_tbd.test.npz'
    loaded = np.load(infile)
    the_variables = loaded['variables']
    the_parameters = np.zeros((the_variables.shape[0], 3), dtype=np.float32)
    the_auxiliaries = loaded['aux']
    print('[INFO] Loaded the variables with shape {0}'.format(the_variables.shape))
    print('[INFO] Loaded the auxiliary info with shape {0}'.format(the_auxiliaries.shape))
  except:
    print('[ERROR] Failed to load data from file: {0}'.format(infile))

  if do_skim:
    the_variables = the_variables[:1000]
    the_parameters = the_parameters[:1000]
    the_auxiliaries = the_auxiliaries[:1000]

  encoder = Encoder(the_variables, the_parameters, adjust_scale=2)
  x = encoder.get_x()
  aux = the_auxiliaries  # jobid, ievt, highest_part_pt, highest_track_pt
  assert np.isfinite(x).all()
  return x, aux

def create_model():
  # This returns a tensor
  inputs = Input(shape=(nvariables,), dtype='float32')

  # a layer instance is callable on a tensor, and returns a tensor
  x = Dense(64, activation='relu', kernel_initializer='glorot_uniform', kernel_regularizer=regularizers.l2(0.0000))(inputs)
  #x = Dropout(0.2)(x)
  x = Dense(32, activation='relu', kernel_initializer='glorot_uniform', kernel_regularizer=regularizers.l2(0.0000))(x)
  #x = Dropout(0.2)(x)
  x = Dense(16, activation='relu', kernel_initializer='glorot_uniform', kernel_regularizer=regularizers.l2(0.0000))(x)
  #x = Dropout(0.2)(x)
  regr = Dense(1, activation='linear', kernel_initializer='glorot_uniform', name='regr')(x)
  discr = Dense(1, activation='sigmoid', kernel_initializer='glorot_uniform', name='discr')(x)

  # This creates a model that includes
  # the Input layer, three Dense layers and the Output layer
  model = Model(inputs=inputs, outputs=[regr, discr])

  # Set loss and optimizers
  binary_crossentropy = losses.binary_crossentropy
  adam = optimizers.Adam(lr=0.001)
  #adam = optimizers.Adam(lr=0.001, amsgrad=True)

  # Compile
  model.compile(optimizer=adam,
    loss={'regr': masked_huber_loss, 'discr': binary_crossentropy},
    loss_weights={'regr': 1.0, 'discr': 0.2},
    metrics={'regr': ['acc', 'mse', 'mae'], 'discr': ['acc',]})
  return model

def save_model(model):
  # Store model to file
  model.summary()
  model.save('model.h5')
  model.save_weights('model_weights.h5')

  # Store model to json
  import json
  with open('model.json', 'w') as outfile:
    json.dump(model.to_json(), outfile)
  return


# Training routine
def train(model, x, y, x_adv, aux_adv, batch_size=None, epochs=1, verbose=1, callbacks=None,
          validation_split=0., shuffle=True, class_weight=None, sample_weight=None):

  # Validate user data.
  x, y, sample_weights = model._standardize_user_data(
    x, y,
    sample_weight=sample_weight,
    class_weight=class_weight,
    batch_size=batch_size)
  ins = x + y + sample_weights

  # Prepare validation data.
  do_validation = False
  if validation_split and 0. < validation_split < 1.:
    do_validation = True
    if hasattr(x[0], 'shape'):
      split_at = int(x[0].shape[0] * (1. - validation_split))
    else:
      split_at = int(len(x[0]) * (1. - validation_split))
    x, val_x = (slice_arrays(x, 0, split_at),
                slice_arrays(x, split_at))
    y, val_y = (slice_arrays(y, 0, split_at),
                slice_arrays(y, split_at))
    sample_weights, val_sample_weights = (
      slice_arrays(sample_weights, 0, split_at),
      slice_arrays(sample_weights, split_at))
    val_ins = val_x + val_y + val_sample_weights
  else:
    val_ins = []


  # ____________________________________________________________________________
  # Fit

  num_train_samples = x[0].shape[0]
  index_array = np.arange(num_train_samples)

  #sess = K.get_session()
  #tf_x = K.placeholder(shape=(None, x[0].shape[1]))
  #tf_y = K.placeholder(shape=(None, y[0].shape[1]))
  #sess.run(tf.initialize_all_variables())

  # Loop over epochs
  for epoch in xrange(epochs):
    print('Epoch {0}'.format(epoch))

    if shuffle:
      np.random.shuffle(index_array)

    batches = make_batches(num_train_samples, batch_size)

    # Loop over batches
    for batch_index, (batch_start, batch_end) in enumerate(batches):
      print('.. batch {0}'.format(batch_index))

      batch_ids = index_array[batch_start:batch_end]
      ins_batch = slice_arrays(ins, batch_ids)
      assert isinstance(ins_batch, list) and len(ins_batch) == 1 + 2 + 2

      # Add noise
      if add_noise:
        noise = x_adv[np.random.randint(0, x_adv.shape[0], ins_batch[0].shape[0])]
        noise_reg = np.zeros_like(ins_batch[1]) + 100.  # mask_value is set to 100
        noise_discr = np.zeros_like(ins_batch[2])
        noise_reg_w = np.ones_like(ins_batch[3])
        noise_discr_w = np.ones_like(ins_batch[3])
        ins_noise = [noise, noise_reg, noise_discr, noise_reg_w, noise_discr_w]
        ins_batch = merge_arrays(ins_batch, ins_noise)

      model._make_train_function()
      f = model.train_function
      outs = f(ins_batch)



# ______________________________________________________________________________
# Main

if __name__ == '__main__':
  assert keras.backend.backend() == 'tensorflow'

  x_train, x_test, y_train, y_test, w_train, w_test, x_mask_train, x_mask_test = muon_data()

  # Add output nodes
  y_train = [y_train, np.ones_like(y_train)]
  y_test = [y_test, np.ones_like(y_test)]

  x_adv, aux_adv = pileup_data()

  model = create_model()

  #history = train(model, x_train, y_train, epochs=5, validation_split=0.1, batch_size=256, verbose=1)
  train(model, x_train, y_train, x_adv, aux_adv, epochs=5, validation_split=0.1, batch_size=256, verbose=1)
  save_model(model)

