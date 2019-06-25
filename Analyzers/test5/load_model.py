#!/usr/bin/env python

import os
import sys
os.environ['KERAS_BACKEND'] = 'tensorflow'

def usage():
  print('usage: python {0} FILE'.format(sys.argv[0]))
  print('')
  print('arguments:')
  print('  FILE    a model JSON file, e.g. \'model.json\'')


# ______________________________________________________________________________
if __name__ == "__main__":
  if len(sys.argv) < 2:
    usage()
    sys.exit(1)

  model_file = sys.argv[1]
  model_weights_file = model_file.replace('model', 'model_weights').replace('.json', '.h5')
  infile_muon = model_file.replace('model', 'histos_tba').replace('.json', '.npz')

  # Load data
  import numpy as np
  np.random.seed(2023)
  with np.load(infile_muon) as loaded:
    the_variables = loaded['variables']
    the_parameters = loaded['parameters']
  print('Loaded the variables with shape {0} and the parameters with shape {1}'.format(the_variables.shape, the_parameters.shape))

  # Data preprocessing
  from nn_encode import nlayers, nvariables, nvariables_input, nparameters_input, Encoder
  reg_pt_scale = 100.
  reg_dxy_scale = 0.4

  encoder = Encoder(the_variables, the_parameters, reg_pt_scale=reg_pt_scale, reg_dxy_scale=reg_dxy_scale)
  x, y = encoder.get_x(), encoder.get_y()

  # Split dataset in training and testing
  from sklearn.model_selection import train_test_split
  test_size = 0.3

  x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=test_size)
  print('Loaded # of training and testing events: {0}'.format((x_train.shape[0], x_test.shape[0])))

  # Load model
  from nn_models import load_my_model, update_keras_custom_objects
  update_keras_custom_objects()
  loaded_model = load_my_model(name=model_file, weights_name=model_weights_file)
  print('Loaded model.')

  # Run
  print('Applying the model. This will take a while ...')
  y_train_true = y_train[:,np.newaxis].copy()
  y_train_pred, y_train_discr = loaded_model.predict(x_train)
  #y_train_true /= reg_pt_scale
  #y_train_pred /= reg_pt_scale

  # Print
  N = 10
  print np.array2string(x_train[:N], separator=", ", max_line_width=90)
  print np.array2string(np.hstack((y_train_pred[:N], y_train_discr[:N])), separator=", ", max_line_width=90)
