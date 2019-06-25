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
  infile_pileup = model_file.replace('model', 'histos_tbd').replace('.json', '.npz')

  # Load data
  import numpy as np
  np.random.seed(2023)
  with np.load(infile_pileup) as loaded:
    the_variables = loaded['variables']
    the_parameters = loaded['parameters']
    the_aux = loaded['aux']
  print('Loaded the variables with shape {0} and the parameters with shape {1}'.format(the_variables.shape, the_parameters.shape))
  print('Loaded the auxiliary PU info with shape {0}'.format(the_aux.shape))

  # Data preprocessing
  from nn_encode import nlayers, nvariables, nvariables_input, nparameters_input, Encoder
  reg_pt_scale = 100.
  reg_dxy_scale = 0.4

  encoder = Encoder(the_variables, the_parameters, reg_pt_scale=reg_pt_scale, reg_dxy_scale=reg_dxy_scale)
  x, y, aux = encoder.get_x(), encoder.get_y(), the_aux

  # Split dataset in training and testing
  from sklearn.model_selection import train_test_split
  test_job = 159

  split = aux[:,0].astype(np.int32) < test_job
  pu_x_train, pu_x_test, pu_y_train, pu_y_test, pu_aux_train, pu_aux_test = x[split], x[~split], y[split], y[~split], aux[split], aux[~split]
  print('Loaded # of training and testing events (PU): {0}'.format((pu_x_train.shape[0], pu_x_test.shape[0])))

  # Load model
  from nn_models import load_my_model, update_keras_custom_objects
  update_keras_custom_objects()
  loaded_model = load_my_model(name=model_file, weights_name=model_weights_file)
  print('Loaded model.')

  # Run
  print('Applying the model. This will take a while ...')
  pu_y_train_true = pu_y_train[:,np.newaxis].copy()
  pu_y_train_pred, pu_y_train_discr = loaded_model.predict(pu_x_train)
  #pu_y_train_true /= reg_pt_scale
  #pu_y_train_pred /= reg_pt_scale

  # Purge
  discr_pt_cut_low = 4.
  mask_train = ((pu_y_train_true == 0.) | (np.abs(1.0/pu_y_train_true) < discr_pt_cut_low/reg_pt_scale)) & (pu_y_train_discr > 0.7733)
  mask_train = mask_train[...,0]

  mask = np.zeros_like(y, dtype=np.bool)  # 'mask' is for the entire x, y, aux arrays
  mask[split] = mask_train
  assert((mask_train==1).sum() == (mask==1).sum())
  print('Removed {0} tracks from training.'.format((mask_train==1).sum()))

  outfile = infile_pileup.replace('histos_tbd', 'histos_tbd_purged')
  out_parameters = the_parameters[~mask]
  out_variables = the_variables[~mask]
  out_aux = the_aux[~mask]
  np.savez_compressed(outfile, parameters=out_parameters, variables=out_variables, aux=out_aux)
  print('Output: {0}'.format(outfile))
