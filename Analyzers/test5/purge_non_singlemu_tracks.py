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
  infile_displ = model_file.replace('model', 'histos_tba_displ1').replace('.json', '.npz')
  #infile_displ = model_file.replace('model', 'histos_tba_displ2').replace('.json', '.npz')

  # Load data
  import numpy as np
  np.random.seed(2023)
  with np.load(infile_displ) as loaded:
    the_variables = loaded['variables']
    the_parameters = loaded['parameters']
  print('Loaded the variables with shape {0} and the parameters with shape {1}'.format(the_variables.shape, the_parameters.shape))

  # Data preprocessing
  from nn_encode import nlayers, nvariables, nvariables_input, nparameters_input, Encoder
  reg_pt_scale = 100.
  reg_dxy_scale = 0.4

  encoder = Encoder(the_variables, the_parameters, reg_pt_scale=reg_pt_scale, reg_dxy_scale=reg_dxy_scale)
  x, y, x_mask = encoder.get_x(), encoder.get_y(), encoder.get_x_mask()

  from numba import njit

  @njit
  def pass_singlemu(x_mask):
    valid = ~x_mask
    b1 = valid[np.array([0,1,5,9,11])].any()  # ME1/1, ME1/2, RE1/2, GE1/1, ME0
    b2 = valid[np.array([2,6,10])].any()      # ME2, RE2, GE2/1
    b3 = valid[np.array([3,7])].any()         # ME3, RE3
    b4 = valid[np.array([4,8])].any()         # ME4, RE4
    mode = (b1 << 3) | (b2 << 2) | (b3 << 1) | (b4 << 0)
    passed = (mode == 11) | (mode == 13) | (mode == 14) | (mode == 15)
    return passed

  @njit
  def pass_singlemu_apply(x_mask, a):
    for i in range(x_mask.shape[0]):
      a[i] = pass_singlemu(x_mask[i])
    return

  #tmp = np.apply_along_axis(pass_singlemu, 1, x_mask)
  tmp = np.zeros((x_mask.shape[0],), dtype=np.bool)
  pass_singlemu_apply(x_mask, tmp)

  outfile = infile_displ.replace('histos_tba_displ1', 'histos_tba_displ1_purged')
  #outfile = infile_displ.replace('histos_tba_displ2', 'histos_tba_displ2_purged')
  out_parameters = the_parameters[tmp]
  out_variables = the_variables[tmp]
  print('Output the variables with shape {0} and the parameters with shape {1}'.format(out_variables.shape, out_parameters.shape))
  np.savez_compressed(outfile, parameters=out_parameters, variables=out_variables)
  print('Output: {0}'.format(outfile))
