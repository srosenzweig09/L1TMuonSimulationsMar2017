def muon_data():
  try:
    print('[INFO] Loading muon data ...')
    loaded = np.load(infile_muon)
    the_variables = loaded['variables']
    the_parameters = loaded['parameters']
    print('[INFO] Loaded the variables with shape {0}'.format(the_variables.shape))
    print('[INFO] Loaded the parameters with shape {0}'.format(the_parameters.shape))
  except:
    print('[ERROR] Failed to load data from file: {0}'.format(infile_muon))

  encoder = Encoder(the_variables, the_parameters, nlayers=nlayers, adjust_scale=3)
  x, y, w, x_mask = encoder.get_x(), encoder.get_y(), encoder.get_w(), encoder.get_x_mask()
  assert np.isfinite(x).all()

  # Split dataset in training and testing
  x_train, x_test, y_train, y_test, w_train, w_test, x_mask_train, x_mask_test = train_test_split(x, y, w, x_mask, test_size=0.3)
  return x_train, x_test, y_train, y_test, w_train, w_test, x_mask_train, x_mask_test


def pileup_data():
  try:
    print('[INFO] Loading pileup data ...')
    loaded = np.load(infile_pileup)
    the_variables = loaded['variables']
    the_parameters = np.zeros((the_variables.shape[0], 3), dtype=np.float32)
    the_auxiliaries = loaded['aux']
    print('[INFO] Loaded the variables with shape {0}'.format(the_variables.shape))
    print('[INFO] Loaded the auxiliary info with shape {0}'.format(the_auxiliaries.shape))
  except:
    print('[ERROR] Failed to load data from file: {0}'.format(infile_pileup))

  sel = the_auxiliaries[:,2] > discr_pt_cut
  the_variables = the_variables[~sel]
  the_parameters = the_parameters[~sel]
  the_auxiliaries = the_auxiliaries[~sel]

  encoder = Encoder(the_variables, the_parameters, nlayers=nlayers, adjust_scale=3)
  x, y, w, x_mask = encoder.get_x(), encoder.get_y(), encoder.get_w(), encoder.get_x_mask()
  aux = the_auxiliaries  # jobid, ievt, highest_part_pt, highest_track_pt
  assert np.isfinite(x).all()

  # Split dataset in training and testing
  split = the_auxiliaries[:,0] < 50.
  x_train, x_test, aux_train, aux_test = x[~split], x[split], aux[~split], aux[split]
  return x_train, x_test, aux_train, aux_test
