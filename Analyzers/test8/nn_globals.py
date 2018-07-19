from nn_logging import getLogger

logger = getLogger()

# ______________________________________________________________________________
# Globals
nlayers = 12  # 5 (CSC) + 4 (RPC) + 3 (GEM)

nvariables = (nlayers * 6) + 8

discr_pt_cut = 14.

reg_pt_scale = 100.

discr_loss_weight = 1.

add_noise = True

infile_muon = '../test7/histos_tba.14.npz'

infile_pileup = '../test7/histos_tbd.14.npz'

# ______________________________________________________________________________
# Import all the libs
import os
import sys
from datetime import datetime
os.environ['KERAS_BACKEND'] = 'tensorflow'
OLD_STDOUT = sys.stdout

import numpy as np
np.random.seed(2023)
logger.info('Using numpy {0}'.format(np.__version__))

import keras
from keras import backend as K
from keras.models import Sequential, Model
from keras.layers import Dense, Activation, Dropout, Input, BatchNormalization
from keras import initializers, regularizers, optimizers, losses
#K.set_epsilon(1e-08)
logger.info('Using keras {0}'.format(keras.__version__))

import tensorflow as tf
logger.info('Using tensorflow {0}'.format(tf.__version__))

import sklearn
#from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
logger.info('Using sklearn {0}'.format(sklearn.__version__))

import matplotlib.pyplot as plt
from matplotlib import colors
#%matplotlib inline

