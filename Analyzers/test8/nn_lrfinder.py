"""
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_4_0
cd CMSSW_10_4_0/src
cmsenv
virtualenv venv
source venv/bin/activate
pip install -U pip
pip install scikit-optimize
cd ../..

# More info about sklearn/skopt
#   http://scikit-learn.org/stable/modules/cross_validation.html
#   http://scikit-learn.org/stable/modules/grid_search.html
#   https://github.com/scikit-optimize/scikit-optimize/blob/master/examples/sklearn-gridsearchcv-replacement.ipynb
"""


# ______________________________________________________________________________
from nn_globals import *

from nn_encode import nlayers, nvariables

from nn_data import (muon_data_split, pileup_data_split, mix_training_inputs)

from nn_models import (create_model_sequential_bn2, lr_decay, modelbestcheck, modelbestcheck_weights)

from nn_training import train_model

from nn_pruning import prune_model


# ______________________________________________________________________________
import skopt
logger.info('Using skopt {0}'.format(skopt.__version__))

infile_muon = 'histos_tba.23.npz'
infile_pileup = 'histos_tbd.23.npz'


# ______________________________________________________________________________
# Import muon data
# 'x' is the array of input variables, 'y' is the q/pT
x_train, x_test, y_train, y_test, w_train, w_test, x_mask_train, x_mask_test, x_road_train, x_road_test = \
      muon_data_split(infile_muon, reg_pt_scale=reg_pt_scale, test_size=0.31)

## Use ShuffleSplit as the CV iterator
#from nn_data import muon_data
#from sklearn.model_selection import ShuffleSplit
#x, y, w, x_mask, x_road = muon_data(infile_muon, reg_pt_scale=reg_pt_scale, correct_for_eta=False)
#cv = ShuffleSplit(n_splits=1, test_size=0.31)

# ______________________________________________________________________________
# Create LRFinder

from keras_lr_finder import LRFinder
epochs = 20
batch_size = 4096

model = create_model_sequential_bn2(nvariables=nvariables, lr=learning_rate, clipnorm=gradient_clip_norm,
                                    l1_reg=l1_reg, l2_reg=l2_reg,
                                    nodes1=30, nodes2=25, nodes3=20)

lr_finder = LRFinder(model)
lr_finder.find(x_train, y_train, 0.0001, 1, batch_size=batch_size, epochs=epochs)
print len(lr_finder.lrs)

#lr_finder.plot_loss(n_skip_beginning=100, n_skip_end=50)

#lr_finder.plot_loss_change(sma=100, n_skip_beginning=100, n_skip_end=50, y_lim=(-0.1, 0.1))

logger.info('DONE')
