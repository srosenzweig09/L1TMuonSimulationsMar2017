import numpy as np

nlayers = 12  # 5 (CSC) + 4 (RPC) + 3 (GEM)

nvariables = (nlayers * 6) + 8

class Encoder(object):

  def __init__(self, x, y, adjust_scale=0):
    if x is not None and y is not None:
      assert(x.shape[1] == (nlayers * 7) + 3)
      assert(y.shape[1] == 3)
      assert(x.shape[0] == y.shape[0])

      self.nentries = x.shape[0]
      self.x_orig  = x
      self.y_orig  = y
      self.x_copy  = x.copy()
      self.y_copy  = y.copy()

      # Get views
      self.x_phi   = self.x_copy[:, nlayers*0:nlayers*1]
      self.x_theta = self.x_copy[:, nlayers*1:nlayers*2]
      self.x_bend  = self.x_copy[:, nlayers*2:nlayers*3]
      self.x_time  = self.x_copy[:, nlayers*3:nlayers*4]
      self.x_ring  = self.x_copy[:, nlayers*4:nlayers*5]
      self.x_fr    = self.x_copy[:, nlayers*5:nlayers*6]
      self.x_mask  = self.x_copy[:, nlayers*6:nlayers*7].astype(np.bool)  # this makes a copy
      self.x_road  = self.x_copy[:, nlayers*7:nlayers*8]  # ipt, ieta, iphi
      self.y_pt    = self.y_copy[:, 0]  # q/pT
      self.y_phi   = self.y_copy[:, 1]
      self.y_eta   = self.y_copy[:, 2]

      # Make event weight
      #self.w       = np.ones(self.y_pt.shape, dtype=np.float32)
      self.w       = np.abs(self.y_pt)/0.2 + 1.0

      # Straightness & zone
      self.x_straightness = self.x_road[:, 0][:, np.newaxis]
      self.x_zone         = self.x_road[:, 1][:, np.newaxis]

      # Subtract median phi from hit phis
      #self.x_phi_median    = self.x_road[:, 2] * 32 - 16  # multiply by 'quadstrip' unit (4 * 8)
      self.x_phi_median    = self.x_road[:, 2] * 16 - 8  # multiply by 'doublestrip' unit (2 * 8)
      self.x_phi_median    = self.x_phi_median[:, np.newaxis]
      self.x_phi          -= self.x_phi_median

      # Subtract median theta from hit thetas
      self.x_theta_median  = np.nanmedian(self.x_theta[:,:5], axis=1)  # CSC only
      self.x_theta_median[np.isnan(self.x_theta_median)] = np.nanmedian(self.x_theta[np.isnan(self.x_theta_median)], axis=1)  # use all
      self.x_theta_median  = self.x_theta_median[:, np.newaxis]
      self.x_theta        -= self.x_theta_median

      # Standard scales
      # + Remove outlier hits by checking hit thetas
      if adjust_scale == 0:  # do not adjust
        x_theta_tmp = np.abs(self.x_theta) > 10000.0
      elif adjust_scale == 1:  # use mean and std
        self.x_mean  = np.nanmean(self.x_copy, axis=0)
        self.x_std   = np.nanstd(self.x_copy, axis=0)
        self.x_std   = self._handle_zero_in_scale(self.x_std)
        self.x_copy -= self.x_mean
        self.x_copy /= self.x_std
        x_theta_tmp = np.abs(self.x_theta) > 1.0
      elif adjust_scale == 2:  # adjust by hand
        theta_cuts    = np.array((6., 6., 6., 6., 6., 12., 12., 12., 12., 9., 9., 9.), dtype=np.float32)
        x_theta_tmp   = np.abs(self.x_theta) > theta_cuts
        self.x_phi   *= 0.000991  # GE1/1 dphi linear correlation with q/pT
        self.x_theta *= (1/12.)   # 12 integer theta units
        self.x_bend  *= 0.188082  # ME1/2 bend linear correlation with q/pT
        x_ring_tmp    = self.x_ring.astype(np.int32)
        x_ring_tmp    = (x_ring_tmp == 1) | (x_ring_tmp == 4)
        self.x_ring[x_ring_tmp] = 0  # ring 1,4 -> 0
        self.x_ring[~x_ring_tmp] = 1 # ring 2,3 -> 1
        x_fr_tmp      = self.x_fr.astype(np.int32)
        x_fr_tmp      = (x_fr_tmp == 0)
        self.x_fr[x_fr_tmp] = 0
        self.x_fr[~x_fr_tmp] = 1
      elif adjust_scale == 3:  # adjust by hand #2
        #theta_cuts    = np.array((6., 6., 6., 6., 6., 12., 12., 12., 12., 9., 9., 9.), dtype=np.float32)
        theta_cuts    = np.array((6., 6., 6., 6., 6., 10., 10., 10., 10., 8., 8., 8.), dtype=np.float32)
        x_theta_tmp   = np.abs(self.x_theta) > theta_cuts
        self.x_bend[:, 5:9] = 0  # do not use RPC bend
        x_ring_tmp    = self.x_ring.astype(np.int32)
        x_ring_tmp    = (x_ring_tmp == 1) | (x_ring_tmp == 4)
        self.x_ring[x_ring_tmp] = 0  # ring 1,4 -> 0
        self.x_ring[~x_ring_tmp] = 1 # ring 2,3 -> 1
        x_fr_tmp      = self.x_fr.astype(np.int32)
        x_fr_tmp      = (x_fr_tmp == 0)
        self.x_fr[x_fr_tmp] = 0
        self.x_fr[~x_fr_tmp] = 1
        s = [ 0.005907,  0.012209, -0.015324, -0.011308, -0.008402,  0.013371,
             -0.0267  , -0.00999 , -0.007698,  0.004138, -0.018402,  0.005181,
              0.603693,  0.580338,  1.456922,  1.50118 ,  1.06335 ,  0.21944 ,
              0.291049,  0.346278,  0.386362,  0.510937,  0.600389,  0.719824,
              0.826546,  0.498369,  1.65504 , -1.189355, -1.215245,  1.278912,
              1.276016, -1.002222, -1.025343,  0.720984,  1.459154,  1.      ,
              1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  0.590936,
              0.605932,  0.722654,  0.71606 ,  1.      ,  1.      ,  0.881102,
              1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,
              1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,
              1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,
              1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,
              1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,
              1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,
              1.      ,  1.      ,  1.      ]
        self.x_copy *= s

      # Remove outlier hits by checking hit thetas
      self.x_phi  [x_theta_tmp] = np.nan
      self.x_theta[x_theta_tmp] = np.nan
      self.x_bend [x_theta_tmp] = np.nan
      self.x_time [x_theta_tmp] = np.nan
      self.x_ring [x_theta_tmp] = np.nan
      self.x_fr   [x_theta_tmp] = np.nan
      self.x_mask [x_theta_tmp] = 1.0

      # Add variables: straightness, zone, theta_median and mode variables
      self.x_straightness = (self.x_straightness - 6.) / 6.  # scaled to [-1,1]
      self.x_zone         = (self.x_zone - 0.) / 5.  # scaled to [0,1]
      self.x_theta_median = (self.x_theta_median - 3.) / 83.  # scaled to [0,1]
      hits_to_station = np.array((5,1,2,3,4,1,2,3,4,5,2,5), dtype=np.int32)  # '5' denotes ME1/1
      assert(len(hits_to_station) == nlayers)
      self.x_mode_vars = np.zeros((self.nentries, 5), dtype=np.bool)
      self.x_mode_vars[:,0] = np.any(self.x_mask[:,hits_to_station == 5] == 0, axis=1)
      self.x_mode_vars[:,1] = np.any(self.x_mask[:,hits_to_station == 1] == 0, axis=1)
      self.x_mode_vars[:,2] = np.any(self.x_mask[:,hits_to_station == 2] == 0, axis=1)
      self.x_mode_vars[:,3] = np.any(self.x_mask[:,hits_to_station == 3] == 0, axis=1)
      self.x_mode_vars[:,4] = np.any(self.x_mask[:,hits_to_station == 4] == 0, axis=1)

      # Remove NaN
      #np.nan_to_num(self.x_copy, copy=False)
      self.x_copy[np.isnan(self.x_copy)] = 0.0

  # Copied from scikit-learn
  def _handle_zero_in_scale(self, scale):
    scale[scale == 0.0] = 1.0
    return scale

  def get_x(self):
    #x_new = self.x_phi
    x_new = np.hstack((self.x_phi, self.x_theta, self.x_bend, self.x_time, self.x_ring, self.x_fr, self.x_straightness, self.x_zone, self.x_theta_median, self.x_mode_vars))
    return x_new

  def get_x_mask(self):
    x_mask = self.x_mask.copy()
    return x_mask

  def get_y(self):
    y_new = self.y_pt.copy()
    return y_new

  def get_w(self):
    w_new = self.w.copy()
    return w_new

  def save_encoder(self, filepath):
    np.savez_compressed(filepath, x_mean=self.x_mean, x_std=self.x_std)

  def load_endcoder(self, filepath):
    loaded = np.load(filepath)
    self.x_mean = loaded['x_mean']
    self.x_std = loaded['x_std']



import tensorflow as tf
from keras import backend as K

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

def masked_huber_loss(y_true, y_pred, delta=1.345):
  x = K.abs(y_true - y_pred)
  squared_loss = 0.5*K.square(x)
  absolute_loss = delta * (x - 0.5*delta)
  #xx = K.switch(x < delta, squared_loss, absolute_loss)
  xx = tf.where(x < delta, squared_loss, absolute_loss)  # needed for tensorflow

  mask_value = 100.
  mask = K.not_equal(y_true, mask_value)
  mask = K.cast(mask, K.floatx())
  xx *= mask
  xx /= K.mean(mask)
  return K.mean(xx, axis=-1)

#def masked_huber_loss(y_true, y_pred, delta=1.345):
#  mask_value = 100.
#  mask_alpha = 0.02
#  mask_target = 0.5 * reg_pt_scale
#  mask = K.equal(y_true, mask_value)
#
#  #x = K.abs(y_true - y_pred)
#  x = tf.where(mask, mask_alpha * K.abs(mask_target - K.abs(y_pred)), K.abs(y_true - y_pred))
#  squared_loss = 0.5*K.square(x)
#  absolute_loss = delta * (x - 0.5*delta)
#  #xx = K.switch(x < delta, squared_loss, absolute_loss)
#  xx = tf.where(x < delta, squared_loss, absolute_loss)  # needed for tensorflow
#  return K.mean(xx, axis=-1)


# ______________________________________________________________________________
# Binary crossentropy
def masked_binary_crossentropy(y_true, y_pred, from_logits=False):
  target, output = y_true, y_pred

  # transform back to logits
  if not from_logits:
    output = K.clip(output, K.epsilon(), 1 - K.epsilon())
    output = K.log(output / (1 - output))

  xx =  tf.nn.sigmoid_cross_entropy_with_logits(labels=target, logits=output)
  #xx =  tf.nn.weighted_cross_entropy_with_logits(targets=target, logits=output, pos_weight=0.5)  # pos_weight < 1 decreases the false positive count

  mask_value = 100.
  mask = K.not_equal(y_true, mask_value)
  mask = K.cast(mask, K.floatx())
  xx *= mask
  xx /= K.mean(mask)
  return K.mean(xx, axis=-1)


# ______________________________________________________________________________
# Learning rate decay by epoch number
from keras.callbacks import LearningRateScheduler

def lr_schedule(epoch):
  if (epoch % 10) == 0:
    lr = K.get_value(model.optimizer.lr)
    K.set_value(model.optimizer.lr, lr*0.95)
    print("lr changed to {}".format(lr*0.95))
  return K.get_value(model.optimizer.lr)

lr_decay = LearningRateScheduler(lr_schedule)


# ______________________________________________________________________________
# Custom objects
from keras.utils.generic_utils import get_custom_objects

get_custom_objects().update({'masked_huber_loss': masked_huber_loss, 'masked_binary_crossentropy': masked_binary_crossentropy, 'NewLeakyReLU': NewLeakyReLU, 'NewTanh': NewTanh})
