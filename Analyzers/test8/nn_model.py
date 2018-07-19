import os
os.environ['KERAS_BACKEND'] = 'tensorflow'

from keras import backend as K
import tensorflow as tf

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
def create_model(lr=0.00113):
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
def create_model_sequential(lr=0.001):
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
def save_model(model, name='model'):
  # Store model to file
  import h5py
  model.summary()
  model.save(name + '.h5')
  model.save_weights(name + '_weights.h5')

  # Store model to json
  import json
  with open(name + '.json', 'w') as outfile:
    outfile.write(model.to_json())
  return

