import tensorflow as tf

# ______________________________________________________________________________
# from tensorflow/models: official/mnist/mnist.py

from cnn_globals import rows_per_zone

def create_model(params={}, use_cnn_keras=False):
  """Model to recognize digits in the MNIST dataset.

  Network structure is equivalent to:
  https://github.com/tensorflow/tensorflow/blob/r1.5/tensorflow/examples/tutorials/mnist/mnist_deep.py
  and
  https://github.com/tensorflow/models/blob/master/tutorials/image/mnist/convolutional.py

  But uses the tf.keras API.

  Args:
    data_format: Either 'channels_first' or 'channels_last'. 'channels_first' is
      typically faster on GPUs while 'channels_last' is typically faster on
      CPUs. See
      https://www.tensorflow.org/performance/performance_guide#data_formats

  Returns:
    A tf.keras.Model.
  """

  # Set parameters
  data_format = params.get('data_format', 'channels_last')
  n_rows = params.get('n_rows', 28)
  n_columns = params.get('n_columns', 28)
  n_channels = params.get('n_channels', 1)
  n_classes = params.get('n_classes', 10)
  dropout = params.get('dropout', 0.4)
  learning_rate = params.get('learning_rate', 0.01)

  if data_format == 'channels_first':
    input_shape = [n_channels, n_rows, n_columns]
  else:
    assert data_format == 'channels_last'
    input_shape = [n_rows, n_columns, n_channels]

  if use_cnn_keras:
    in_training_mode = None  # in Keras, it's figured out correctly.
  else:
    in_training_mode = True

  # ____________________________________________________________________________
  l = tf.keras.layers

  xavier_init = tf.keras.initializers.glorot_uniform()
  he_init = tf.keras.initializers.he_normal()

  inputs = l.Input(shape=input_shape, dtype='float32')

  inputs_zone0 = l.Lambda(lambda x: x[:, 0*rows_per_zone:(0+1)*rows_per_zone])(inputs)
  inputs_zone1 = l.Lambda(lambda x: x[:, 1*rows_per_zone:(1+1)*rows_per_zone])(inputs)
  inputs_zone2 = l.Lambda(lambda x: x[:, 2*rows_per_zone:(2+1)*rows_per_zone])(inputs)
  inputs_zone3 = l.Lambda(lambda x: x[:, 3*rows_per_zone:(3+1)*rows_per_zone])(inputs)
  inputs_zone4 = l.Lambda(lambda x: x[:, 4*rows_per_zone:(4+1)*rows_per_zone])(inputs)
  inputs_zone5 = l.Lambda(lambda x: x[:, 5*rows_per_zone:(5+1)*rows_per_zone])(inputs)
  inputs_zone6 = l.Lambda(lambda x: x[:, 6*rows_per_zone:(6+1)*rows_per_zone])(inputs)

  def first_conv2d_layer_fn(inputs):
    x = l.Conv2D(
            12,
            (11,63),
            strides=(11,4),
            padding='same',
            activation='relu')(inputs)
    return x

  x_zone0 = first_conv2d_layer_fn(inputs_zone0)
  x_zone1 = first_conv2d_layer_fn(inputs_zone1)
  x_zone2 = first_conv2d_layer_fn(inputs_zone2)
  x_zone3 = first_conv2d_layer_fn(inputs_zone3)
  x_zone4 = first_conv2d_layer_fn(inputs_zone4)
  x_zone5 = first_conv2d_layer_fn(inputs_zone5)
  x_zone6 = first_conv2d_layer_fn(inputs_zone6)

  x_all_zones = [x_zone0, x_zone1, x_zone2, x_zone3, x_zone4, x_zone5, x_zone6]
  x = l.Concatenate()(x_all_zones)

  x = l.Conv2D(
          48,
          (5,5),
          strides=(1,1),
          padding='same',
          activation='relu')(x)
  #x = l.BatchNormalization(momentum=0.99, epsilon=1e-3)(x, training=in_training_mode)
  x = l.MaxPooling2D(
          (1,5),
          strides=(1,3),
          padding='valid')(x)
  x = l.MaxPooling2D(
          (1,3),
          strides=(1,2),
          padding='valid')(x)
  x = l.Flatten()(x)
  x = l.Dense(40, use_bias=False, activation=None, kernel_initializer=he_init)(x)
  x = l.BatchNormalization(momentum=0.99, epsilon=1e-3)(x, training=in_training_mode)
  x = l.Activation('relu')(x)
  #x = l.Dropout(dropout)(x)
  x = l.Dense(12, use_bias=False, activation=None, kernel_initializer=he_init)(x)
  x = l.BatchNormalization(momentum=0.99, epsilon=1e-3)(x, training=in_training_mode)
  x = l.Activation('relu')(x)
  #x = l.Dropout(dropout)(x)
  if use_cnn_keras:
    x = l.Dense(n_classes, activation='softmax')(x)  # in Keras, get softmax outputs instead of logits
  else:
    x = l.Dense(n_classes)(x)

  outputs = x

  model = tf.keras.Model(inputs=inputs, outputs=outputs)

  if use_cnn_keras:
    # in Keras, compile with Keras optimizer and loss
    model.compile(optimizer=tf.keras.optimizers.Adam(lr=learning_rate),
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])
    model.summary()
  else:
    pass
  return model


# ______________________________________________________________________________
# Save/Load models

import h5py
import json

def save_my_model(model, name='model'):
  # Store model to file
  #model.summary()
  model.save(name + '.h5')
  model.save_weights(name + '_weights.h5')
  # Store model to json
  with open(name + '.json', 'w') as outfile:
    outfile.write(model.to_json())
  return

def load_my_model(name='model', weights_name='model_weights'):
  with open(name + '.json', 'r') as f:
    json_string = json.dumps(json.load(f))
    model = tf.keras.models.model_from_json(json_string)
  #model = tf.keras.models.load_model(name + '.h5')
  model.load_weights(weights_name + '.h5')
  return model
