import numpy as np
import tensorflow as tf

from cnn_estimator import (superstrip_size, n_zones, rows_per_zone, n_rows, n_columns, n_channels, n_classes, dropout, learning_rate)

def parse_image_fn(pixels, channels):
  n = pixels.shape[0]
  assert(pixels.shape == (n,2))
  #assert(channels.shape == (n,n_channels))

  mask = (pixels[:,0] != -99)
  image_shape = (n_rows, n_columns, n_channels)
  image = np.zeros(image_shape, dtype=channels.dtype)
  #image[pixels[mask,0], pixels[mask,1]] = channels[mask]
  image[pixels[mask,0], pixels[mask,1]] = np.column_stack((channels[mask,0],1-channels[mask,0]))
  return image

def parse_label_fn(labels):
  assert(labels.shape == (3,))
  lb = labels[:1]
  lb = tf.keras.utils.to_categorical(lb[0], num_classes=n_classes)  # in Keras, use one-hot encoding
  return lb


# ______________________________________________________________________________
def create_model(params={}):
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

  if data_format == 'channels_first':
    input_shape = [n_channels, n_rows, n_columns]
  else:
    assert data_format == 'channels_last'
    input_shape = [n_rows, n_columns, n_channels]

  # ____________________________________________________________________________
  l = tf.keras.layers

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
  x = l.BatchNormalization(momentum=0.9,epsilon=1e-4)(x)
  x = l.MaxPooling2D(
          (1,5),
          strides=(1,3),
          padding='valid')(x)
  x = l.MaxPooling2D(
          (1,3),
          strides=(1,2),
          padding='valid')(x)
  x = l.Flatten()(x)
  x = l.Dense(40, use_bias=False)(x)
  x = l.BatchNormalization(momentum=0.9,epsilon=1e-4)(x)
  x = l.Activation('relu')(x)
  #x = l.Dropout(dropout)(x)
  x = l.Dense(12, use_bias=False)(x)
  x = l.BatchNormalization(momentum=0.9,epsilon=1e-4)(x)
  x = l.Activation('relu')(x)
  #x = l.Dropout(dropout)(x)
  x = l.Dense(n_classes, activation='softmax')(x)  # in Keras, get softmax outputs instead of logits

  outputs = x

  model = tf.keras.Model(inputs=inputs, outputs=outputs)

  model.compile(optimizer=tf.keras.optimizers.Adam(lr=learning_rate),
                loss='categorical_crossentropy',
                metrics=['accuracy'])

  model.summary()
  return model


# ______________________________________________________________________________
def define_reiam_flags():
  """Define flags which will be used for MNIST models."""
  from cnn_utils import define_reiam_base_flags, clear_flags, unparse_flags, set_defaults

  define_reiam_base_flags()
  set_defaults(data_dir='./reiam_data',
               model_dir='./reiam_model',
               batch_size=50,
               num_epochs=10)


# ______________________________________________________________________________
def run_reiam(flags_obj, data):
  """Run MNIST training and eval loop.
  Args:
    flags_obj: An object containing parsed flag values.
  """

  assert(tf.keras.backend.backend() == 'tensorflow')
  assert(tf.keras.backend.image_data_format() == 'channels_last')

  (images_px_train, images_px_test, images_ch_train, images_ch_test, labels_train, labels_test, parameters_train, parameters_test) = data

  images = map(lambda (px, ch): parse_image_fn(px, ch), zip(images_px_train, images_ch_train))
  labels = map(lambda lb: parse_label_fn(lb), labels_train)
  images = np.asarray(images)
  labels = np.asarray(labels)

  params={
      #'data_format': data_format,
      #'multi_gpu': multi_gpu,
      'n_rows': n_rows,
      'n_columns': n_columns,
      'n_channels': n_channels,
      'n_classes': n_classes,
      'dropout': dropout,
      'learning_rate': learning_rate,
  }

  model = create_model(params)

  # ____________________________________________________________________________
  import datetime
  import sys

  from nn_logging import getLogger
  logger = getLogger()

  from nn_training import TrainingLog

  start_time = datetime.datetime.now()
  logger.info('Begin training ...')

  # Redirect sys.stdout
  log = TrainingLog()
  sys.stdout = log

  history = model.fit(images, labels, epochs=flags_obj.num_epochs, batch_size=flags_obj.batch_size, validation_split=0.1, verbose=1)

  # Restore sys.stdout
  sys.stdout = log.stdout

  logger.info('Done training. Time elapsed: {0} sec'.format(str(datetime.datetime.now() - start_time)))

  from nn_models import save_my_model
  save_my_model(model, name='model_cnn')

  # ____________________________________________________________________________
  images = map(lambda (px, ch): parse_image_fn(px, ch), zip(images_px_test, images_ch_test))
  labels = map(lambda lb: parse_label_fn(lb), labels_test)
  images = np.asarray(images)
  labels = np.asarray(labels)

  loss_and_metrics = model.evaluate(images, labels, batch_size=flags_obj.batch_size, verbose=0)
  logger.info('loss and metrics: {0}'.format(loss_and_metrics))
  return model
