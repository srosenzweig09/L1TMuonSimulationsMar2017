import numpy as np
import tensorflow as tf

from cnn_estimator import (superstrip_size, n_zones, n_rows, n_columns, n_channels, n_classes, dropout, learning_rate)

def parse_image_fn(pixels, channels):
  n = pixels.shape[0]
  assert(pixels.shape == (n,2))
  assert(channels.shape == (n,n_channels))

  mask = (pixels[:,0] != -99)
  image_shape = (n_rows, n_columns, n_channels)
  image = np.zeros(image_shape, dtype=channels.dtype)
  image[pixels[mask,0], pixels[mask,1]] = channels[mask]
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

  l = tf.keras.layers

  # The model consists of a sequential chain of layers, so tf.keras.Sequential
  # (a subclass of tf.keras.Model) makes for a compact description.
  model = tf.keras.Sequential(
      [
          #l.Reshape(
          #    target_shape=input_shape,
          #    input_shape=(n_rows * n_columns,)),
          l.Conv2D(
              32,
              5,
              input_shape=input_shape,
              padding='same',
              data_format=data_format,
              activation='relu'),
          l.MaxPooling2D(
              (2, 2),
              (2, 2),
              padding='same',
              data_format=data_format),
          l.Conv2D(
              64,
              5,
              padding='same',
              data_format=data_format,
              activation='relu'),
          l.MaxPooling2D(
              (2, 2),
              (2, 2),
              padding='same',
              data_format=data_format),
          l.Flatten(),
          l.Dense(1024, activation='relu'),
          #l.Dropout(dropout),
          l.Dense(n_classes, activation='softmax'),  # in Keras, get softmax outputs instead of logits
      ])

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
               batch_size=50,  # use 32?
               num_epochs=1)


# ______________________________________________________________________________
def run_reiam(flags_obj, data):
  """Run MNIST training and eval loop.
  Args:
    flags_obj: An object containing parsed flag values.
  """

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

  from nn_models import save_my_model

  from nn_training import TrainingLog

  start_time = datetime.datetime.now()
  logger.info('Begin training ...')

  # Redirect sys.stdout
  log = TrainingLog()
  sys.stdout = log

  history = model.fit(images, labels, epochs=flags_obj.num_epochs, batch_size=flags_obj.batch_size, verbose=1)

  # Restore sys.stdout
  sys.stdout = log.stdout

  save_my_model(model, name='model_cnn')

  logger.info('Done training. Time elapsed: {0} sec'.format(str(datetime.datetime.now() - start_time)))

  # ____________________________________________________________________________
  images = map(lambda (px, ch): parse_image_fn(px, ch), zip(images_px_test, images_ch_test))
  labels = map(lambda lb: parse_label_fn(lb), labels_test)
  images = np.asarray(images)
  labels = np.asarray(labels)

  loss_and_metrics = model.evaluate(images, labels, batch_size=flags_obj.batch_size, verbose=0)
  logger.info('loss and metrics: {0}'.format(loss_and_metrics))
  return
