# Based on official TensorFlow implementation of MNIST
# https://github.com/tensorflow/models/tree/master/official/mnist

import numpy as np
import tensorflow as tf

# ______________________________________________________________________________
# Globals

superstrip_size = 16

n_zones = 7

rows_per_zone = 11

n_rows = n_zones * rows_per_zone

n_columns = 5040 // superstrip_size

n_channels = 2

n_classes = 21

dropout = 0.2

learning_rate = 1e-3


# ______________________________________________________________________________
def parse_image_fn(pixels, channels):
  n = pixels.shape[0]
  assert(pixels.shape == (n,2))
  #assert(channels.shape == (n,n_channels))

  bad_pixel = -99
  mask = tf.not_equal(pixels[:,0], bad_pixel)

  indices = tf.boolean_mask(pixels, mask)
  updates = tf.boolean_mask(channels, mask)
  updates = tf.stack([updates[:,0], 1-updates[:,0]], axis=1)
  image_shape = (n_rows, n_columns, n_channels)
  image = tf.scatter_nd(indices, updates, image_shape)
  return image

def parse_label_fn(labels):
  assert(labels.shape == (3,))
  lb = labels[0]
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
  #x = l.BatchNormalization(momentum=0.9,epsilon=1e-4)(x)
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
  #x = l.BatchNormalization(momentum=0.9,epsilon=1e-4)(x)
  x = l.Activation('relu')(x)
  #x = l.Dropout(dropout)(x)
  x = l.Dense(12, use_bias=False)(x)
  #x = l.BatchNormalization(momentum=0.9,epsilon=1e-4)(x)
  x = l.Activation('relu')(x)
  #x = l.Dropout(dropout)(x)
  x = l.Dense(n_classes)(x)

  outputs = x

  model = tf.keras.Model(inputs=inputs, outputs=outputs)

  #model.compile(optimizer=tf.keras.optimizers.Adam(lr=learning_rate),
  #              loss='categorical_crossentropy',
  #              metrics=['accuracy'])

  #model.summary()
  return model


# ______________________________________________________________________________
def model_fn(features, labels, mode, params):
  """The model_fn argument for creating an Estimator."""

  # Set parameters
  learning_rate = params['learning_rate']

  model = create_model(params)
  image = features
  if isinstance(image, dict):
    image = features['image']

  if mode == tf.estimator.ModeKeys.PREDICT:
    try:
      logits = model(image, training=False)
    except TypeError:  # no keyword argument 'training' in tensorflow 1.5
      tf.keras.backend.set_learning_phase(False)
      logits = model(image)
    predictions = {
        'classes': tf.argmax(logits, axis=1),
        'probabilities': tf.nn.softmax(logits),
    }
    export_outputs = {
        'classify': tf.estimator.export.PredictOutput(predictions)
    }
    # For mode == ModeKeys.PREDICT: required fields are predictions.
    return tf.estimator.EstimatorSpec(
        mode=tf.estimator.ModeKeys.PREDICT,
        predictions=predictions,
        export_outputs=export_outputs)

  if mode == tf.estimator.ModeKeys.TRAIN:
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate)

    # If we are running multi-GPU, we need to wrap the optimizer.
    if params.get('multi_gpu'):
      optimizer = tf.contrib.estimator.TowerOptimizer(optimizer)

    try:
      logits = model(image, training=True)
    except TypeError:  # no keyword argument 'training' in tensorflow 1.5
      tf.keras.backend.set_learning_phase(True)
      logits = model(image)
    loss = tf.nn.sparse_softmax_cross_entropy_with_logits(labels=labels, logits=logits)
    loss = tf.reduce_mean(loss)
    accuracy = tf.metrics.accuracy(
        labels=labels, predictions=tf.argmax(logits, axis=1))
    accuracy_at_k = tf.metrics.mean(tf.to_float(tf.nn.in_top_k(
        targets=labels, predictions=logits, k=2)))

    # Name tensors to be logged with LoggingTensorHook.
    tf.identity(learning_rate, 'learning_rate')
    tf.identity(loss, 'cross_entropy')
    tf.identity(accuracy[1], name='train_accuracy')
    tf.identity(accuracy_at_k[1], name='train_accuracy_at_k')

    # Save accuracy scalar to Tensorboard output.
    tf.summary.scalar('train_accuracy', accuracy[1])
    tf.summary.scalar('train_accuracy_at_k', accuracy_at_k[1])

    # For mode == ModeKeys.TRAIN: required fields are loss and train_op
    return tf.estimator.EstimatorSpec(
        mode=tf.estimator.ModeKeys.TRAIN,
        loss=loss,
        train_op=optimizer.minimize(loss, global_step=tf.train.get_or_create_global_step()))

  if mode == tf.estimator.ModeKeys.EVAL:
    try:
      logits = model(image, training=False)
    except TypeError:  # no keyword argument 'training' in tensorflow 1.5
      tf.keras.backend.set_learning_phase(False)
      logits = model(image)
    loss = tf.nn.sparse_softmax_cross_entropy_with_logits(labels=labels, logits=logits)
    loss = tf.reduce_mean(loss)
    accuracy = tf.metrics.accuracy(
        labels=labels, predictions=tf.argmax(logits, axis=1))
    accuracy_at_k = tf.metrics.mean(tf.to_float(tf.nn.in_top_k(
        targets=labels, predictions=logits, k=2)))

    # Save accuracy scalar to Tensorboard output.
    tf.summary.scalar('eval_accuracy', accuracy[1])
    tf.summary.scalar('eval_accuracy_at_k', accuracy_at_k[1])

    # For mode == ModeKeys.EVAL: required field is loss.
    return tf.estimator.EstimatorSpec(
        mode=tf.estimator.ModeKeys.EVAL,
        loss=loss,
        eval_metric_ops={'accuracy': accuracy, 'accuracy_at_k': accuracy_at_k})


# ______________________________________________________________________________
# Patch the function _get_features_and_labels_from_input_fn()

# from tensorflow/python/training/basic_session_run_hooks.py
class FeedFnHook(tf.train.SessionRunHook):
  """Runs `feed_fn` and sets the `feed_dict` accordingly."""

  def __init__(self, feed_fn):
    """Initializes a `FeedFnHook`.
    Args:
      feed_fn: function that takes no arguments and returns `dict` of `Tensor`
        to feed.
    """
    self.feed_fn = feed_fn

  def before_run(self, run_context):  # pylint: disable=unused-argument
    return tf.train.SessionRunArgs(
        fetches=None, feed_dict=self.feed_fn())

# from tensorflow/python/estimator/estimator.py
class _DatasetInitializerHook(tf.train.SessionRunHook):
  """Creates a SessionRunHook that initializes the passed iterator."""

  def __init__(self, iterator, feed_fn):
    self._iterator = iterator
    self.feed_fn = feed_fn

  def begin(self):
    self._initializer = self._iterator.initializer

  def after_create_session(self, session, coord):
    del coord
    session.run(self._initializer, feed_dict=self.feed_fn())

# from tensorflow/python/estimator/estimator.py
def _get_features_and_labels_from_input_fn(self, input_fn, mode):
  """Extracts the `features` and labels from return values of `input_fn`."""
  result = self._call_input_fn(input_fn, mode)
  input_hooks = []
  if isinstance(result, tf.data.Dataset):
    iterator = result.make_initializable_iterator()
    #input_hooks.append(_DatasetInitializerHook(iterator))
    if mode == tf.estimator.ModeKeys.TRAIN:
      input_hooks.append(_DatasetInitializerHook(iterator, self.train_input_hook.feed_fn))
    else:  # mode == tf.estimator.ModeKeys.EVAL
      input_hooks.append(_DatasetInitializerHook(iterator, self.eval_input_hook.feed_fn))
    result = iterator.get_next()
  if isinstance(result, (list, tuple)):
    if len(result) != 2:
      raise ValueError(
          'input_fn should return (feautures, labels) as a len 2 tuple.')
    return result[0], result[1], input_hooks
  return result, None, input_hooks


# ______________________________________________________________________________
def define_reiam_flags():
  """Define flags which will be used for MNIST models."""
  from cnn_utils import define_reiam_base_flags, clear_flags, unparse_flags, set_defaults

  define_reiam_base_flags()
  set_defaults(batch_size=50,
               num_epochs=5,
               data_dir='./reiam_data',
               model_dir='./reiam_model',
               benchmark_logger_type='BenchmarkFileLogger',
               benchmark_log_dir='./reiam_model',
               hooks='LoggingTensorHook,ProfilerHook,ExamplesPerSecondHook,LoggingMetricHook')


# ______________________________________________________________________________
def save_model(reiam_classifier, model_name='model_cnn'):
  from tensorflow.python.training import saver as saver_lib
  checkpoint_path = saver_lib.latest_checkpoint(reiam_classifier._model_dir)
  if not checkpoint_path:
    raise ValueError('Could not find trained model in model_dir: {}.'.
                     format(reiam_classifier._model_dir))
  reader = tf.train.NewCheckpointReader(checkpoint_path)

  #keys = []
  #for key in reader.get_variable_to_shape_map():
  #  if key == 'global_step':
  #    continue
  #  keys.append(key)
  #print len(keys), keys

  keras_model = reiam_classifier._keras_model
  names = [weight.name for layer in keras_model.layers for weight in layer.weights]
  weights = keras_model.get_weights()
  assert(len(names) == len(weights))

  for i, (name, weight) in enumerate(zip(names, weights)):
    key = str(name).strip(':0')
    arr = reader.get_tensor(key)
    assert(weights[i].shape == arr.shape)
    weights[i] = arr
  keras_model.set_weights(weights)

  from nn_models import save_my_model
  save_my_model(keras_model, name=model_name)


# ______________________________________________________________________________
def run_reiam(flags_obj, data):
  """Run MNIST training and eval loop.
  Args:
    flags_obj: An object containing parsed flag values.
  """

  (images_px_train, images_px_test, images_ch_train, images_ch_test, labels_train, labels_test, parameters_train, parameters_test) = data

  # ____________________________________________________________________________
  class ModelHelpers(object):
    def __init__(self):
      from cnn_utils import past_stop_threshold, apply_clean
      self.past_stop_threshold = past_stop_threshold
      self.apply_clean = apply_clean

  model_helpers = ModelHelpers()

  class HooksHelper(object):
    def __init__(self):
      from cnn_utils import get_train_hooks
      self.get_train_hooks = get_train_hooks

  hooks_helper = HooksHelper()

  model_helpers.apply_clean(flags_obj)
  model_function = model_fn

  # Get number of GPUs as defined by the --num_gpus flags and the number of
  # GPUs available on the machine.
  num_gpus = flags_obj.num_gpus
  multi_gpu = num_gpus > 1

  if multi_gpu:
    # Validate that the batch size can be split into devices.
    distribution_utils.per_device_batch_size(flags_obj.batch_size, num_gpus)

    # There are two steps required if using multi-GPU: (1) wrap the model_fn,
    # and (2) wrap the optimizer. The first happens here, and (2) happens
    # in the model_fn itself when the optimizer is defined.
    model_function = tf.contrib.estimator.replicate_model_fn(
        model_fn, loss_reduction=tf.losses.Reduction.MEAN,
        devices=["/device:GPU:%d" % d for d in range(num_gpus)])

  data_format = flags_obj.data_format
  if data_format is None:
    data_format = ('channels_first'
                   if tf.test.is_built_with_cuda() else 'channels_last')

  params = {
      'data_format': data_format,
      'multi_gpu': multi_gpu,
      'n_rows': n_rows,
      'n_columns': n_columns,
      'n_channels': n_channels,
      'n_classes': n_classes,
      'dropout': dropout,
      'learning_rate': learning_rate,
  }

  # Create Estimator
  sess_config = tf.ConfigProto(
      log_device_placement=flags_obj.log_device_placement,
      intra_op_parallelism_threads=flags_obj.num_intra_threads,
      inter_op_parallelism_threads=flags_obj.num_inter_threads,
      allow_soft_placement=True)

  model_config = tf.estimator.RunConfig(
      model_dir=flags_obj.model_dir,
      tf_random_seed=flags_obj.tf_random_seed,
      save_summary_steps=100,
      save_checkpoints_secs=1200,
      session_config=sess_config)

  reiam_classifier = tf.estimator.Estimator(
      model_fn=model_function,
      model_dir=flags_obj.model_dir,
      config=model_config,
      params=params)

  reiam_classifier._keras_model = create_model(params)

  # Set up training and evaluation input functions.
  def get_train_input_fn_and_hook():
    feed_fn_hook = tf.train.FeedFnHook(feed_fn=None)

    def input_fn():
      with tf.name_scope('train_data'):
        get_shape_ph = lambda x: [None] + list(x[1:])
        images_px_ph         = tf.placeholder(images_px_train.dtype, get_shape_ph(images_px_train.shape))
        images_ch_ph         = tf.placeholder(images_ch_train.dtype, get_shape_ph(images_ch_train.shape))
        labels_ph            = tf.placeholder(labels_train.dtype, get_shape_ph(labels_train.shape))
        parameters_ph        = tf.placeholder(parameters_train.dtype, get_shape_ph(parameters_train.shape))
        feed_dict            = {images_px_ph: images_px_train, images_ch_ph: images_ch_train, labels_ph: labels_train}
        feed_fn_hook.feed_fn = lambda: feed_dict

        dataset1 = tf.data.Dataset.from_tensor_slices((images_px_ph, images_ch_ph))
        dataset1 = dataset1.map(map_func=parse_image_fn)
        dataset2 = tf.data.Dataset.from_tensor_slices((labels_ph))
        dataset2 = dataset2.map(map_func=parse_label_fn)
        dataset = tf.data.Dataset.zip((dataset1, dataset2))

        if flags_obj.shuffle:
          dataset = dataset.shuffle(buffer_size=flags_obj.shuffle_buffer_size)
        dataset = dataset.batch(batch_size=flags_obj.batch_size)
        if flags_obj.prefetch:
          dataset = dataset.prefetch(buffer_size=flags_obj.prefetch_buffer_size)
        dataset = dataset.repeat(flags_obj.epochs_between_evals)
        return dataset
    return input_fn, feed_fn_hook

  def get_eval_input_fn_and_hook():
    feed_fn_hook = tf.train.FeedFnHook(feed_fn=None)

    def input_fn():
      with tf.name_scope('eval_data'):
        get_shape_ph = lambda x: [None] + list(x[1:])
        images_px_ph         = tf.placeholder(images_px_test.dtype, get_shape_ph(images_px_test.shape))
        images_ch_ph         = tf.placeholder(images_ch_test.dtype, get_shape_ph(images_ch_test.shape))
        labels_ph            = tf.placeholder(labels_test.dtype, get_shape_ph(labels_test.shape))
        parameters_ph        = tf.placeholder(parameters_test.dtype, get_shape_ph(parameters_test.shape))
        feed_dict            = {images_px_ph: images_px_test, images_ch_ph: images_ch_test, labels_ph: labels_test}
        feed_fn_hook.feed_fn = lambda: feed_dict

        dataset1 = tf.data.Dataset.from_tensor_slices((images_px_ph, images_ch_ph))
        dataset1 = dataset1.map(map_func=parse_image_fn)
        dataset2 = tf.data.Dataset.from_tensor_slices((labels_ph))
        dataset2 = dataset2.map(map_func=parse_label_fn)
        dataset = tf.data.Dataset.zip((dataset1, dataset2))

        if flags_obj.shuffle:
          dataset = dataset.shuffle(buffer_size=flags_obj.shuffle_buffer_size)
        dataset = dataset.batch(batch_size=flags_obj.batch_size)
        if flags_obj.prefetch:
          dataset = dataset.prefetch(buffer_size=flags_obj.prefetch_buffer_size)
        dataset = dataset.repeat(flags_obj.epochs_between_evals)
        return dataset
    return input_fn, feed_fn_hook

  train_input_fn, train_input_hook = get_train_input_fn_and_hook()

  eval_input_fn, eval_input_hook = get_eval_input_fn_and_hook()

  # Set up hook that outputs training logs every 100 steps.
  train_hooks = hooks_helper.get_train_hooks(
      flags_obj.hooks, model_dir=flags_obj.model_dir,
      batch_size=flags_obj.batch_size)

  eval_hooks = []

  # Patch the function _get_features_and_labels_from_input_fn()
  import types
  reiam_classifier.train_input_fn = train_input_fn
  reiam_classifier.train_input_hook = train_input_hook
  reiam_classifier.eval_input_fn = eval_input_fn
  reiam_classifier.eval_input_hook = eval_input_hook
  reiam_classifier._get_features_and_labels_from_input_fn = types.MethodType(_get_features_and_labels_from_input_fn, reiam_classifier)

  # Train and evaluate model.
  for epoch in range(flags_obj.num_epochs // flags_obj.epochs_between_evals):
    reiam_classifier.train(input_fn=train_input_fn, hooks=train_hooks)
    eval_results = reiam_classifier.evaluate(input_fn=eval_input_fn, hooks=eval_hooks)
    #print('Epoch %i/%i evaluation results:\n\t%s\n' % (epoch+1, flags_obj.num_epochs, eval_results))
    if model_helpers.past_stop_threshold(flags_obj.stop_threshold, eval_results['accuracy']):
      break

  ## Train and evaluate model.
  #train_spec = tf.estimator.TrainSpec(input_fn=train_input_fn, hooks=train_hooks, max_steps=None)
  #eval_spec = tf.estimator.EvalSpec(input_fn=eval_input_fn, hooks=eval_hooks, steps=None)
  #for epoch in range(flags_obj.num_epochs // flags_obj.epochs_between_evals):
  #  tf.estimator.train_and_evaluate(reiam_classifier, train_spec, eval_spec)

  # ____________________________________________________________________________
  ## Export the model
  #if flags_obj.export_dir is not None:
  #  image = tf.placeholder(tf.float32, [None, n_rows, n_columns, n_channels])
  #  input_fn = tf.estimator.export.build_raw_serving_input_receiver_fn({
  #      'image': image,
  #  })
  #  reiam_classifier.export_savedmodel(flags_obj.export_dir, input_fn)

  # Export the model
  save_model(reiam_classifier)

  return reiam_classifier
