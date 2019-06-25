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
  constant_graph_file = model_file.replace('model', 'model_graph').replace('.json', '.pb')

  # Load model
  from nn_models import load_my_model, update_keras_custom_objects
  update_keras_custom_objects()
  loaded_model = load_my_model(name=model_file, weights_name=model_weights_file)

  # Save as a constant graph
  import tensorflow as tf
  from keras import backend as K
  sess = K.get_session()

  _input = loaded_model.input.op.name
  outputs = [node.op.name for node in loaded_model.outputs]
  print('input: {}'.format(_input))
  print('outputs: {}'.format(outputs))  # names of output operations you want to use later
  constant_graph = tf.graph_util.convert_variables_to_constants(
      sess, sess.graph.as_graph_def(), outputs)
  tf.train.write_graph(constant_graph, './', constant_graph_file, as_text=False)
  print('constant graph: {}'.format(constant_graph_file))
