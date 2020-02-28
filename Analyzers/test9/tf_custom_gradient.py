# Copyright 2017 The TensorFlow Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================
"""Decorator to overrides the gradient for a function.

   See https://www.tensorflow.org/api_docs/python/tf/custom_gradient
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import tensorflow as tf

from tensorflow.python.framework import ops
#from tensorflow.python.ops import array_ops
#from tensorflow.python.ops import gen_array_ops
from tensorflow.python.util import nest
from tensorflow.python.util import tf_decorator


# From https://github.com/tensorflow/tensorflow/blob/v1.7.0/tensorflow/python/ops/custom_gradient.py
def tf_custom_gradient(f):
  """Decorator to define a function with a custom gradient.
  """

  def decorated(*args, **kwargs):
    """Decorated function with custom gradient."""

    if kwargs:
      raise ValueError(
          "The custom_gradient decorator currently suports keywords "
          "arguments only when eager execution is enabled.")
    name = "CustomGradient-%s" % ops.uid()
    args = [ops.convert_to_tensor(x) for x in args]
    result, grad_fn = f(*args)
    flat_result = nest.flatten(result)
    all_tensors = flat_result + args

    @ops.RegisterGradient(name)
    def internal_grad_fn(unused_op, *result_grads):
      gradients = nest.flatten(grad_fn(*result_grads[:len(flat_result)]))
      # Need to return one value per input to the IdentityN, so pad the
      # gradients of the inputs of the custom_gradient function with the
      # gradients of the outputs as well.
      return ([None] * len(flat_result)) + gradients

    with ops.get_default_graph().gradient_override_map({"IdentityN": name}):
      all_tensors = tf.identity_n(all_tensors)
    return nest.pack_sequence_as(
        structure=result, flat_sequence=all_tensors[:len(flat_result)])

  return tf_decorator.make_decorator(f, decorated)

# From https://github.com/tensorflow/tensorflow/blob/v1.7.0/tensorflow/python/ops/gradients_test.py
def testCustomGradientTrivial():

  @tf_custom_gradient
  def MyIdentity(x):

    def Grad(dy):
      return [3 * dy]

    return x, Grad

  with ops.Graph().as_default():
    x = tf.constant(3.)
    y = MyIdentity(MyIdentity(x))
    dy = tf.gradients(y, x)[0]
    with tf.Session():
      assert(9. == dy.eval())

def testCustomGradient():

  @tf_custom_gradient
  def MyMultiply(x1, x2):
    result = x1 * x2

    def Grad(dy):
      # Switched the ordering here.
      return [dy * x1, dy * x2]

    return result, Grad

  with ops.Graph().as_default():
    x1 = tf.constant(3.)
    x2 = tf.constant(5.)
    y = MyMultiply(x1, x2)
    dy = tf.gradients(y, [x1, x2])
    with tf.Session() as sess:
      assert([3., 5.] == sess.run(dy))

def testCustomGradientErrors():

  @tf_custom_gradient
  def F(x):

    def Grad(_):
      raise RuntimeError("x")

    return x, Grad

  with ops.Graph().as_default():
    x = tf.constant(1.0)
    y = F(x)
    tf.gradients(y, x)


# ______________________________________________________________________________
if __name__ == "__main__":

  testCustomGradientTrivial()
  testCustomGradient()
  try:
    testCustomGradientErrors()
  except Exception as e:
    assert(issubclass(type(e), RuntimeError))
