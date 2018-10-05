import numpy as np

class IncrementalStats(object):
  def __init__(self, n_features=1, dtype=np.float32):
    self._means = np.zeros(n_features, dtype=dtype)
    self._variances = np.zeros(n_features, dtype=dtype)
    self._n = 0
  
  def add(self, x):
    x = np.atleast_1d(x)
    assert(x.shape == self._means.shape)

    self._n += 1
    delta = (x - self._means)
    self._means += delta / float(self._n)
    self._variances += delta * (x - self._means)

  def mean(self):
    return self._means

  def variance(self, ddof=0):
    if self._n < 2:
      return np.nan
    else:
      return self._variances / float(self._n - ddof)
  
  def std(self, ddof=0):
    if self._n < 2:
      return np.nan
    else:
      return np.sqrt(self.variance(ddof=ddof))
