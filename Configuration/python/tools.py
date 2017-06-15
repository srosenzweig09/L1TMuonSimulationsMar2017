import os

def loadFromFile(filename):
  if not os.path.exists(filename):
    raise RuntimeError('Bad filename: %s' % filename)
  lines = tuple(open(filename))
  lines = [line.strip() for line in lines if not line.lstrip().startswith('#')]  # remove comment lines
  return lines

