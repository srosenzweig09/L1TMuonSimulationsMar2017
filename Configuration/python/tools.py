import os

def loadFromFile(filename, fmt=''):
  """
  'filename' is the filename of the text file which contains the list of all the filenames.
  'fmt' is a formatting string to change the filenames. (default: '')
  """
  if not os.path.exists(filename):
    raise RuntimeError('Bad filename: %s' % filename)
  lines = tuple(open(filename))
  lines = [line.strip() for line in lines if not line.lstrip().startswith('#')]  # remove comment lines
  if fmt:  lines = [fmt % line for line in lines]
  return lines

