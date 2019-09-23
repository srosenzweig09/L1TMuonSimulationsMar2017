def parse_args():
  import sys
  rargs = sys.argv[1:]
  largs = []
  line = []

  # Get the arguments
  while rargs:
    arg = rargs.pop(0)
    if '\"' in arg:
      arg = '\'' + arg + '\''
    if arg == '--':
      break
    elif arg[0:2] == '--':
      if line:
        largs.append(' '.join(line))
        del line[:]
      line.append(arg)
    elif arg[:1] == '-' and len(arg) > 1:
      if line:
        largs.append(' '.join(line))
        del line[:]
      line.append(arg)
    else:  # assume it is the value
      line.append(arg)

  if line:
    largs.append(' '.join(line))
    del line[:]

  # Sort the optional args
  # (refer to https://github.com/cms-sw/cmssw/blob/master/Configuration/Applications/python/Options.py)
  priorities = [
    '-s',
    '--step',
    '--conditions',
    '--eventcontent',
    '--datatier',
    '--era',
    '--beamspot',
    '--geometry',
    '--magfield',
    '--pileup_input',
    '--pileup',
    '--customise_commands',
    '--customise',
  ]

  rargs = largs[1:]
  largs = largs[:1]  # keep the positional argument

  for p in priorities:
    i = 0
    while i < len(rargs):
      arg = rargs[i]
      if arg.startswith(p):
        largs.append(arg)
        del rargs[i]
      else:
        i += 1

  # Add the remaining args
  largs += rargs

  # Print
  printme = ' \\\n    '.join(largs)
  printme = 'cmsDriver.py ' + printme
  print('')
  print(printme)
  return


# ______________________________________________________________________________
if __name__ == '__main__':
  parse_args()
