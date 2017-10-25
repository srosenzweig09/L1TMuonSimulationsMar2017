#!/bin/bash

#etas=( 1.2 1.3 )
etas=( 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 )

maxevents=5000

# Generate .py
for eta in "${etas[@]}"
do
  mkdir -p "$eta"
  sed "s@__MAXEVENTS__@$maxevents@g" pset_SingleMuon_Toy_0M.py.template | sed "s@__ETA__@$eta@g" > $eta/pset_SingleMuon_Toy_0M.py
done

# Run .py
for eta in "${etas[@]}"
do
  echo "cd $eta/; cmsRun pset_SingleMuon_Toy_0M.py >& pset_SingleMuon_Toy_0M.log & cd -"
done

