#scp histos_tba.14.npz histos_tbd.14.npz jlow@15.226.54.11:/scratch/CMS/L1MuonTrigger/P2_10_1_5/SingleMuon_Toy_2GeV/
#scp histos_tba.15.npz histos_tbd.15.npz jlow@15.226.54.11:/scratch/CMS/L1MuonTrigger/P2_10_1_5/SingleMuon_Toy_2GeV/
#scp histos_tba.16.npz histos_tbd.16.npz jlow@15.226.54.11:/scratch/CMS/L1MuonTrigger/P2_10_1_5/SingleMuon_Toy_2GeV/

tar czf hpe_default.tgz nn_*.py mykeras5.ipynb tdrstyle.mplstyle
scp hpe_default.tgz jlow@15.226.54.11:/home/uf/jlow/jftest2/slurm/
