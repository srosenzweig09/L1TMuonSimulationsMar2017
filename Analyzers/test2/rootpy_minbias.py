import numpy as np
np.random.seed(2023)

from rootpy.plotting import Hist, Hist2D, Graph, Efficiency, Legend, Canvas
from rootpy.tree import Tree, TreeModel, FloatCol, IntCol, ShortCol
from rootpy.io import root_open

from ROOT import TF1


# ______________________________________________________________________________
infile = root_open('minbiasmuonanalyzer_all.root')
h2 = infile["minbiasmuonanalyzer/muon_invPt_vs_eta"]

nevents = infile.minbiasmuonanalyzer.nevents
nevents = nevents.GetBinContent(2)

pt_rebin = 4
eta_rebin = 2


make_plots = True

if make_plots:
  from drawer import *
  mydrawer = MyDrawer()
  options = mydrawer.options

  bookkeeping = []
  kBlue = 600

  #gROOT.SetBatch(0)


  # ____________________________________________________________________________
  if False:
    h1 = infile["minbiasmuonanalyzer/muon_ptmin2_absEtaMin1.24_absEtaMax2.5_dxy"]
    h1.Sumw2()
    h1.SetMarkerSize(0)
    h1.SetLineWidth(2)
    h1.SetLineColor(kBlue)
    h1.GetXaxis().SetRangeUser(0,50)
    h1.Draw()
    gPad.SetLogy()
    gPad.Print("minbias_muon_dxy.png")

    h1 = infile["minbiasmuonanalyzer/muon_ptmin2_absEtaMin1.24_absEtaMax2.5_dz"]
    h1.Sumw2()
    h1.SetMarkerSize(0)
    h1.SetLineWidth(2)
    h1.SetLineColor(kBlue)
    h1.GetXaxis().SetRangeUser(0,50)
    h1.Draw()
    gPad.SetLogy()
    gPad.Print("minbias_muon_dz.png")

    h1 = infile["minbiasmuonanalyzer/muon_absEtaMin1.24_absEtaMax2.5_pt"]
    h1.Sumw2()
    h1.SetMarkerSize(0)
    h1.SetLineWidth(2)
    h1.SetLineColor(kBlue)
    h1.GetXaxis().SetRangeUser(0,50)
    h1.Scale(1.0/h1.Integral())
    h1.Draw()
    gPad.SetLogy()
    gPad.Print("minbias_muon_pt.png")

    h1 = infile["minbiasmuonanalyzer/muon_absEtaMin1.24_absEtaMax2.5_logPt"]
    h1.Sumw2()
    h1.SetMarkerSize(0)
    h1.SetLineWidth(2)
    h1.SetLineColor(kBlue)
    h1.GetXaxis().SetRangeUser(0,6)
    h1.Scale(1.0/h1.Integral())
    h1.Draw()
    gPad.SetLogy()
    gPad.Print("minbias_muon_logPt.png")

    h1 = infile["minbiasmuonanalyzer/muon_absEtaMin1.24_absEtaMax2.5_invPt"]
    h1.Sumw2()
    h1.SetMarkerSize(0)
    h1.SetLineWidth(2)
    h1.SetLineColor(kBlue)
    h1.Rebin(2)
    h1.Scale(1.0/h1.Integral())
    h1.Draw()
    gPad.SetLogy(0)
    gPad.Print("minbias_muon_invPt.png")

    h1 = infile["minbiasmuonanalyzer/muon_absEtaMin1.24_absEtaMax2.5_invPt2"]
    h1.Sumw2()
    h1.SetMarkerSize(0)
    h1.SetLineWidth(2)
    h1.SetLineColor(kBlue)
    h1.Rebin(2)
    h1.Scale(1.0/h1.Integral())
    h1.Draw()
    gPad.SetLogy(0)
    gPad.Print("minbias_muon_invPt2.png")

    h1 = infile["minbiasmuonanalyzer/muon_absEtaMin1.24_absEtaMax2.5_invPt3"]
    h1.Sumw2()
    h1.SetMarkerSize(0)
    h1.SetLineWidth(2)
    h1.SetLineColor(kBlue)
    h1.Rebin(2)
    h1.Scale(1.0/h1.Integral())
    h1.Draw()
    gPad.SetLogy(0)
    gPad.Print("minbias_muon_invPt3.png")

    h1 = infile["minbiasmuonanalyzer/muon_absEtaMin1.24_absEtaMax2.5_invPt4"]
    h1.Sumw2()
    h1.SetMarkerSize(0)
    h1.SetLineWidth(2)
    h1.SetLineColor(kBlue)
    h1.Rebin(2)
    h1.Scale(1.0/h1.Integral())
    h1.Draw()
    gPad.SetLogy(0)
    gPad.Print("minbias_muon_invPt4.png")

    h1 = infile["minbiasmuonanalyzer/muon_absEtaMin1.24_absEtaMax2.5_invPt5"]
    h1.Sumw2()
    h1.SetMarkerSize(0)
    h1.SetLineWidth(2)
    h1.SetLineColor(kBlue)
    h1.Rebin(2)
    h1.Scale(1.0/h1.Integral())
    h1.Draw()
    gPad.SetLogy(0)
    gPad.Print("minbias_muon_invPt5.png")


  # ____________________________________________________________________________
  for i in xrange(13):
    eta = 1.2 + 0.1 * i
    b = h2.GetXaxis().FindBin(eta - 0.001)
    h_py = h2.ProjectionY("_py", b, b+(eta_rebin-1), "")
    h_py.Rebin(pt_rebin)
    h_py.Sumw2()
    h_py.Scale(1.0/nevents/eta_rebin/pt_rebin)
    h_py.SetMaximum(1.6e-6)
    h_py.SetMinimum(0.)

    fa1 = TF1("fa1", "[0]*x + [1]*x*x + [2]*x*x*x")
    h_py.Fit("fa1", "IM", "", 0.04, 0.5)
    h_py.Draw()
    if (1.899 < eta < 1.901):
      bookkeeping.append(fa1.GetParameter(0))
      bookkeeping.append(fa1.GetParameter(1))
      bookkeeping.append(fa1.GetParameter(2))
      bookkeeping.append(h_py.Integral(h_py.FindBin(0), h_py.FindBin(0.3)))

    gPad.Print("minbias_muon_pt_eta%.1f.png" % eta)

  print bookkeeping
