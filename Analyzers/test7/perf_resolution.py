#!/usr/bin/env python

hname2023_f = lambda hname: "emtf2023_" + hname[5:]

donotdelete = []



# ______________________________________________________________________________
if __name__ == '__main__':
  from ROOT import gROOT, gPad, gStyle, TFile, TCanvas, TH1F, TH2F, TPolyLine, TLatex, TColor, TEfficiency, TGraphAsymmErrors

  # ____________________________________________________________________________
  # Setup basic drawer
  gROOT.LoadMacro("tdrstyle.C")
  gROOT.ProcessLine("setTDRStyle();")
  #gROOT.SetBatch(True)
  gStyle.SetPadRightMargin(0.05)
  gStyle.SetPadGridX(True)
  gStyle.SetPadGridY(True)
  gStyle.SetTitleOffset(0.9, "Y")
  gStyle.SetTitleOffset(1.2, "Y")
  gROOT.ForceStyle()

  infile = "histos_tbc_add.18.root"
  tfile = TFile.Open(infile)


  # ____________________________________________________________________________
  # emtf_l1pt_vs_genpt
  hname = "emtf_l1pt_vs_genpt"
  h2a = tfile.Get(hname)
  h2a.Draw("COLZ")
  gPad.Print("figures_winter/" + hname + ".png")

  hname = hname2023_f(hname)
  h2b = tfile.Get(hname)
  h2b.Draw("COLZ")
  gPad.Print("figures_winter/" + hname + ".png")

  donotdelete.append([h2a, h2b])

  # ____________________________________________________________________________
  # emtf_l1ptres_vs_genpt
  def doit():
    frame = h.ProfileX(hname+"_frame", 1, -1, "s")
    gr1 = TGraphAsymmErrors(h.GetNbinsX())
    gr2 = TGraphAsymmErrors(h.GetNbinsX())
    gr1_aspt = TGraphAsymmErrors(h.GetNbinsX())
    gr2_aspt = TGraphAsymmErrors(h.GetNbinsX())
    # Apply gaussian fits
    for i in xrange(h.GetNbinsX()):
      h_py = h.ProjectionY("_py", i+1, i+1)
      h_py.Rebin(3)  # 300 bins -> 100
      if h_py.Integral() < 20:  continue
      r = h_py.Fit("gaus", "SNQ", "", -1, 1.2)
      #r = h_py.Fit("gaus", "SNQ", "", h_py.GetMean() - 0.04*5, h_py.GetMean() + 0.04*5)
      mean, sigma, meanErr, sigmaErr = r.Parameter(1), r.Parameter(2), r.ParError(1), r.ParError(2)
      gr1.SetPoint(i, h.GetXaxis().GetBinCenter(i+1), mean)
      gr1.SetPointError(i, 0, 0, sigma, sigma)
      gr2.SetPoint(i, h.GetXaxis().GetBinCenter(i+1), sigma)
      gr2.SetPointError(i, 0, 0, sigmaErr, sigmaErr)
      gr1_aspt.SetPoint(i, 1.0/h.GetXaxis().GetBinCenter(i+1), mean)
      gr1_aspt.SetPointError(i, 0, 0, sigma, sigma)
      gr2_aspt.SetPoint(i, 1.0/h.GetXaxis().GetBinCenter(i+1), sigma)
      gr2_aspt.SetPointError(i, 0, 0, sigmaErr, sigmaErr)
    # Draw
    h.Draw("COLZ")
    gPad.SetLogx(0)
    gPad.Print("figures_winter/" + hname + ".png")
    #
    frame.Reset()
    frame.SetBins(50, 0, 50)
    frame.GetXaxis().SetTitle("gen p_{T} [GeV]")
    frame.GetYaxis().SetTitle("#Delta(p_{T})/p_{T} bias")
    frame.SetMaximum(0.5)
    frame.SetMinimum(-0.5)
    frame.SetStats(0)
    frame.Draw()
    gr1_aspt.SetLineColor(col)
    gr1_aspt.SetMarkerColor(col)
    gr1_aspt.Draw("p")
    gPad.SetLogx()
    gPad.Print("figures_winter/" + hname + "_bias.png")
    #
    frame.GetXaxis().SetTitle("gen p_{T} [GeV]")
    frame.GetYaxis().SetTitle("#Delta(p_{T})/p_{T} resolution")
    frame.SetMaximum(0.6)
    frame.SetMinimum(0.0)
    frame.SetStats(0)
    frame.Draw()
    gr2_aspt.SetLineColor(col)
    gr2_aspt.SetMarkerColor(col)
    gr2_aspt.Draw("p")
    #gr2_aspt.Fit("pol1", "", "", 10, 40)
    gPad.SetLogx()
    gPad.Print("figures_winter/" + hname + "_res.png")
    #
    h.cache = [frame, gr1, gr2, gr1_aspt, gr2_aspt]

  hname = "emtf_l1ptres_vs_genpt"
  h = tfile.Get(hname)
  col = 632  # kRed
  doit()

  hname = hname2023_f(hname)
  h = tfile.Get(hname)
  col = 600  # kBlue
  doit()

  donotdelete.append([h2a, h2b])

