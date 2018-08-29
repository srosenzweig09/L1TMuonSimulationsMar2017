#!/usr/bin/env python

hname2023_f = lambda hname: "emtf2023_" + hname[5:]

donotdelete = []



# ______________________________________________________________________________
if __name__ == '__main__':
  from ROOT import gROOT, gPad, gStyle, TFile, TCanvas, TH1F, TH2F, TPolyLine, TLatex, TColor, TEfficiency, TLine

  # ____________________________________________________________________________
  # Setup basic drawer
  gROOT.LoadMacro("tdrstyle.C")
  gROOT.ProcessLine("setTDRStyle();")
  #gROOT.SetBatch(True)
  gStyle.SetMarkerStyle(1)
  gStyle.SetEndErrorSize(0)
  gStyle.SetPadGridX(True)
  gStyle.SetPadGridY(True)
  gROOT.ForceStyle()

  infile = "histos_tbc_add.17.root"
  tfile = TFile.Open(infile)

  tline = TLine()
  tline.SetLineColor(920+2)  # kGray+2
  tline.SetLineStyle(2)


  # ____________________________________________________________________________
  # emtf_eff_vs_genpt_l1pt20
  hname = "emtf_eff_vs_genpt_l1pt20"
  h1a_denom = tfile.Get(hname + "_denom")
  h1a_numer = tfile.Get(hname + "_numer")
  h1b_denom = tfile.Get(hname2023_f(hname) + "_denom")
  h1b_numer = tfile.Get(hname2023_f(hname) + "_numer")

  print h1a_denom.GetEntries(), h1a_numer.GetEntries()
  print h1b_denom.GetEntries(), h1b_denom.GetEntries()

  h1a_eff = TEfficiency(h1a_numer, h1a_denom)
  h1a_eff.SetStatisticOption(0)  # kFCP
  h1a_eff.SetConfidenceLevel(0.682689492137)  # one sigma
  h1a_eff.SetMarkerColor(632)  # kRed
  h1a_eff.SetLineColor(632)  # kRed
  h1a_eff.SetLineWidth(2)

  h1b_eff = TEfficiency(h1b_numer, h1b_denom)
  h1b_eff.SetStatisticOption(0)  # kFCP
  h1b_eff.SetConfidenceLevel(0.682689492137)  # one sigma
  h1b_eff.SetMarkerColor(600)  # kBlue
  h1b_eff.SetLineColor(600)  # kBlue
  h1b_eff.SetLineWidth(2)

  gr = h1a_eff.CreateGraph()
  gr.Draw("ap")
  frame = gr.GetHistogram()
  frame = frame.Clone(hname + "_frame")
  frame.GetYaxis().SetTitle("#varepsilon")
  frame.SetMinimum(0.0)
  frame.SetMaximum(1.2)
  frame.Draw()
  xmin, xmax = frame.GetXaxis().GetXmin(), frame.GetXaxis().GetXmax()
  tline.DrawLine(xmin, 1.0, xmax, 1.0)
  h1a_eff.Draw("same")
  h1b_eff.Draw("same")

  gPad.Print("figures_winter/" + hname + ".png")
  donotdelete.append([frame, h1a_eff, h1b_eff])

  if True:
    outfile = TFile.Open("figures_winter/" + hname + ".root", "RECREATE")
    for obj in [frame, h1a_eff, h1b_eff]:
      obj.Write()
    outfile.Close()


  # ____________________________________________________________________________
  # emtf_eff_vs_geneta_l1pt20
  hname = "emtf_eff_vs_geneta_l1pt20"
  h1a_denom = tfile.Get(hname + "_denom")
  h1a_numer = tfile.Get(hname + "_numer")
  h1b_denom = tfile.Get(hname2023_f(hname) + "_denom")
  h1b_numer = tfile.Get(hname2023_f(hname) + "_numer")

  print h1a_denom.GetEntries(), h1a_numer.GetEntries()
  print h1b_denom.GetEntries(), h1b_denom.GetEntries()

  h1a_eff = TEfficiency(h1a_numer, h1a_denom)
  h1a_eff.SetStatisticOption(0)  # kFCP
  h1a_eff.SetConfidenceLevel(0.682689492137)  # one sigma
  h1a_eff.SetMarkerColor(632)  # kRed
  h1a_eff.SetLineColor(632)  # kRed
  h1a_eff.SetLineWidth(2)

  h1b_eff = TEfficiency(h1b_numer, h1b_denom)
  h1b_eff.SetStatisticOption(0)  # kFCP
  h1b_eff.SetConfidenceLevel(0.682689492137)  # one sigma
  h1b_eff.SetMarkerColor(600)  # kBlue
  h1b_eff.SetLineColor(600)  # kBlue
  h1b_eff.SetLineWidth(2)

  gr = h1a_eff.CreateGraph()
  gr.Draw("ap")
  frame = gr.GetHistogram()
  frame = frame.Clone(hname + "_frame")
  frame.GetYaxis().SetTitle("#varepsilon")
  frame.SetMinimum(0.0)
  frame.SetMaximum(1.2)
  frame.Draw()
  xmin, xmax = frame.GetXaxis().GetXmin(), frame.GetXaxis().GetXmax()
  tline.DrawLine(xmin, 1.0, xmax, 1.0)
  h1a_eff.Draw("same")
  h1b_eff.Draw("same")

  gPad.Print("figures_winter/" + hname + ".png")
  donotdelete.append([frame, h1a_eff, h1b_eff])


  # ____________________________________________________________________________
  # emtf_eff_vs_geneta_genpt30_l1pt20
  hname = "emtf_eff_vs_geneta_genpt30_l1pt20"
  h1a_denom = tfile.Get(hname + "_denom")
  h1a_numer = tfile.Get(hname + "_numer")
  h1b_denom = tfile.Get(hname2023_f(hname) + "_denom")
  h1b_numer = tfile.Get(hname2023_f(hname) + "_numer")

  print h1a_denom.GetEntries(), h1a_numer.GetEntries()
  print h1b_denom.GetEntries(), h1b_denom.GetEntries()

  h1a_eff = TEfficiency(h1a_numer, h1a_denom)
  h1a_eff.SetStatisticOption(0)  # kFCP
  h1a_eff.SetConfidenceLevel(0.682689492137)  # one sigma
  h1a_eff.SetMarkerColor(632)  # kRed
  h1a_eff.SetLineColor(632)  # kRed
  h1a_eff.SetLineWidth(2)

  h1b_eff = TEfficiency(h1b_numer, h1b_denom)
  h1b_eff.SetStatisticOption(0)  # kFCP
  h1b_eff.SetConfidenceLevel(0.682689492137)  # one sigma
  h1b_eff.SetMarkerColor(600)  # kBlue
  h1b_eff.SetLineColor(600)  # kBlue
  h1b_eff.SetLineWidth(2)

  gr = h1a_eff.CreateGraph()
  gr.Draw("ap")
  frame = gr.GetHistogram()
  frame = frame.Clone(hname + "_frame")
  frame.GetYaxis().SetTitle("#varepsilon")
  frame.SetMinimum(0.0)
  frame.SetMaximum(1.2)
  frame.Draw()
  xmin, xmax = frame.GetXaxis().GetXmin(), frame.GetXaxis().GetXmax()
  tline.DrawLine(xmin, 1.0, xmax, 1.0)
  h1a_eff.Draw("same")
  h1b_eff.Draw("same")

  gPad.Print("figures_winter/" + hname + ".png")
  donotdelete.append([frame, h1a_eff, h1b_eff])


  # ____________________________________________________________________________
  # emtf_eff_vs_genpt_l1pt99
  palette = ("#333333", "#377eb8", "#e41a1c", "#984ea3", "#ff7f00", "#4daf4a")
  palette = map(lambda x: TColor.GetColor(x), palette)
  hname = "emtf_eff_vs_genpt_l1pt%i"

  i = 0
  effs = []

  for l in (0, 10, 20, 30, 40, 50):
    denom = tfile.Get((hname % (l)) + "_denom")
    numer = tfile.Get((hname % (l)) + "_numer")
    eff = TEfficiency(numer, denom)
    eff.SetStatisticOption(0)  # kFCP
    eff.SetConfidenceLevel(0.682689492137)  # one sigma
    eff.SetMarkerColor(palette[i])  # kRed
    eff.SetLineColor(palette[i])  # kRed
    eff.SetLineWidth(2)
    effs.append(eff)

    if l == 0:
      gr = eff.CreateGraph()
      gr.Draw("ap")
      frame = gr.GetHistogram()
      frame = frame.Clone(hname + "_frame")
      frame.GetYaxis().SetTitle("#varepsilon")
      frame.SetMinimum(0.0)
      frame.SetMaximum(1.2)
      frame.Draw()
      xmin, xmax = frame.GetXaxis().GetXmin(), frame.GetXaxis().GetXmax()
      tline.DrawLine(xmin, 1.0, xmax, 1.0)
      #eff.Draw("same")
    else:
      eff.Draw("same")
    if l == 50:
      gPad.Print("figures_winter/" + (hname % (99)) + ".png")

    i += 1


  # ____________________________________________________________________________
  # emtf2023_eff_vs_genpt_l1pt99
  palette = ("#333333", "#377eb8", "#e41a1c", "#984ea3", "#ff7f00", "#4daf4a")
  palette = map(lambda x: TColor.GetColor(x), palette)
  hname = "emtf2023_eff_vs_genpt_l1pt%i"

  i = 0
  effs = []

  for l in (0, 10, 20, 30, 40, 50):
    denom = tfile.Get((hname % (l)) + "_denom")
    numer = tfile.Get((hname % (l)) + "_numer")
    eff = TEfficiency(numer, denom)
    eff.SetStatisticOption(0)  # kFCP
    eff.SetConfidenceLevel(0.682689492137)  # one sigma
    eff.SetMarkerColor(palette[i])  # kRed
    eff.SetLineColor(palette[i])  # kRed
    eff.SetLineWidth(2)
    effs.append(eff)

    if l == 0:
      gr = eff.CreateGraph()
      gr.Draw("ap")
      frame = gr.GetHistogram()
      frame = frame.Clone(hname + "_frame")
      frame.GetYaxis().SetTitle("#varepsilon")
      frame.SetMinimum(0.0)
      frame.SetMaximum(1.2)
      frame.Draw()
      xmin, xmax = frame.GetXaxis().GetXmin(), frame.GetXaxis().GetXmax()
      tline.DrawLine(xmin, 1.0, xmax, 1.0)
      #eff.Draw("same")
    elif l == 10:
      pass
    else:
      #if l == 20: eff.CreateGraph().Print("all")
      eff.Draw("same")
    if l == 50:
      gPad.Print("figures_winter/" + (hname % (99)) + ".png")

    i += 1
