"""
Usage: python -i pattern.py --pt 17.6 --phi 43 43 43 43 --imgname patt1
"""


# ______________________________________________________________________________
if __name__ == '__main__':
  from ROOT import gROOT, gPad, gStyle, TCanvas, TH1F, TH2F, TPolyLine, TLatex
  from array import array
  import argparse

  # ____________________________________________________________________________
  # Setup parser
  parser = argparse.ArgumentParser()
  parser.add_argument("--phi", type=int, nargs='*', help="Input zone phis (default: %(default)s)")
  parser.add_argument("--pt", type=float, help="Input pt (default: %(default)s)")
  parser.add_argument("--imgname", type=str, help="Image name (default: %(default)s)")
  options = parser.parse_args()
  assert len(options.phi) == 4

  # ____________________________________________________________________________
  # Setup basic drawer
  gROOT.LoadMacro("tdrstyle.C")
  gROOT.ProcessLine("setTDRStyle();")
  gROOT.SetBatch(True)
  gStyle.SetTickLength(0., "XYZ")
  gStyle.SetLabelSize(0., "XYZ")
  gStyle.SetPadTopMargin(0.05)
  gStyle.SetPadBottomMargin(0.05)
  gStyle.SetPadLeftMargin(0.05)
  gStyle.SetPadRightMargin(0.05)
  gStyle.SetFrameLineColor(0)

  # Frame
  nbinsx, xlow, xhi = 32, -1, 31
  nbinsy, ylow, yhi = 5, 0, 5
  frame = TH2F("frame", "", nbinsx, xlow, xhi, nbinsy, ylow, yhi)

  # Empty boxes
  tpolylines = []
  for lay in xrange(4):
    for strip in xrange(31):
      x = [strip - 0.5, strip - 0.5, strip + 0.5, strip + 0.5, strip - 0.5]
      y = [(4-lay) - 0.5, (4-lay) + 0.5, (4-lay) + 0.5, (4-lay) - 0.5, (4-lay) - 0.5]
      tpolyline = TPolyLine(len(x), array('d', x), array('d', y))
      tpolylines.append(tpolyline)

  # Find 'key' zone_phi
  key_zone_phi = -1
  min_distance = 999
  for lay in range(1,4):
    distance = 0
    for zone_phi in options.phi[1:4]:
      if zone_phi != -1:
        distance += abs(zone_phi - options.phi[lay])
    if min_distance > distance:
      min_distance = distance
      key_zone_phi = options.phi[lay]

  # Red boxes
  tpolylines_1 = []
  for zone_phi, lay in zip(options.phi, range(4)):
    strip = zone_phi - key_zone_phi + 15
    x = [strip - 0.5, strip - 0.5, strip + 0.5, strip + 0.5, strip - 0.5]
    y = [(4-lay) - 0.5, (4-lay) + 0.5, (4-lay) + 0.5, (4-lay) - 0.5, (4-lay) - 0.5]
    tpolyline = TPolyLine(len(x), array('d', x), array('d', y))
    tpolylines_1.append(tpolyline)

  #strips = TH2F("strips", "", nbinsx, xlow+0.5, xhi+0.5, nbinsy, ylow+0.5, yhi+0.5)
  #for zone_phi, lay in zip(options.phi, range(4)):
  #  strip = zone_phi - key_zone_phi + 15
  #  strips.Fill(strip, (4-lay))

  # Draw
  c1 = TCanvas("c1", "c1", 900, 150)
  frame.SetStats(0)
  frame.Draw()

  for l in tpolylines:
    l.Draw()

  for l in tpolylines_1:
    l.SetFillColor(632)  # kRed
    l.Draw("f")
    l.Draw()

  #strips.SetFillColor(632)  # kRed
  #strips.Draw("same col")

  # Text
  tlatex = TLatex()
  tlatex.SetNDC()
  tlatex.SetTextFont(42)
  tlatex.SetTextSize(0.2)
  text = "#color[632]{p_{T} = %.2f}" % options.pt
  tlatex.DrawLatex(0.8, 0.171, text)

  c1.Print("plots_pattern/"+options.imgname+".png")
