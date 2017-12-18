input_list = []

with open("log3", "r") as f:
  for line in f:
    line = line.strip()
    t = eval(line)
    input_list.append(t)

assert len(input_list) > 1000


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
  #assert len(options.phi) == 4

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
  #gStyle.SetNumberContours(100)
  gStyle.SetNumberContours(50)
  #gStyle.SetPalette(55)  # kRainbow
  #gStyle.SetPalette(57)  # kBird

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

  # Draw
  c1 = TCanvas("c1", "c1", 900, 150)
  frame.SetStats(0)
  frame.Draw()

  # Loop over events
  ievt = 0
  for i in xrange(100):
    strips = TH2F("strips_%i" % i, "", nbinsx, xlow+0.5, xhi+0.5, nbinsy, ylow+0.5, yhi+0.5)

    options.pt = input_list[ievt][0]
    options.pt_list = []

    print "[INFO] ievt %i range %i pt %f" % (ievt, (50 * (1<<i)), options.pt)

    xrange_i = xrange(50 * (1<<i))

    for j in xrange_i:
      if ievt == len(input_list):
        break

      options.pt = input_list[ievt][0]
      options.phi = input_list[ievt][1]

      options.pt_list.append(options.pt)

      # Find 'key' zone_phi
      key_zone_phi = -1
      min_distance = 999
      for lay in range(1,4):  # exclude station 1
        distance = 0
        for zone_phi in options.phi[1:4]:  # exclude station 1
          if zone_phi != -1:
            distance += abs(zone_phi - options.phi[lay])
        if min_distance > distance:
          min_distance = distance
          key_zone_phi = options.phi[lay]

      for zone_phi, lay in zip(options.phi, range(4)):
        if zone_phi != -1:
          strip = zone_phi - key_zone_phi + 15
          strips.Fill(strip, (4-lay))

      ievt += 1
      continue  # end j loop

    if ievt == len(input_list):
      break

    # Draw "same"
    frame.Draw()
    strips.Draw("same col")
    for l in tpolylines:
      l.Draw()

    # Text
    tlatex = TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(42)
    tlatex.SetTextSize(0.2)
    text = "#color[632]{<p_{T}> = %.2f, N=%i}" % (sum(options.pt_list)/float(len(options.pt_list)), len(options.pt_list))
    tlatex.DrawLatex(0.65, 0.171, text)

    options.imgname = "patt%i" % (i+1)
    c1.Print("plots_pattern_merge/"+options.imgname+".png")
    continue  # end i loop

  #scp -r -C plots_pattern_merge/ jlow@melrose.ihepa.ufl.edu:/home/jlow/L1MuonTrigger/P2_CMSSW_9_2_3_patch1/src/L1TMuonSimulations/Analyzers/test4/
