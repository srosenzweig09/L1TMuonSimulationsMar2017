import numpy as np

input_list = []

with open("log4", "r") as f:
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
  for i in xrange(6):
    h_strips = TH2F("strips_%i" % i, "", nbinsx, xlow+0.5, xhi+0.5, nbinsy, ylow+0.5, yhi+0.5)
    h_strips_cov90 = TH2F("strips_cov90_%i" % i, "", nbinsx, xlow+0.5, xhi+0.5, nbinsy, ylow+0.5, yhi+0.5)
    s1_strips = []
    s2_strips = []
    s3_strips = []
    s4_strips = []

    options.pt_list = []

    if i == 0:
      pt_select = lambda x: 20 <= x < 100
    elif i == 1:
      pt_select = lambda x: 10 <= x < 20
    elif i == 2:
      pt_select = lambda x: 5 <= x < 10
    elif i == 3:
      pt_select = lambda x: 3 <= x < 5
    elif i == 4:
      pt_select = lambda x: 2 <= x < 3
    elif i == 5:
      pt_select = lambda x: 1 <= x < 2


    for ievt in xrange(len(input_list)):
      options.pt = input_list[ievt][0]
      options.q = input_list[ievt][1]
      options.phi = input_list[ievt][2]

      ## Find 'key' zone_phi
      #key_zone_phi = -1
      #min_distance = 999
      #for lay in range(1,4):  # exclude station 1
      #  distance = 0
      #  for zone_phi in options.phi[1:4]:  # exclude station 1
      #    if zone_phi != -1:
      #      distance += abs(zone_phi - options.phi[lay])
      #  if min_distance > distance:
      #    min_distance = distance
      #    key_zone_phi = options.phi[lay]

      ## Find 'key' zone_phi
      #key_zone_phi = -1
      #zone_phi = options.phi[1]  # use station 2
      ##zone_phi = options.phi[2]  # use station 3
      ##zone_phi = options.phi[3]  # use station 4
      #if zone_phi != -1:
      #  key_zone_phi = zone_phi

      # Find 'key' zone_phi
      key_zone_phi = -1
      if options.phi[1] != -1:
        key_zone_phi = options.phi[1]
      elif options.phi[2] != -1:
        key_zone_phi = options.phi[2]
      elif options.phi[3] != -1:
        key_zone_phi = options.phi[3]

      ## Find station 1 zone_phi
      #s1_zone_phi = -1
      #zone_phi = options.phi[0]
      #if zone_phi != -1:
      #  s1_zone_phi = zone_phi - key_zone_phi + 15

      # Fill histogram
      if key_zone_phi != -1 and pt_select(options.pt):
        for zone_phi, lay in zip(options.phi, range(4)):
          if zone_phi != -1:
            strip = zone_phi - key_zone_phi + 15
            if options.q == 1:
              strip = -strip + 30
            if i == 0: # debug
              print ievt, lay, strip
            h_strips.Fill(strip, (4-lay))
            if lay == 0:
              s1_strips.append(strip)
            elif lay == 1:
              s2_strips.append(strip)
            elif lay == 2:
              s3_strips.append(strip)
            elif lay == 3:
              s4_strips.append(strip)
        options.pt_list.append(options.pt)
      continue  # end ievt loop

    # __________________________________________________________________________
    s1_strips_a, s1_strips_b = np.percentile(s1_strips, [5,95])
    s2_strips_a, s2_strips_b = np.percentile(s2_strips, [5,95])
    s3_strips_a, s3_strips_b = np.percentile(s3_strips, [5,95])
    s4_strips_a, s4_strips_b = np.percentile(s4_strips, [5,95])

    s1_strips_a, s1_strips_b = -999, 999
    #s2_strips_a, s2_strips_b = -999, 999  # use station 2
    #s3_strips_a, s3_strips_b = -999, 999  # use station 3
    #s4_strips_a, s4_strips_b = -999, 999  # use station 4

    for strip in s1_strips:
      if s1_strips_a <= strip <= s1_strips_b:
        lay = 0
        h_strips_cov90.Fill(strip, (4-lay))
    for strip in s2_strips:
      if s2_strips_a <= strip <= s2_strips_b:
        lay = 1
        h_strips_cov90.Fill(strip, (4-lay))
    for strip in s3_strips:
      if s3_strips_a <= strip <= s3_strips_b:
        lay = 2
        h_strips_cov90.Fill(strip, (4-lay))
    for strip in s4_strips:
      if s4_strips_a <= strip <= s4_strips_b:
        lay = 3
        h_strips_cov90.Fill(strip, (4-lay))

    # Draw "same"
    frame.Draw()
    h_strips_cov90.Draw("same col")
    for l in tpolylines:
      l.Draw()

    # Text
    tlatex = TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(42)
    tlatex.SetTextSize(0.2)
    text = "#color[632]{<p_{T}> = %.2f, N=%i}" % (sum(options.pt_list)/float(len(options.pt_list)), len(options.pt_list))
    tlatex.DrawLatex(0.65, 0.171, text)

    if i == 15:  # debug
      h_strips.Print("all")

    options.imgname = "patt%i" % i
    c1.Print("plots_pattern_merge/"+options.imgname+".png")
    continue  # end i loop

  #scp -r -C plots_pattern_merge/ jlow@melrose.ihepa.ufl.edu:/home/jlow/L1MuonTrigger/P2_CMSSW_9_2_3_patch1/src/L1TMuonSimulations/Analyzers/test4/
