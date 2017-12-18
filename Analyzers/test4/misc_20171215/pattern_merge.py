import numpy as np

input_list = []
histos = []

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
  c2 = TCanvas("c2", "c2", 600, 600)
  c2.SetTopMargin(0.05)
  c2.SetBottomMargin(0.13)
  c2.SetLeftMargin(0.15)
  c2.SetRightMargin(0.03)
  frame.SetStats(0)
  frame.Draw()

  # Loop over events
  for i in xrange(9):
    h_strips = TH2F("strips_%i" % i, "", nbinsx, xlow+0.5, xhi+0.5, nbinsy, ylow+0.5, yhi+0.5)
    h_strips_cov90 = TH2F("strips_cov90_%i" % i, "", nbinsx, xlow+0.5, xhi+0.5, nbinsy, ylow+0.5, yhi+0.5)
    h_pt_vs_s1 = TH2F("pt_vs_strips_%i_s1" % i, "; s1 strip; p_{T} [GeV]", nbinsx, xlow+0.5, xhi+0.5, 40, 0., 20.)
    h_pt_vs_s2 = TH2F("pt_vs_strips_%i_s2" % i, "; s2 strip; p_{T} [GeV]", nbinsx, xlow+0.5, xhi+0.5, 40, 0., 20.)
    h_pt_vs_s3 = TH2F("pt_vs_strips_%i_s3" % i, "; s3 strip; p_{T} [GeV]", nbinsx, xlow+0.5, xhi+0.5, 40, 0., 20.)
    h_pt_vs_s4 = TH2F("pt_vs_strips_%i_s4" % i, "; s4 strip; p_{T} [GeV]", nbinsx, xlow+0.5, xhi+0.5, 40, 0., 20.)
    s1_strips = []
    s2_strips = []
    s3_strips = []
    s4_strips = []

    options.pt_list = []

    #s1_select = lambda x: x == i  # dumb version

    if i == 0:
      s1_select = lambda x:  (x - 15) == 0
    elif i == 1:
      s1_select = lambda x:  -(x - 15) == 1
    elif i == 2:
      s1_select = lambda x:  2 <= -(x - 15) <= 3
    elif i == 3:
      s1_select = lambda x:  4 <= -(x - 15) <= 6
    elif i == 4:
      s1_select = lambda x:  7 <= -(x - 15) <= 10
    elif i == 5:
      s1_select = lambda x:  (x - 15) == 1
    elif i == 6:
      s1_select = lambda x:  2 <= (x - 15) <= 3
    elif i == 7:
      s1_select = lambda x:  4 <= (x - 15) <= 6
    elif i == 8:
      s1_select = lambda x:  7 <= (x - 15) <= 10

    if i == 0:
      part_select = lambda pt, q: pt > 10
    elif i == 1:
      part_select = lambda pt, q: pt > 6 and q == -1
    elif i == 2:
      part_select = lambda pt, q: pt > 4 and q == -1
    elif i == 3:
      part_select = lambda pt, q: pt > 3 and q == -1
    elif i == 4:
      part_select = lambda pt, q: pt > 2.5 and q == -1
    elif i == 5:
      part_select = lambda pt, q: pt > 6 and q == 1
    elif i == 6:
      part_select = lambda pt, q: pt > 4 and q == 1
    elif i == 7:
      part_select = lambda pt, q: pt > 3 and q == 1
    elif i == 8:
      part_select = lambda pt, q: pt > 2.5 and q == 1


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

      # Find station 1 zone_phi
      s1_zone_phi = -1
      zone_phi = options.phi[0]
      if zone_phi != -1:
        s1_zone_phi = zone_phi - key_zone_phi + 15

      # Fill histogram
      if key_zone_phi != -1 and s1_select(s1_zone_phi) and part_select(options.pt, options.q):
        for zone_phi, lay in zip(options.phi, range(4)):
          if zone_phi != -1:
            strip = zone_phi - key_zone_phi + 15
            if i == 0: # debug
              print ievt, lay, strip
            h_strips.Fill(strip, (4-lay))
            if lay == 0:
              s1_strips.append(strip)
              h_pt_vs_s1.Fill(strip, options.pt)
            elif lay == 1:
              s2_strips.append(strip)
              h_pt_vs_s2.Fill(strip, options.pt)
            elif lay == 2:
              s3_strips.append(strip)
              h_pt_vs_s3.Fill(strip, options.pt)
            elif lay == 3:
              s4_strips.append(strip)
              h_pt_vs_s4.Fill(strip, options.pt)
        options.pt_list.append(options.pt)
      continue  # end ievt loop

    # __________________________________________________________________________
    s1_strips_a, s1_strips_b = np.percentile(s1_strips, [12.5,87.5])
    s2_strips_a, s2_strips_b = np.percentile(s2_strips, [12.5,87.5])
    s3_strips_a, s3_strips_b = np.percentile(s3_strips, [12.5,87.5])
    s4_strips_a, s4_strips_b = np.percentile(s4_strips, [12.5,87.5])

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

    # Draw pT vs strips
    if False:
      c2.cd()
      for h in [h_pt_vs_s1, h_pt_vs_s2, h_pt_vs_s3, h_pt_vs_s4]:
        h.SetTickLength(0.03, "XYZ")
        h.SetLabelSize(0.05, "XYZ")
        h.Draw("col")
        c2.Print("plots_pattern_merge/" + h.GetName() + ".png")

    # Draw "same"
    c1.cd()
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


  # ____________________________________________________________________________
  if True:
    # Finalize
    cfi_strings = [""] * 9

    for i in xrange(9):
      h_strips_new = TH2F("strips_new_%i" % i, "", nbinsx, xlow+0.5, xhi+0.5, nbinsy, ylow+0.5, yhi+0.5)

      if i == 0:
        s1_strips = [15]
        s2_strips = [15]
        s3_strips = [15]
        s4_strips = [15]
      elif i == 1 or i == 5:
        s1_strips = [14]
        s2_strips = [15]
        s3_strips = [15,16]
        s4_strips = [14,15,16]
      elif i == 2 or i == 6:
        s1_strips = [12,13]
        s2_strips = [15]
        s3_strips = [14,15,16]
        s4_strips = [13,14,15,16]
      elif i == 3 or i == 7:
        s1_strips = [9,10,11]
        s2_strips = [15]
        s3_strips = [14,15,16]
        s4_strips = [12,13,14,15,16,17]
      elif i == 4 or i == 8:
        s1_strips = [5,6,7,8]
        s2_strips = [15]
        s3_strips = [13,14,15,16,17]
        s4_strips = [11,12,13,14,15,16,17]

      if i == 5 or i == 6 or i == 7 or i == 8:
        s1_strips = [-strip + 30 for strip in s1_strips]
        s2_strips = [-strip + 30 for strip in s2_strips]
        s3_strips = [-strip + 30 for strip in s3_strips]
        s4_strips = [-strip + 30 for strip in s4_strips]

      # debug
      #print i, "s1", s1_strips
      #print i, "s2", [x-8 for x in s2_strips]
      #print i, "s3", [x-8 for x in s3_strips]
      #print i, "s4", [x-8 for x in s4_strips]

      # Print for cfi
      if i == 0:
        cfi_i = 0
        cfi_straight = 4 - (cfi_i + 1) / 2
        cfi_string = "%i,%i:%i,%i:%i,%i:%i,%i:%i" % (cfi_straight, s1_strips[-1], s1_strips[0], s2_strips[-1]-8, s2_strips[0]-8, s3_strips[-1]-8, s3_strips[0]-8, s4_strips[-1]-8, s4_strips[0]-8)
      elif i <= 4:
        cfi_i = i * 2 - 1
        cfi_straight = 4 - (cfi_i + 1) / 2
        cfi_string = "%i,%i:%i,%i:%i,%i:%i,%i:%i" % (cfi_straight, s1_strips[-1], s1_strips[0], s2_strips[-1]-8, s2_strips[0]-8, s3_strips[-1]-8, s3_strips[0]-8, s4_strips[-1]-8, s4_strips[0]-8)
      elif i <= 8:
        cfi_i = (i-4) * 2 - 0
        cfi_straight = 4 - (cfi_i + 1) / 2
        cfi_string = "%i,%i:%i,%i:%i,%i:%i,%i:%i" % (cfi_straight, s1_strips[0], s1_strips[-1], s2_strips[0]-8, s2_strips[-1]-8, s3_strips[0]-8, s3_strips[-1]-8, s4_strips[0]-8, s4_strips[-1]-8)
      cfi_strings[cfi_i] = cfi_string


      for strip in s1_strips:
        lay = 0
        h_strips_new.Fill(strip, (4-lay))
      for strip in s2_strips:
        lay = 1
        h_strips_new.Fill(strip, (4-lay))
      for strip in s3_strips:
        lay = 2
        h_strips_new.Fill(strip, (4-lay))
      for strip in s4_strips:
        lay = 3
        h_strips_new.Fill(strip, (4-lay))

      # Draw "same"
      c1.cd()
      frame.Draw()
      gStyle.SetPalette(100)  # kSolar
      h_strips_new.SetMinimum(1)
      h_strips_new.Draw("same col")
      for l in tpolylines:
        l.Draw()

      options.imgname = "new_patt%i" % i
      c1.Print("plots_pattern_merge/"+options.imgname+".png")
      continue  # end i loop

    # Print for cfi
    print "process.simEmtfDigis.spPRParams16.PatternDefinitions = ["
    for cfi_string in cfi_strings:
      print "    \"" + cfi_string + "\","
    print "]"
    print "process.simEmtfDigis.spPRParams16.UseSymmetricalPatterns = False"
