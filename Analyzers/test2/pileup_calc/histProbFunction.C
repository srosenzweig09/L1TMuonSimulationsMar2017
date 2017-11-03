{
  TFile* _file0 = TFile::Open("MyDataPileupHistogram.root");
  Int_t nbins = 80;
  Int_t xlow = 33;
  Int_t xhi = 75;
  //
  TH1D* pileup = (TH1D*) _file0->Get("pileup");
  TH1D* pileup_test = (TH1D*) pileup->Clone("pileup_test");
  pileup_test->SetDirectory(0);
  // Suppress low pileup
  assert(pileup_test->GetNbinsX() == nbins);
  for (int i=0; i<xlow; ++i) {
    pileup_test->SetBinContent(i+1, 0);
    pileup_test->SetBinError(i+1, 0);
  }
  // Suppress high pileup
  for (int i=xhi+1; i<nbins; ++i) {
    pileup_test->SetBinContent(i+1, 0);
    pileup_test->SetBinError(i+1, 0);
  }
  // Normalize
  pileup_test->Scale(1.0/pileup_test->Integral());
  // Print 
  for (int i=0; i<xhi; ++i) {
    if (i==0)  std::cout << "0,";
    std::cout << i+1 << ",";
  }
  std::cout << std::endl;
  for (int i=0; i<xhi; ++i) {
    if (i==0)  std::cout << "0,";
    std::cout << pileup_test->GetBinContent(i+1) << ",";
  }
  std::cout << std::endl;
  //
  _file0->Close();
}
