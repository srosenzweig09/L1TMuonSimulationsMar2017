{
  TFile* _file0 = TFile::Open("MyDataPileupHistogram.root", "UPDATE");
  //
  TH1D* pileup = (TH1D*) _file0->Get("pileup");
  TH1D* pileup_test = (TH1D*) pileup->Clone("pileup_test");
  // Suppress low pile-up
  assert(pileup_test->GetNbinsX() == 50);
  for (int i=0; i<15; ++i) {
    pileup_test->SetBinContent(i, 0);
    pileup_test->SetBinError(i, 0);
    if (i == 14) {
      pileup_test->SetBinContent(i, 5000);
      pileup_test->SetBinError(i, 5000);
    }
  }
  // Normalize
  pileup_test->Scale(1.0/pileup_test->Integral());
  // Print 
  for (int i=0; i<50; ++i) {
    std::cout << i << ",";
  }
  std::cout << std::endl;
  for (int i=0; i<50; ++i) {
    std::cout << pileup_test->GetBinContent(i+1) << ",";
  }
  std::cout << std::endl;
  //
  _file0->Close();
}
