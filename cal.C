void cal()
{
  double E0, cost, phi, depth;
  double E_final, w;

  TChain * tree = new TChain("flux");
  tree->Add("data/flat_*.root");

  tree->SetBranchAddress("E0", &E0);
  tree->SetBranchAddress("Ef", &E_final);
  tree->SetBranchAddress("cost", &cost);
  tree->SetBranchAddress("phi", &phi);
  tree->SetBranchAddress("depth", &depth);
  tree->SetBranchAddress("w", &w);

  const int nBins = 1000;
  const double Emin = 1.0; // GeV
  const double Emax = 1e6; // GeV

  // 2) build the array of log-spaced edges
  double binEdges[nBins + 1];
  double logMin = std::log10(Emin);
  double logMax = std::log10(Emax);
  for (int i = 0; i <= nBins; ++i) {
    double t = double(i) / nBins; // 0 → 1
    binEdges[i] = std::pow(10, logMin + t * (logMax - logMin));
  }

  // 3) create the histogram with variable bins
  TH1D * h_dJdE = new TH1D("h_dJdE", "", nBins, binEdges);
  TH1D * h_dJdtheta = new TH1D("h_dJdtheta", "", 90, 0, 90);

  int N = tree->GetEntries();
  cout << Form("total number of toy sample: %d", N) << endl;

  double S_I = 0;
  double S_E = 0;

  for (int i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    if (E_final <= 0.0) continue;

    double weight = w / N;

    S_I += weight;
    S_E += weight * E_final;

    h_dJdE->Fill(E_final, weight);
    h_dJdtheta->Fill(acos(cost)*180/TMath::Pi(), weight);
  }

  double Jm = 2 * TMath::Pi() * S_I;
  double Em = 2 * TMath::Pi() * S_E / Jm;

  cout << Form("%e %.2f", Jm, Em) << endl;

  h_dJdE->Scale(1.0, "width");
  h_dJdtheta->Scale(1.0, "width");

  TCanvas * can1 = new TCanvas("can1", "", 1600, 1200);
  gPad->SetLogx();
  gPad->SetLogy();
  h_dJdE->Draw("L");

  TCanvas * can2 = new TCanvas("can2", "", 1600, 1200);
  h_dJdtheta->Draw("L");
}