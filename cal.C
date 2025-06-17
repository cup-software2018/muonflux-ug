void cal()
{
  double E0, theta, phi, depth;
  double E_final, w;

  TChain * tree = new TChain("flux");
  tree->Add("data/flat_d1_4.root");
  
  tree->SetBranchAddress("E0", &E0);
  tree->SetBranchAddress("Ef", &E_final);
  tree->SetBranchAddress("theta", &theta);
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
  TH1D * h_dJdphi = new TH1D("h_dJdphi", "", 90, 0, 360);

  TH2D * h_dJdthetadphi =
      new TH2D("h_dJdthetadphi", "", 90, 0, 90, 360, 0, 360);
  TProfile2D * h_dEdthetadphi =
      new TProfile2D("dEdthetadphi", "", 90, 0, 90, 360, 0, 360);

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

    double t = theta * 180 / TMath::Pi();
    double p = phi * 180 / TMath::Pi();

    h_dJdE->Fill(E_final, weight);
    h_dJdtheta->Fill(t, weight);
    h_dJdphi->Fill(p, weight);
    h_dJdthetadphi->Fill(t, p, weight);
    h_dEdthetadphi->Fill(t, p, E_final, weight);
  }

  double Jm = 2 * TMath::Pi() * S_I;
  double Em = 2 * TMath::Pi() * S_E / Jm;

  cout << Form("%e %.2f", Jm, Em) << endl;

  h_dJdE->Scale(1.0, "width");
  h_dJdtheta->Scale(1.0, "width");
  h_dJdphi->Scale(1.0, "width");

  TCanvas * can1 = new TCanvas("can1", "", 1600, 1200);
  gPad->SetLogx();
  gPad->SetLogy();
  h_dJdE->Draw("L");

  TCanvas * can2 = new TCanvas("can2", "", 1600, 1200);
  h_dJdtheta->Draw("L");

  TCanvas * can3 = new TCanvas("can3", "", 1600, 1200);
  h_dJdphi->Draw("L");

  TCanvas * can4 = new TCanvas("can4", "", 1600, 1600);
  h_dJdthetadphi->Draw("colz");

  TCanvas * can5 = new TCanvas("can5", "", 1600, 1600);
  h_dEdthetadphi->Draw("colz");
}