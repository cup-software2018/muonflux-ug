R__LOAD_LIBRARY(libGaisser)

void plot_gaisser()
{
  double dloge = 0.01;
  double logemin = -1;
  double logemax = 4;

  double degree[4] = {0, 60, 75, 80};

  TGraph * grp[4];
  for (size_t i = 0; i < 4; i++) {
    grp[i] = new TGraph();
  }

  int np = 0;
  for (double loge = logemin; loge <= logemax; loge += dloge) {
    double me = pow(10, loge);

    for (int i = 0; i < 4; i++) {
      double theta = degree[i] * TMath::Pi() / 180.;
      double gval = pow(me, 2.7) * ModifiedGaisser(me, theta);
      grp[i]->SetPoint(np, me, gval);
    }
    np += 1;
  }

  TCanvas * can1 = new TCanvas("can1", "", 1600, 1200);
  gPad->SetLogx();
  gPad->SetLogy();
  grp[0]->GetYaxis()->SetRangeUser(pow(10, -7), pow(10, 0));
  for (int i = 0; i < 4; i++) {
    const char * opt = (i == 0) ? "AL" : "L";
    grp[i]->Draw(opt);
  }
}