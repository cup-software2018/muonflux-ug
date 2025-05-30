#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"

using namespace std;

extern "C" {
void initialize_music_(const char * data1, const char * data2,
                       const char * data3, int, int, int);
void muon_transport_(double & x0, double & y0, double & z0, double & cx0,
                     double & cy0, double & cz0, double & emuin0,
                     double & depth0, double & tmu0);
}

double AngleCorrection(double cost);
double ModifiedGaisser(double E_mu, double theta_rad);

int main(int argc, char ** argv)
{
  // initialize music
  const char * data1 = "/home/cupsoft/Works/Muon/muon_music/music-eloss-sr.dat";
  const char * data2 = "/home/cupsoft/Works/Muon/muon_music/music-double-diff-rock.dat";
  const char * data3 = "/home/cupsoft/Works/Muon/muon_music/music-cross-sections-sr.dat";
  initialize_music_(data1, data2, data3, 100, 100, 100);

  gRandom->SetSeed(0);

  const double vdepth = atof(argv[1]);

  const double pi = 3.14159265358979323846;
  const double km2cm = 100. * 1000.;

  const double emu_min = 0.106;
  const double emu_max = 1.0e+06;

  double E0, cost, phi, X;
  double E_final, weight;

  TFile f(argv[2], "recreate");
  TTree * tree = new TTree("flux", "flux");
  tree->Branch("E0", &E0);
  tree->Branch("Ef", &E_final);
  tree->Branch("cost", &cost);
  tree->Branch("phi", &phi);
  tree->Branch("depth", &X);
  tree->Branch("w", &weight);

  const int N = 10000; // number of try

  double L = log10(emu_min);
  double U = log10(emu_max);

  double dL = U - L;
  double C = log(10.0) * dL;

  for (int i = 0; i < N; ++i) {
    double loge = gRandom->Uniform(L, U);
    E0 = pow(10, loge);
    cost = gRandom->Rndm();
    phi = 2 * pi * gRandom->Rndm();

    double dN0_dE0dOm = ModifiedGaisser(E0, cost);
    double w0 = E0 * dN0_dE0dOm;

    X = vdepth / cost;

    double x = 0, y = 0, z = 0;
    double cx = 0, cy = 0, cz = 0;
    double tmu = 0;

    double me = E0;
    double d = X * km2cm / 2.65;

    muon_transport_(x, y, z, cx, cy, cz, me, d, tmu);

    E_final = me;
    weight = C * w0;

    tree->Fill();
  }

  tree->Write();
  f.Close();

  return 0;
}

double AngleCorrection(double cost)
{
  double x = cost;
  double p[5] = {0.102573, -0.068287, 0.958633, 0.0407253, 0.817285};
  double deno = 1.0 + p[0] * p[0] + p[1] + p[3];
  double nume = x * x + p[0] * p[0] + p[1] * pow(x, p[2]) + p[3] * pow(x, p[4]);
  return sqrt(nume / deno);
}

double ModifiedGaisser(double E_mu, double theta_rad)
{
  double cost = cos(theta_rad);
  double cost_star = AngleCorrection(cost);

  const double A = 0.14; // normalization (cm⁻² s⁻¹ sr⁻¹ GeV⁻¹)
  const double gamma = 2.7;
  const double eps_pi = 115.0;    // pion critical energy (GeV)
  const double eps_K = 850.0;     // kaon critical energy (GeV)
  const double r_K_pi = 0.054;    // kaon-to-pion ratio
  const double E_mu_decay = 3.64; // muon decay scale (GeV)

  // Modified spectrum includes low-energy decay correction
  double E_corr = E_mu * (1.0 + E_mu_decay / (E_mu * cost_star));
  double spectrum = pow(E_corr, -gamma);

  double term_pi = 1.0 / (1.0 + 1.1 * E_mu * cost_star / eps_pi);
  double term_K = r_K_pi / (1.0 + 1.1 * E_mu * cost_star / eps_K);

  return A * spectrum * (term_pi + term_K); // (cm²·s·sr·GeV)⁻¹
}