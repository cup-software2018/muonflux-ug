#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <array>
#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <unistd.h> // for chdir()
#include <vector>

#include "TFile.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"

using namespace std;

extern "C" {
// Three CHARACTER*(*) args => three pointers + three size_t lengths
void initialize_music_(const char * file1, const char * file2,
                       const char * file3, size_t len1, size_t len2,
                       size_t len3);
// your muon_transport_ stays as-is (10 doubles by ref)
void muon_transport_(double & x0, double & y0, double & z0, double & cx0,
                     double & cy0, double & cz0, double & emuin0,
                     double & depth0, double & tmu0);
}

double AngleCorrection(double cost);
double ModifiedGaisser(double E_mu, double theta_rad);

// CGAL typedefs
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_2;
using Info = double; // our z
using Vb = CGAL::Triangulation_vertex_base_with_info_2<Info, K>;
using Fb = CGAL::Triangulation_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using Delaunay = CGAL::Delaunay_triangulation_2<K, Tds>; 

std::array<double, 3> barycentric_weights(const Point & P, const Point & A,
                                          const Point & B, const Point & C);
double interpolate_z(const Delaunay & T, double x, double y);

void print_usage(const char * prog)
{
  std::cout << "Usage: " << prog << " -i input.txt -o output.root -n samples\n"
            << "  -i, --input   mountain profile data file (x y z)\n"
            << "  -o, --output  output ROOT file\n"
            << "  -n, --number  number of samples\n"
            << "  -w, --workdir  working directory before running\n"
            << "  -h, --help    show this help\n";
}

int main(int argc, char ** argv)
{
  std::string input_file, output_file, work_dir(".");
  size_t n_samples = 10000;

  static struct option long_opts[] = {{"input", required_argument, 0, 'i'},
                                      {"output", required_argument, 0, 'o'},
                                      {"number", required_argument, 0, 'n'},
                                      {"workdir", required_argument, 0, 'w'},
                                      {"help", no_argument, 0, 'h'},
                                      {0, 0, 0, 0}};

  int opt;
  while ((opt = getopt_long(argc, argv, "i:o:n:w:h", long_opts, nullptr)) !=
         -1) {
    switch (opt) {
      case 'i': input_file = optarg; break;
      case 'o': output_file = optarg; break;
      case 'n': n_samples = std::stoul(optarg); break;
      case 'w': work_dir = optarg; break;
      case 'h': print_usage(argv[0]); return 0;
      default: print_usage(argv[0]); return 1;
    }
  }

  if (input_file.empty() || output_file.empty()) {
    print_usage(argv[0]);
    return 1;
  }

  // initialize music
  string data1 = work_dir + "/music-eloss-sr.dat";
  string data2 = work_dir + "/music-double-diff-rock.dat";
  string data3 = work_dir + "/music-cross-sections-sr.dat";
  initialize_music_(data1.c_str(), data2.c_str(), data3.c_str(), data1.size(),
                    data2.size(), data3.size());

  gRandom->SetSeed(0);

  const double km2cm = 100. * 1000.;
  const double emu_min = 0.106;
  const double emu_max = 1.0e+07;

  double E0, theta, phi, X;
  double E_final, weight;

  TFile f(output_file.c_str(), "recreate");
  TTree * tree = new TTree("flux", "flux");
  tree->Branch("E0", &E0);
  tree->Branch("Ef", &E_final);
  tree->Branch("theta", &theta);
  tree->Branch("phi", &phi);
  tree->Branch("depth", &X);
  tree->Branch("w", &weight);

  // Load (x,y,z) and track bounding box
  std::ifstream fin(input_file.c_str());
  if (!fin) {
    std::cerr << "Cannot open " << input_file << "\n";
    return 1;
  }
  std::cout << input_file << " opened" << std::endl;

  std::vector<std::pair<Point, Info>> pts;
  double x, y, z;
  double xmin = std::numeric_limits<double>::infinity();
  double ymin = std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();

  int nline = 0;
  while (fin >> x >> y >> z) {
    x /= 1000.;
    y /= 1000.;
    z /= 1000.;
    pts.emplace_back(Point(x, y), z);
    xmin = std::min(xmin, x);
    xmax = std::max(xmax, x);
    ymin = std::min(ymin, y);
    ymax = std::max(ymax, y);

    ++nline;
  }
  fin.close();

  std::cout << nline << " lines read" << std::endl;

  Delaunay tri;
  tri.insert(pts.begin(), pts.end());

  const int N = n_samples; // number of try

  double L = log10(emu_min);
  double U = log10(emu_max);

  double dL = U - L;
  double C = log(10.0) * dL;

  int generated = 0;
  while (generated < N) {
    // position on mountain
    double x0 = gRandom->Uniform(xmin, xmax);
    double y0 = gRandom->Uniform(ymin, ymax);
    double z0 = interpolate_z(tri, x0, y0);
    if (z0 == 0.0) continue;

    X = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);
    theta = std::acos(z0 / X);
    phi = std::atan2(y0, x0);

    // muon energy
    double loge = gRandom->Uniform(L, U);
    E0 = pow(10, loge);

    double dN0_dE0dOm = ModifiedGaisser(E0, theta);
    double w0 = E0 * dN0_dE0dOm;

    // muon transport by music
    double cx = 0, cy = 0, cz = 0;
    double tmu = 0;

    double me = E0;
    double d = X * km2cm;

    muon_transport_(x0, y0, z0, cx, cy, cz, me, d, tmu);

    E_final = me;
    weight = C * w0;

    tree->Fill();

    ++generated;
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

// Compute barycentric coords λ0,λ1,λ2 for P wrt triangle (A, B, C)
std::array<double, 3> barycentric_weights(const Point & P, const Point & A,
                                          const Point & B, const Point & C)
{
  double x = P.x(), y = P.y();
  double x1 = A.x(), y1 = A.y();
  double x2 = B.x(), y2 = B.y();
  double x3 = C.x(), y3 = C.y();

  double denom = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
  double λ0 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denom;
  double λ1 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denom;
  double λ2 = 1.0 - λ0 - λ1;
  return {λ0, λ1, λ2};
}

double interpolate_z(const Delaunay & T, double x, double y)
{
  Point q(x, y);
  auto fh = T.locate(q);
  if (T.is_infinite(fh)) return 0.0;

  // triangle vertices
  Point A = fh->vertex(0)->point();
  Point B = fh->vertex(1)->point();
  Point C = fh->vertex(2)->point();
  auto w = barycentric_weights(q, A, B, C);

  // weighted sum of z
  return w[0] * fh->vertex(0)->info() + w[1] * fh->vertex(1)->info() +
         w[2] * fh->vertex(2)->info();
}