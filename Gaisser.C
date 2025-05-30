#include <cmath>

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

    const double A = 0.14;         // normalization (cm⁻² s⁻¹ sr⁻¹ GeV⁻¹)
    const double gamma = 2.7;
    const double eps_pi = 115.0;   // pion critical energy (GeV)
    const double eps_K = 850.0;    // kaon critical energy (GeV)
    const double r_K_pi = 0.054;   // kaon-to-pion ratio
    const double E_mu_decay = 3.64; // muon decay scale (GeV)

    // Modified spectrum includes low-energy decay correction
    double E_corr = E_mu * (1.0 + E_mu_decay / (E_mu * cost_star));
    double spectrum = pow(E_corr, -gamma);

    double term_pi = 1.0 / (1.0 + 1.1 * E_mu * cost_star / eps_pi);
    double term_K  = r_K_pi / (1.0 + 1.1 * E_mu * cost_star / eps_K);

    return A * spectrum * (term_pi + term_K); // (cm²·s·sr·GeV)⁻¹
}
