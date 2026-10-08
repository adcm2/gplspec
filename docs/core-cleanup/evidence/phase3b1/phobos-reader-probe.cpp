#include <gplspec/GeneralModels>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
  std::vector<int> degree, order;
  std::vector<double> cosine, cosine_sigma, sine, sine_sigma;
#ifdef GPLSPEC_TYPED_INPUT
  const auto rows =
      GPLSpec::ModelInput::ReadPhobosHarmonics("modeldata/TABLEA1.DAT");
  for (const auto &row : rows.rows) {
    degree.push_back(row.degree);
    order.push_back(row.order);
    cosine.push_back(row.cosine);
    cosine_sigma.push_back(row.cosine_sigma);
    sine.push_back(row.sine);
    sine_sigma.push_back(row.sine_sigma);
  }
#else
  std::ifstream file("modeldata/TABLEA1.DAT");
  int l, m;
  double a, sa, b, sb;
  while (file >> l >> m >> a >> sa >> b >> sb) {
    degree.push_back(l);
    order.push_back(m);
    cosine.push_back(a);
    cosine_sigma.push_back(sa);
    sine.push_back(b);
    sine_sigma.push_back(sb);
  }
#endif
  const std::size_t count = degree.size();
  const std::size_t inferred_lmax =
      static_cast<std::size_t>((std::sqrt(1.0 + 8.0 * count) - 3.0) / 2.0);
  const std::size_t lmax = 100;
  const std::size_t lmodmax = 20;
  const std::size_t used_max = std::min({inferred_lmax, lmax, lmodmax});
  std::cout << std::setprecision(17);
  std::cout << "META," << count << ',' << inferred_lmax << ',' << used_max << '\n';
  for (std::size_t i = 0; i < count; ++i) {
    std::complex<double> coefficient{0.0, 0.0};
    bool retained = false;
    for (std::size_t l = 0, idx = 0; l <= used_max && !retained; ++l) {
      for (std::size_t m = 0; m <= l; ++m, ++idx) {
        if (idx == i) {
          const double factor = m == 0 ? 1.0 : 0.5;
          coefficient = factor * std::complex<double>(cosine[i], -sine[i]);
          retained = true;
          break;
        }
      }
    }
    std::cout << "R," << i << ',' << degree[i] << ',' << order[i] << ','
              << cosine[i] << ',' << cosine_sigma[i] << ',' << sine[i] << ','
              << sine_sigma[i] << ',' << retained << ',' << coefficient.real()
              << ',' << coefficient.imag() << '\n';
  }
}
