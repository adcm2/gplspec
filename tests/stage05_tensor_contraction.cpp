#include "../gplspec/src/Shared_Utilities.h"

#include <array>
#include <complex>
#include <stdexcept>

namespace {
using Complex = std::complex<double>;
using Tensor = Eigen::Matrix3cd;
using Components = std::array<Complex, 3>;

void original_contraction(const Tensor &a, Complex gm, Complex g0, Complex gp,
                          Components &q) {
  q[0] -= a(0, 0) * gp;
  q[0] += a(0, 1) * g0;
  q[0] -= a(0, 2) * gm;
  q[1] -= a(1, 0) * gp;
  q[1] += a(1, 1) * g0;
  q[1] -= a(1, 2) * gm;
  q[2] -= a(2, 0) * gp;
  q[2] += a(2, 1) * g0;
  q[2] -= a(2, 2) * gm;
}

Complex sample(int seed) {
  return {0.125 * ((seed % 31) - 15), -0.0625 * ((seed % 23) - 11)};
}
}  // namespace

int main() {
  for (int sample_index = 0; sample_index < 1000; ++sample_index) {
    Tensor tensor;
    for (int row = 0; row < 3; ++row)
      for (int col = 0; col < 3; ++col)
        tensor(row, col) = sample(sample_index * 9 + row * 3 + col + 1);
    const Complex gm = sample(sample_index + 3001);
    const Complex g0 = sample(sample_index + 4001);
    const Complex gp = sample(sample_index + 5001);
    Components expected{sample(sample_index + 6001), sample(sample_index + 7001),
                        sample(sample_index + 8001)};
    Components actual = expected;

    original_contraction(tensor, gm, g0, gp, expected);
    GPLSpec::detail::ContractCanonicalTensorVector(
        tensor, gm, g0, gp, actual[0], actual[1], actual[2]);
    if (actual != expected)
      throw std::runtime_error("canonical contraction differs from original statements");
  }
}
