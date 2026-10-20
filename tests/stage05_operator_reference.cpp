#include <gplspec/All>

#include <complex>
#include <iomanip>
#include <iostream>

int main() {
  using Scalar = std::complex<double>;
  using Operator = MatrixReplacement<Scalar>;
  using Grid = Operator::Grid;

  constexpr int lmax = 2;
  Grid grid(lmax, 2);
  auto quadrature =
      GaussQuad::GaussLobattoLegendreQuadrature1D<double>(4);
  const double length = 6371000.0;
  const double time = 3600.0;
  const double mass = 5.972e24;
  const std::string model_path =
      std::string(GPLSPEC_SOURCE_DIR) + "/modeldata/PREM500";
  GeneralEarthModels::spherical_1D model(model_path, quadrature, 0.5, 1.2);
  Operator op(model, grid);

  Eigen::VectorXcd x(op.cols()), y(op.cols());
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    x(i) = Scalar(0.125 * (i + 1), -0.0625 * (i + 2));
    y(i) = Scalar(-0.03125 * (i + 3), 0.09375 * (i + 1));
  }

  const Eigen::VectorXcd ax = op * x;
  const auto op_without_boundary = op.WithBoundary(false);
  const Eigen::VectorXcd ax_without_boundary = op_without_boundary * x;
  const Eigen::VectorXcd ay = op * y;
  const Scalar xay = x.dot(ay);
  const Scalar ax_y = ax.dot(y);
  const auto force = Gravity_Tools::FindForce(model, lmax);
  const auto &mesh = model.Node_Information();
  std::vector<std::size_t> interface_indices{0};
  for (std::size_t element = 1; element < mesh.NumberOfElements(); ++element) {
    if (mesh.LayerNumber(element) != mesh.LayerNumber(element - 1))
      interface_indices.push_back(element);
  }

  std::cout << std::setprecision(17);
  std::cout << "rows," << op.rows() << '\n';
  for (Eigen::Index i = 0; i < ax.size(); ++i) {
    std::cout << "Ax," << i << ',' << ax(i).real() << ',' << ax(i).imag() << '\n';
    const Scalar boundary_delta = ax(i) - ax_without_boundary(i);
    std::cout << "boundary_delta," << i << ',' << boundary_delta.real() << ','
              << boundary_delta.imag() << '\n';
  }
  for (Eigen::Index i = 0; i < force.size(); ++i)
    std::cout << "force," << i << ',' << force(i).real() << ','
              << force(i).imag() << '\n';
  std::cout << "xAy," << xay.real() << ',' << xay.imag() << '\n';
  std::cout << "Ax_y," << ax_y.real() << ',' << ax_y.imag() << '\n';
  std::cout << "origin,0," << mesh.ElementLowerRadius(0) << ','
            << ax(0).real() << ',' << ax(0).imag() << '\n';
  for (std::size_t layer = 1; layer < interface_indices.size(); ++layer) {
    const auto radial_node = interface_indices[layer];
    const auto spectral_index = radial_node * (model.Poly_Order()) *
                                GSHTrans::GSHIndices<GSHTrans::All>(
                                    lmax, lmax, 0).Size();
    std::cout << "interface," << layer << ','
              << mesh.NodeRadius(radial_node, 0)
              << ',' << ax(spectral_index).real() << ','
              << ax(spectral_index).imag() << '\n';
  }
  return 0;
}
