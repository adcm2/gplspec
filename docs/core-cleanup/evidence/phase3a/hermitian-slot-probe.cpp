#include <gplspec/All>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>

using Complex = std::complex<double>;

template <class Operator>
void examine(const char* kind, Operator& op, bool boundary) {
  if (boundary) op.IncludeBoundary(); else op.NoBoundary();
  Eigen::VectorXcd x(op.cols()), y(op.cols());
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    x(i) = Complex(.125*(i+1), -.0625*(i+2));
    y(i) = Complex(-.03125*(i+3), .09375*(i+1));
  }
  int lmax = 2, stride = 9;
  for (int mode = 0; mode < 5; ++mode) {
    auto u = x.eval(), v = y.eval();
    for (Eigen::Index i = 0; i < u.size(); ++i) {
      int slot = i % stride;
      bool excluded = (mode >= 1 && mode <= 3 && slot == 8) ||
                      (mode >= 2 && mode <= 3 && slot == 4) ||
                      (mode >= 3 && i < stride && slot > 0);
      if (excluded) u(i) = v(i) = 0;
    }
    auto au = Eigen::VectorXcd(op*u), av = Eigen::VectorXcd(op*v);
    auto defect = u.dot(av) - au.dot(v);
    double scale = u.norm()*av.norm() + au.norm()*v.norm();
    double leak = 0;
    double total = au.squaredNorm();
    for (Eigen::Index i = 0; i < au.size(); ++i) {
      int slot = i % stride;
      bool excluded = (mode >= 1 && mode <= 3 && slot == 8) ||
                      (mode >= 2 && mode <= 3 && slot == 4) ||
                      (mode >= 3 && i < stride && slot > 0);
      if (excluded) leak += std::norm(au(i));
    }
    std::cout << kind << ',' << boundary << ',' << mode << ','
              << defect.real() << ',' << defect.imag() << ','
              << std::abs(defect)/scale << ',' << std::sqrt(leak/total) << '\n';
  }
}

struct ZeroMapping {
  auto RadialMapping(int) const { return [](double,double,double){return 0.0;}; }
};
struct SmoothMapping {
  auto RadialMapping(int) const {
    return [](double r,double theta,double phi){
      return .1*r*(1.-r)*std::sin(theta)*std::cos(phi);
    };
  }
};

void support(const char* name, const Eigen::VectorXcd& v) {
  double total = v.squaredNorm(), plus = 0., both = 0., origin = 0.;
  for (Eigen::Index i=0; i<v.size(); ++i) {
    double s = std::norm(v(i));
    if (i%9==8) plus += s;
    if (i%9==4 || i%9==8) both += s;
    if (i<9 && i%9>0) origin += s;
  }
  std::cout << "support," << name << ',' << v.size() << ','
            << std::sqrt(total) << ',' << std::sqrt(plus/total) << ','
            << std::sqrt(both/total) << ',' << std::sqrt(origin/total) << '\n';
}

int main() {
  std::cout << std::setprecision(17);
  using Op1 = MatrixReplacement<Complex>;
  Op1::Grid grid(2,2);
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(4);
  std::string source = GPLSPEC_SOURCE_DIR;
  GeneralEarthModels::spherical_1D one(source+"/modeldata/PREM500",q,.5,1.2);
  Op1 op1(one,grid);
  examine("one",op1,true);
  examine("one",op1,false);
  auto material = GeneralEarthModels::SimpleModels::spherical_model::HomogeneousSphere(1.0,1.0,6371000.,3600.,5.972e24);
  ZeroMapping zero;
  GeneralEarthModels::Density3D three(material,PlanetaryModel::TomographyZeroModel(),zero,3,2,.5,1.2);
  GeneralEarthModels::MappingPerturbation perturbation(three,source+"/tests/reference/stage05-file-perturbation-coefficients.dat",0,2);
  MatrixReplacement3D<Complex> op3(three,perturbation);
  examine("three",op3,true);
  examine("three",op3,false);
  MatrixReplacement3D<Complex> base3(three);
  examine("three_base",base3,true);
  SmoothMapping mapping;
  GeneralEarthModels::Density3D mapped(material,PlanetaryModel::TomographyZeroModel(),mapping,2,2,.5,1.2);
  MatrixReplacement3D<Complex> mappedop(mapped);
  examine("mapped_base",mappedop,true);
  auto force = Gravity_Tools::FindForce(mapped);
  support("mapped_force",force);
  Eigen::SparseMatrix<Complex> matrix = -mapped.SpectralElementInformation().fullmatrix<Complex>(0,mapped.GSH_Grid().MaxDegree());
  matrix.makeCompressed();
  Eigen::SphericalGeometryPreconditioner<Complex> pre;
  pre.compute(mappedop);
  pre.addmatrix(matrix);
  Eigen::VectorXcd guess = pre.solve(force);
  support("mapped_guess",guess);
  support("mapped_op_force",Eigen::VectorXcd(mappedop*force));
  support("mapped_op_guess",Eigen::VectorXcd(mappedop*guess));
  using Solver = Eigen::ConjugateGradient<MatrixReplacement3D<Complex>,
                                        Eigen::Lower | Eigen::Upper,
                                        Eigen::SphericalGeometryPreconditioner<Complex>>;
  Solver solver;
  solver.compute(mappedop);
  solver.preconditioner().addmatrix(matrix);
  solver.setTolerance(1e-6);
  Eigen::VectorXcd solution = solver.solveWithGuess(force,guess);
  support("mapped_solution",solution);
  support("mapped_op_solution",Eigen::VectorXcd(mappedop*solution));
  std::cout << "solve," << solver.iterations() << ',' << solver.error()
            << ',' << (force - mappedop*solution).norm()/force.norm() << '\n';
}
