#include <gplspec/All>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>

using Complex = std::complex<double>;
using Model = GeneralEarthModels::Density3D;
using Op = MatrixReplacement3D<Complex>;
using Pre = Eigen::SphericalGeometryPreconditioner<Complex>;
using CG = Eigen::ConjugateGradient<Op, Eigen::Lower | Eigen::Upper, Pre>;

void line(const std::string& name, double v) { std::cout << name << ',' << std::setprecision(17) << v << '\n'; }

int main(int argc, char** argv) {
  struct HomogeneousMap {
    auto RadialMapping(int) const {
      return [](double r,double,double){return .2*r;};
    }
  };
  auto base = GeneralEarthModels::SimpleModels::spherical_model::HomogeneousSphere(
      1.,1.,6371000.,3600.,5.972e24);
  Model model(base,PlanetaryModel::TomographyZeroModel(),HomogeneousMap{},5,2,.1,1.4);
  Op a(model);
  auto b = Gravity_Tools::FindForce(model);
  Eigen::SparseMatrix<Complex> pmat = -model.SpectralElementInformation().fullmatrix<Complex>(0,2);
  pmat.makeCompressed();
  CG solver;
  solver.compute(a);
  solver.preconditioner().addmatrix(pmat);
  solver.setTolerance(1e-12);
  const Eigen::VectorXcd guess = solver.preconditioner().solve(b);
  line("rows",a.rows());
  line("b_norm",b.norm());
  line("guess_norm",guess.norm());
  const Eigen::VectorXcd ag1 = a*guess;
  const Eigen::VectorXcd ag2 = a*guess;
  const Eigen::VectorXcd direct1 = b - a*guess;
  const Eigen::VectorXcd materialized1 = b-ag1;
  const Eigen::VectorXcd materialized2 = b-ag2;
  line("ag1_norm",ag1.norm());
  line("matvec_repeat_difference",(ag1-ag2).norm()/b.norm());
  line("guess_direct_residual",direct1.norm()/b.norm());
  line("guess_materialized_residual1",materialized1.norm()/b.norm());
  line("guess_materialized_residual2",materialized2.norm()/b.norm());
  line("residual_eval_difference",(direct1-materialized1).norm()/b.norm());
  const Eigen::VectorXcd solution = solver.solveWithGuess(b,guess);
  line("iterations",solver.iterations());
  line("solver_info",solver.info());
  line("solver_error",solver.error());
  line("solution_guess_difference",(solution-guess).norm()/guess.norm());
  const Eigen::VectorXcd as = a*solution;
  line("solution_materialized_residual",(b-as).norm()/b.norm());
  line("solution_direct_residual",Eigen::VectorXcd(b-a*solution).norm()/b.norm());
  line("solution_action_difference",(as-ag1).norm()/b.norm());
  const Eigen::VectorXcd correction = solver.preconditioner().solve(b-as);
  line("preconditioned_true_residual_correction_fraction",correction.norm()/solution.norm());
  line("correction_action_relative_mismatch",(a*correction-(b-as)).norm()/(b-as).norm());
  const Eigen::VectorXcd refined = solution+correction;
  line("one_refinement_true_residual",Eigen::VectorXcd(b-a*refined).norm()/b.norm());
  Eigen::VectorXcd v(a.cols());
  for (Eigen::Index i=0; i<v.size(); ++i) v(i)=Complex(std::sin(.01*(i+1)),std::cos(.02*(i+1)));
  const double a_lower_bound = (a*v).norm()/v.norm();
  line("a_2norm_lower_bound",a_lower_bound);
  line("backward_error_upper_bound_using_a_lower_bound",(b-as).norm()/(b.norm()+a_lower_bound*solution.norm()));
}
