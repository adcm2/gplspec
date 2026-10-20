#include <gplspec/All>
#include <iostream>
#include <complex>
#include <array>
#include <fstream>
#include <cmath>

struct Map {
  double h;
  auto RadialMapping(int) const {
    double a=h;
    return [a](double r,double theta,double){return a*std::cos(theta)*r;};
  }
};
int main(){
  using C=std::complex<double>;
  constexpr int p=5,lmax=30;
  constexpr double maxstep=.1, outer=1.5, h=.001, tol=1e-12;
  const double L=EarthModels::EarthConstants<double>().LengthNorm();
  const double T=EarthModels::EarthConstants<double>().TimeNorm();
  const double M=EarthModels::EarthConstants<double>().MassNorm();
  auto sphere=GeneralEarthModels::SimpleModels::spherical_model::HomogeneousSphere(1.,1.,L,T,M);
  Map map{h}, zero{0.};
  GeneralEarthModels::Density3D model(sphere,PlanetaryModel::TomographyZeroModel(),zero,p,lmax,maxstep,outer);
  MatrixReplacement3D<C> A(model);
  Eigen::VectorXcd rhs=Gravity_Tools::FindForce(model);
  Eigen::VectorXcd b=Gravity_Tools::FindBoundaryPerturbationForce(model,map);
#ifdef BASELINE_PROBE
  Eigen::SparseMatrix<C> P=-model.SpectralElementInformation().fullmatrix<C>(0,lmax);
#else
  Eigen::SparseMatrix<C> P=Gravity_Tools::detail::SphericalPreconditionerMatrix(model,lmax);
#endif
  P.makeCompressed();
  using Solver=Eigen::BiCGSTAB<MatrixReplacement3D<C>,Eigen::SphericalGeometryPreconditioner<C>>;
  Solver solver; solver.compute(A); solver.preconditioner().addmatrix(P); solver.setTolerance(tol);
  Eigen::VectorXcd guess=solver.preconditioner().solve(rhs);
  Eigen::VectorXcd x=solver.solveWithGuess(rhs,guess);
  auto base_status=solver.info();
  auto base_info=std::array<double,4>{double(solver.iterations()),solver.error(),(rhs-A*x).norm(),rhs.norm()};
  Eigen::VectorXcd guessb=solver.preconditioner().solve(b);
  Eigen::VectorXcd xb=solver.solveWithGuess(b,guessb);
  auto bound_status=solver.info();
  Eigen::VectorXcd boundary_residual=b-A*xb;
  auto bound_info=std::array<double,4>{double(solver.iterations()),solver.error(),boundary_residual.norm(),b.norm()};
  std::cout<<"base status="<<base_status<<" it="<<base_info[0]<<" error="<<base_info[1]<<" absres="<<base_info[2]<<" rhs="<<base_info[3]<<"\n";
  std::cout<<"boundary status="<<bound_status<<" it="<<bound_info[0]<<" error="<<bound_info[1]<<" absres="<<bound_info[2]<<" rhs="<<bound_info[3]<<" guessnorm="<<guessb.norm()<<" solnorm="<<xb.norm()<<"\n";
#ifdef BASELINE_PROBE
  std::ofstream dump("/tmp/gplspec-phase3-resume/classical-components-baseline.csv");
#else
  std::ofstream dump("/tmp/gplspec-phase3-resume/classical-components-parent.csv");
#endif
  dump<<"index,rhs_re,rhs_im,b_re,b_im,guess_re,guess_im,x_re,x_im,res_re,res_im\n";
  for(Eigen::Index i=0;i<b.size();++i)dump<<i<<","<<rhs[i].real()<<","<<rhs[i].imag()<<","<<b[i].real()<<","<<b[i].imag()<<","<<guessb[i].real()<<","<<guessb[i].imag()<<","<<xb[i].real()<<","<<xb[i].imag()<<","<<boundary_residual[i].real()<<","<<boundary_residual[i].imag()<<"\n";
  const std::size_t stride=(lmax+1)*(lmax+1);
  for(int l=0;l<=lmax;++l){double bn=0,xn=0,rn=0; for(std::size_t r=0;r<rhs.size()/stride;++r)for(int m=-l;m<=l;++m){auto k=r*stride+l*l+l+m;bn+=std::norm(b[k]);xn+=std::norm(xb[k]);auto z=boundary_residual[k];rn+=std::norm(z);} if(bn>1e-12||xn>1e-12)std::cout<<"l="<<l<<" b="<<std::sqrt(bn)<<" x="<<std::sqrt(xn)<<" residual="<<std::sqrt(rn)<<"\n";}
}
