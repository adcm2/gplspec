#include <gplspec/All>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <memory>
#include <algorithm>
using Clock=std::chrono::steady_clock;
struct Map { auto RadialMapping(int) const { return [](double r,double,double){return .2*r;}; } };
int main(int argc,char**argv){
 if(argc!=2) return 2;
 bool mapped=std::string(argv[1])=="mapped";
 using Model=GeneralEarthModels::Density3D;
 std::ostringstream diagnostic;
 auto* out=std::cout.rdbuf(diagnostic.rdbuf());
 std::unique_ptr<Model> model;
 if(mapped){auto sphere=GeneralEarthModels::SimpleModels::spherical_model::HomogeneousSphere(1.,1.,6371000.,3600.,5.972e24); model=std::make_unique<Model>(sphere,PlanetaryModel::TomographyZeroModel(),Map{},5,2,.1,1.4);}
 else model=std::make_unique<Model>(Model::SphericalThreeDimensionalPlanetFromFile("modeldata/prem.200","modeldata/S40RTS_dvs.nc",12,20,.05,1.2));
 double tolerance=mapped?1e-12:1e-14;
 Eigen::VectorXcd final;
 std::vector<double> times;
 double total=0;
 for(int i=-1;i<100&&(i<3||total<.5);++i){
 diagnostic.str(""); diagnostic.clear();
 auto start=Clock::now();
#ifdef ORIGINAL_API
 auto result=Gravity_Tools::FindGravitationalPotential(*model,tolerance);
#else
 auto result=Gravity_Tools::SolveGravitationalPotential(*model,tolerance);
#endif
 double duration=std::chrono::duration<double>(Clock::now()-start).count();
#ifdef ORIGINAL_API
 const auto& mesh=model->Node_Information(); const int p=model->Poly_Order(); const int stride=(model->GSH_Grid().MaxDegree()+1)*(model->GSH_Grid().MaxDegree()+1);
 final.resize((model->Num_Elements()*p+1)*stride);
 for(std::size_t e=0;e<model->Num_Elements();++e) for(int j=0;j<=p;++j) for(int lm=0;lm<stride;++lm) final[(e*p+j)*stride+lm]=result[e][j][lm];
#else
 final=result.solution;
#endif
 if(i>=0){times.push_back(duration);total+=duration;}
 }
 MatrixReplacement3D<std::complex<double>> op(*model);
 Eigen::VectorXcd rhs=Gravity_Tools::FindForce(*model);
 Eigen::VectorXcd residual=rhs-op*final;
 std::cout.rdbuf(out);
 std::cout<<std::setprecision(17)<<"scenario,trial,seconds,coefficient_norm,absolute_residual,relative_residual\n";
 for(std::size_t i=0;i<times.size();++i)std::cout<<argv[1]<<','<<i<<','<<times[i]<<','<<final.norm()<<','<<residual.norm()<<','<<residual.norm()/rhs.norm()<<'\n';
 return final.allFinite()&&residual.allFinite()?0:1;
}
