#include <GaussQuad/All>
#include <gplspec/All>
#include <gplspec/Test>
#include <gplspec/Timer>
#include <PlanetaryModel/All>
#include <TomographyModels/All>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <ranges>
#include <sstream>
#include <cmath>

int
main() {

   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;
   using namespace SimpleModels;
   using FLOAT = double;

   // declaring variables:
   double maxstep = 0.1;
   double ballrad = 1.5;
   int npoly = 5;
   //  int lmax = 1;
   double normrad = 1.0;
   double normrho = 1.0;
   int lmax;
   double h;
   double radius = 1.0;
   double density = 1.0;
   double lengthnorm = EarthModels::EarthConstants<FLOAT>().LengthNorm();
   double timenorm = EarthModels::EarthConstants<FLOAT>().TimeNorm();
   double massnorm = EarthModels::EarthConstants<FLOAT>().MassNorm();
   std::string pathtoprem = "modeldata/prem.200";
   std::string pathtotomo = "modeldata/S40RTS_dvs.nc";

   Timer timer1;

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////
   // declare model and find the potential
   // mapping class
   class inp_map {
    public:
      inp_map() {};
      inp_map(const double h, const double r) : _h{h}, _prad{r} {};
      auto RadialMapping(int i) const {
         auto lambdamap = [hmult = _h, pr = _prad](double r, double theta,
                                                   double phi) {
            double rscale = r / pr;
            // return hmult * std::sin(theta) * std::sin(phi) * r * (1.0 -
            // rscale);
            // return hmult * std::sin(theta) * std::cos(theta) * std::sin(phi)
            // *
            //        r * (1.0 - rscale);
            return hmult * std::cos(theta) * r * (1.0 - rscale);
            // return hmult * r;
         };
         return lambdamap;
      }

    private:
      double _h = 0.0;
      double _prad = 1.0;
   };

   h = 0.0;
   lmax = 20;

   std::cout << "Type in h: ";
   std::cin >> h;
   inp_map map0(h, 1.0);
   inp_map map1(2.0 * h, 1.0);
   // std::cout << "Type in lmax: ";
   // std::cin >> lmax2;

   // make model
   timer1.start();
   // Density3D testtomo = Density3D::SphericalThreeDimensionalPlanetFromFile(
   //     pathtoprem, pathtotomo, npoly, lmax, maxstep, ballrad);
   // Density3D testtomo = Density3D::SphericalHomogeneousPlanet(
   //     normrad, normrho, lengthnorm, timenorm, massnorm, maxstep, ballrad);
   SimpleModels::spherical_model _sphere1D =
       SimpleModels::spherical_model::HomogeneousSphere(
           normrad, normrho, lengthnorm, timenorm, massnorm);
   Density3D testtomo(_sphere1D, PlanetaryModel::TomographyZeroModel(), map0,
                      npoly, lmax, maxstep, ballrad);
   Density3D testtomo2(_sphere1D, PlanetaryModel::TomographyZeroModel(), map1,
                       npoly, lmax, maxstep, ballrad);
   timer1.stop("Time to make model");

   // do PSM solution
   timer1.start();
   auto potential_result = SolvePotential(testtomo, std::pow(10.0, -12.0));
   auto exact_result = SolvePotential(testtomo2, std::pow(10.0, -12.0));
   timer1.stop("Time for 3D");

   // solution using spherical integration method:
   timer1.start();
   auto vec_integral_potential = GravitationalSphericalIntegral(testtomo);
   timer1.stop("Time for integral");

   //    Eigen::VectorXcd testforce = FindBoundaryPerturbationForce(testtomo,
   //    map0); Eigen::VectorXcd testchange =
   //        AdvectiveBoundaryPerturbation(testtomo, map1, testforce);

   // testing perturbation class
   // std::cout << "HELLO\n";

   // find perturbed solution
   MappingPerturbation testpertmap(testtomo, map0);
   auto mapping_result = SolveMappingPerturbation(testtomo, testpertmap,
                                                  std::pow(10.0, -12.0));

   // variation with h:
   //    timer1.start();
   //    int lenh = 100;
   //    std::vector<double> vec_h(lenh);
   //    std::vector<double> vec_zerosol(lenh), vec_exactsol(lenh);
   //    for (int idx = 0; idx < lenh; ++idx) {
   //       vec_h[idx] = 0.2 / lenh * idx;
   //    }
   //    {
   //       double pi_db = 3.1415926535897932;

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   // output test
   std::string pathtofolder1 = "./work/Bench8/exact";
   std::string pathtofolder2 = "./work/Bench8/perturbed";
   std::string pathtofolder3 = "./work/Bench8/unperturbed";
   auto exact_bundle = GPLSpec::FieldBundle::WithSolve(
       exact_result, "clean_bench_8", "exact_mapped_sphere", "CG",
       exact_result.final_solve().requested_tolerance);
   exact_bundle.inputs.push_back(GPLSpec::FieldBundle::InputPath(pathtoprem));
   exact_bundle.inputs.push_back(GPLSpec::FieldBundle::InputPath(pathtotomo));
   exact_bundle.parameters = {{"polynomial_order", std::to_string(npoly)},
       {"degree_max", std::to_string(lmax)}, {"max_step_norm", GPLSpec::FieldBundle::ParameterNumber(maxstep)},
       {"exterior_radius_norm", GPLSpec::FieldBundle::ParameterNumber(ballrad)},
       {"mapping_amplitude", GPLSpec::FieldBundle::ParameterNumber(h)},
       {"length_norm_m", GPLSpec::FieldBundle::ParameterNumber(lengthnorm)},
       {"mass_norm_kg", GPLSpec::FieldBundle::ParameterNumber(massnorm)},
       {"time_norm_s", GPLSpec::FieldBundle::ParameterNumber(timenorm)}};
   GPLSpec::FieldBundle::WriteAngular(
       std::filesystem::path(pathtofolder1) / "field_bundle", testtomo2,
       exact_result.solution, GPLSpec::FieldEvaluation::Quantity::Potential,
       exact_bundle);
   auto map_bundle = GPLSpec::FieldBundle::WithSolve(
       mapping_result, "clean_bench_8", "mapping_perturbation", "CG",
       mapping_result.final_solve().requested_tolerance);
   map_bundle.inputs = exact_bundle.inputs;
   map_bundle.parameters = exact_bundle.parameters;
   GPLSpec::FieldBundle::WriteAngular(
       std::filesystem::path(pathtofolder2) / "field_bundle", testtomo,
       mapping_result.solution, GPLSpec::FieldEvaluation::Quantity::Potential,
       std::move(map_bundle));
   auto base_bundle = GPLSpec::FieldBundle::WithSolve(
       potential_result, "clean_bench_8", "unperturbed", "CG",
       potential_result.final_solve().requested_tolerance);
   base_bundle.inputs = exact_bundle.inputs;
   base_bundle.parameters = exact_bundle.parameters;
   GPLSpec::FieldBundle::WriteAngular(
       std::filesystem::path(pathtofolder3) / "field_bundle", testtomo,
       potential_result.solution, GPLSpec::FieldEvaluation::Quantity::Potential,
       std::move(base_bundle));

   //    std::string pathtofile = "./work/Bench8/zeroradius/MatrixSolution.out";
   //    auto file2 = std::ofstream(pathtofile);
   //    for (int idx = 0; idx < lenh; ++idx) {
   //       file2 << vec_h[idx] << ";" << vec_zerosol[idx] << ";" <<
   //       vec_exactsol[idx]
   //             << std::endl;
   //    }
   //    file2.close();
   //  vec_exactsol2);


   return 0;
}
