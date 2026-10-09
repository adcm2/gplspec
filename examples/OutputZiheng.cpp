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
#include <numbers>
#include <numeric>

int
main() {
   using namespace GeneralEarthModels;
   using namespace Gravity_Tools;
   using namespace SimpleModels;
   using FLOAT = double;
   std::string pathtofile = "modeldata/TABLEA1.DAT";
   Timer timer1;
   timer1.start();
   const auto phobos_data =
       GPLSpec::ModelInput::ReadPhobosHarmonics(pathtofile);
   double physdensity = 1860.0;
   double lengthnorm = 10000.0;
   double usedensity = 1342.0;
   double massnorm = std::pow(lengthnorm, 3.0) * usedensity;
   double timenorm = 3600.0;
   double scaledensity = physdensity / massnorm * std::pow(lengthnorm, 3.0);
   int lMax = 100;
   int npoly = 5;
   Density3D phobos(scaledensity, phobos_data, npoly, lMax, lengthnorm, timenorm,
                    massnorm, 0.1, 1.5);

   {
      const std::size_t model_lmax = phobos_data.InferredMaximumDegree();
      // appending to vec_A if model lmax smaller than lmax:
      auto coefficientnumber =
          GSHTrans::GSHIndices<GSHTrans::NonNegative>(lMax, lMax, 0).Size();
      auto coefficientnumberall =
          GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).Size();
      std::vector<std::complex<double>> vec_lm_rad(coefficientnumber, 0.0);
      // std::vector<std::complex<double>> vec_lm_rad_full(coefficientnumberall,
      // 0.0);
      //   int maxmodell = lModMax;
      {
         std::complex<double> imag1(0.0, 1.0);
         std::size_t idxoverall = 0;
         for (int l = 0; l < model_lmax + 1; ++l) {
            for (int m = 0; m < l + 1; ++m) {
               double multfact = 0.5;
               // multfact = 1.0;
               if (m == 0) {
                  multfact = 1.0;
               }
               vec_lm_rad[idxoverall] =
                   multfact *
                   (phobos_data.rows[idxoverall].cosine -
                    imag1 * phobos_data.rows[idxoverall].sine);
               ++idxoverall;
            }
         }
      }

      // std::cout << "Got spatialsize\n";
      using GRID =
          GSHTrans::GaussLegendreGrid<double, GSHTrans::All, GSHTrans::All>;
      GRID _grid(lMax, 0);
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<double> vec_outerradius(spatialsize, 0.0);
      std::vector<std::complex<double>> vec_outer_comp(spatialsize, 0.0);
      std::vector<std::complex<double>> vec_lm_check(coefficientnumberall, 0.0);
      _grid.InverseTransformation(lMax, 0, vec_lm_rad, vec_outerradius);
      for (int idx = 0; idx < spatialsize; ++idx) {
         vec_outer_comp[idx] = vec_outerradius[idx];
      }
      _grid.ForwardTransformation(lMax, 0, vec_outer_comp, vec_lm_check);

      double theta = std::numbers::pi_v<double> / 4.0;
      auto wigtemp = GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All,
                                      GSHTrans::All, GSHTrans::Single,
                                      GSHTrans::ColumnMajor>(
          model_lmax, model_lmax, model_lmax, theta);
      auto dl = wigtemp[0];

      // angles
      //   auto vec_latitudes = std::vector<double>(719, 0.0);
      auto vec_latitudes = std::vector<double>(359, 0.0);
      auto vec_longitudes = std::vector<double>(720, 0.0);
      std::generate(vec_latitudes.begin(), vec_latitudes.end(),
                    [n = -90.0]() mutable {
                       n += 0.5;
                       return n;
                    });
      std::generate(vec_longitudes.begin(), vec_longitudes.end(),
                    [n = -0.5]() mutable {
                       n += 0.5;
                       return n;
                    });

      std::vector<std::vector<std::complex<double>>> vec_radius(
          vec_latitudes.size(),
          std::vector<std::complex<double>>(vec_longitudes.size(), 0.0));
      double pid180 = 3.141592653589793238462643383279 / 180.0;
      auto i1 = std::complex<double>(0.0, 1.0);
      for (int idx = 0; idx < vec_latitudes.size(); ++idx) {
         double colat = (90.0 - vec_latitudes[idx]) * pid180;
         auto wigtemp = GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All,
                                         GSHTrans::All, GSHTrans::Single,
                                         GSHTrans::ColumnMajor>(
             model_lmax, model_lmax, model_lmax, colat);
         auto dl = wigtemp[0];
         std::size_t idxoverall = 0;
         for (int idx2 = 0; idx2 < vec_longitudes.size(); ++idx2) {
            double lon = vec_longitudes[idx2] * pid180;
            for (int l = 0; l < model_lmax + 1; ++l) {
               for (int m = -l; m < l + 1; ++m) {
                  std::complex<double> ylm =
                      dl[l, m] * std::exp(i1 * double(m) * lon);
                  vec_radius[idx][idx2] += ylm * vec_lm_check[l * l + l + m];
               }
            }
         }
      }

      std::string pathtofile1 = "work/Ziheng/PhobosRadius.out";
      std::ofstream file1(pathtofile1);
      if (!file1) {
         std::cerr << "Error: unable to open output file: " << pathtofile1
                   << "\n";
         return 1;
      }

      file1.setf(std::ios::fixed);
      file1 << std::setprecision(16);

      double maxrad = 0.0, minrad = 1e10;
      for (int idx = 0; idx < vec_latitudes.size(); ++idx) {
         for (int idx2 = 0; idx2 < vec_longitudes.size(); ++idx2) {
            file1 << vec_latitudes[idx] << ";" << vec_longitudes[idx2] << ";"
                  << vec_radius[idx][idx2].real() << "\n";

            maxrad = std::max(maxrad, vec_radius[idx][idx2].real());
            minrad = std::min(minrad, vec_radius[idx][idx2].real());
            if (std::abs(vec_radius[idx][idx2].imag()) > 1e-10) {
               std::cout << "Warning: non-negligible imaginary part at "
                         << vec_latitudes[idx] << " " << vec_longitudes[idx2]
                         << " with value " << vec_radius[idx][idx2].imag()
                         << "\n";
            }
         }
      }
      file1.close();
      auto surface_metadata = GPLSpec::FieldBundle::RunMetadata{};
      surface_metadata.workflow = "output_ziheng";
      surface_metadata.dataset = "phobos_tablea1_surface_0p5deg";
      surface_metadata.inputs.push_back(
          GPLSpec::FieldBundle::InputPath(pathtofile));
      surface_metadata.parameters = {{"latitude_step_deg", "0.5"},
          {"longitude_step_deg", "0.5"},
          {"maximum_input_degree", std::to_string(model_lmax)}};
      GPLSpec::FieldBundle::Writer(
          "work/Ziheng/PhobosRadius.bundle", phobos,
          std::move(surface_metadata)).Write([&](auto &&emit) {
         for (std::size_t ilat = 0; ilat < vec_latitudes.size(); ++ilat) {
            for (std::size_t ilon = 0; ilon < vec_longitudes.size(); ++ilon) {
               GPLSpec::FieldEvaluation::Sample sample;
               sample.quantity = GPLSpec::FieldEvaluation::Quantity::SurfaceRadius;
               sample.representation = GPLSpec::FieldEvaluation::Representation::AngularGrid;
               sample.theta_rad = (90.0 - vec_latitudes[ilat]) * pid180;
               sample.phi_rad = vec_longitudes[ilon] * pid180;
               sample.physical_radius_m = vec_radius[ilat][ilon].real();
               sample.value = {sample.physical_radius_m, 0.0};
               emit(sample);
            }
         }
      });
      double minval =
          *std::min_element(vec_outerradius.begin(), vec_outerradius.end());

      std::cout << "\nMax radius grid: " << maxrad
                << ", min radius grid: " << minrad << " " << minval << "\n\n";
   }
   timer1.stop("Time to construct");

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////
   timer1.start();
   auto phobos_result = SolvePotential(phobos, std::pow(10.0, -12.0));
   auto stdvec_potsol = ToHarmonicField(phobos, phobos_result);
   timer1.stop("Time for gravity");
   std::cout << "Size of stdvec_potsol: " << stdvec_potsol.size() << " "
             << stdvec_potsol[0].size() << " " << stdvec_potsol[0][0].size()
             << "\n";

   std::cout << "Volume: " << phobos.Volume() / std::pow(10.0, 9.0) << "km^3"
             << "\n";
   std::cout << "Density: " << phobos.Mass() / (phobos.Volume() * 1000.0)
             << " g/cm^3\n";

   std::cout << "Outermost radius: "
             << phobos.Node_Information().OuterRadius() * phobos.LengthNorm()
             << " m\n";
   // output outermost:
   {
      std::string pathtofile1 = "work/Ziheng/PHOBOS_POTENTIAL_OUTER.out";
      std::ofstream file1(pathtofile1);
      if (!file1) {
         std::cerr << "Error: unable to open output file: " << pathtofile1
                   << "\n";
         return 1;
      }
      file1.setf(std::ios::fixed);
      file1 << std::setprecision(16);
      std::size_t fidx = 0;
      for (int idxl = 0; idxl < lMax + 1; ++idxl) {
         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
            file1 << idxl << ";" << idxm << ";"
                  << stdvec_potsol.back().back()[fidx].real() << ";"
                  << stdvec_potsol.back().back()[fidx++].imag() << "\n";
         }
      }
      file1.close();
      auto outer_metadata = GPLSpec::FieldBundle::WithSolve(
          phobos_result, "output_ziheng", "phobos_outer_potential_coefficients",
          "CG", phobos_result.final_solve().requested_tolerance);
      outer_metadata.inputs.push_back(
          GPLSpec::FieldBundle::InputPath(pathtofile));
      outer_metadata.parameters = {{"degree_max", std::to_string(lMax)},
          {"polynomial_order", std::to_string(npoly)},
          {"length_norm_m", std::to_string(lengthnorm)},
          {"mass_norm_kg", std::to_string(massnorm)},
          {"time_norm_s", std::to_string(timenorm)}};
      GPLSpec::FieldBundle::Writer(
          "work/Ziheng/PHOBOS_POTENTIAL_OUTER.bundle", phobos,
          std::move(outer_metadata)).Write([&](auto &&emit) {
         const auto last_element = static_cast<std::size_t>(
             phobos.Node_Information().NumberOfElements() - 1);
         GPLSpec::FieldEvaluation::ForEachRadialHarmonicSample(
             phobos, phobos_result.solution,
             GPLSpec::FieldEvaluation::Quantity::Potential,
             [&](const GPLSpec::FieldEvaluation::Sample &sample) {
                if (sample.element == last_element &&
                    sample.radial_node == phobos.Node_Information().PolynomialOrder())
                   emit(sample);
             });
      });
   }
   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////
   {   // Path to the PREM model data file.
      std::string pathtoprem = "modeldata/prem.200.no";

      Density3D testprem =
          Density3D::OneDimensionalPlanetFromFile(pathtoprem, 5, 2, 0.001, 1.2);
      auto vec_integral_potential = GravitationalSphericalIntegral(testprem);

      std::string pathtofile1 = "work/Ziheng/PREM200.NO.INTEGRAL.out";
      std::ofstream file1(pathtofile1);
      if (!file1) {
         std::cerr << "Error: unable to open output file: " << pathtofile1
                   << "\n";
         return 1;
      }

      file1.setf(std::ios::fixed);
      file1 << std::setprecision(16);
      file1 << 0.0 << ";" << vec_integral_potential[0][0].real() << "\n";
      for (int idx = 0; idx < vec_integral_potential.size() - 1; ++idx) {
         file1 << testprem.Node_Information().ElementUpperRadius(idx) << ";"
               << vec_integral_potential[idx + 1][0].real() << "\n";
      }
      file1.close();
      auto integral_metadata = GPLSpec::FieldBundle::RunMetadata{};
      integral_metadata.workflow = "output_ziheng";
      integral_metadata.dataset = "prem200_no_spherical_integral";
      integral_metadata.inputs.push_back(
          GPLSpec::FieldBundle::InputPath(pathtoprem));
      integral_metadata.parameters = {{"polynomial_order", "5"},
          {"degree_max", "2"}, {"max_step_norm", "0.001"},
          {"exterior_radius_norm", "1.2"}};
      GPLSpec::FieldBundle::Writer(
          "work/Ziheng/PREM200.NO.INTEGRAL.bundle", testprem,
          std::move(integral_metadata)).Write([&](auto &&emit) {
         GPLSpec::FieldEvaluation::ForEachElementBoundaryHarmonics(
             testprem, vec_integral_potential,
             GPLSpec::FieldEvaluation::Quantity::Potential, emit);
      });
   }
   // get gravitational field

   //////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////

   /*
      // output test
      std::string pathtofolder1 = "./work/Phobos/Cartesian";
      std::string pathtofolder2 = "./work/Phobos/Spherical";
      std::string pathtofolder3 = "./work/Phobos/Sensitivity";
      std::string pathtofile1 =
          pathtofolder2 + "/MatrixSolutionReferentialRotated.out";
      std::string pathtofile2 = pathtofolder2 + "/MatrixSolutionPhysical.out";
      std::string pathtofile3 = pathtofolder2 + "/MatrixSolutionRotated.out";
      std::string pathtofile4 =
          pathtofolder2 + "/PhysicalDensitySolutionRotated.out";
      std::string pathtofile5 =
          pathtofolder2 + "/ReferentialDensitySolutionRotated.out";

      timer1.start();
      phobos.CartesianOutputAtElement(pathtofolder1, stdvec_potsol);
      timer1.stop("First output");

      timer1.start();
      phobos.PhysicalOutputAtElement(pathtofolder2, stdvec_potsol);
      timer1.stop("Second output");

      timer1.start();
      phobos.ReferentialOutputAtElement(pathtofolder3, stdvec_senskernel, true);
      timer1.stop("Third output");

      timer1.start();
      phobos.ReferentialOutputSlice(pathtofolder3, stdvec_senskernel, true);
      timer1.stop("Fifth output");

      timer1.start();
      phobos.PhysicalOutputSlice(pathtofolder2, stdvec_potsol);
      timer1.stop("Sixth output");

      timer1.start();
      phobos.ReferentialOutputRotated(pathtofile1, vec_ang1, vec_ang2,
                                      stdvec_potsol);
      timer1.stop("Referential rotated");

      timer1.start();
      phobos.PhysicalOutputRotated(pathtofile2, vec_ang12, vec_ang22,
                                   stdvec_potsol);
      timer1.stop("Physical rotated");
   */
   return 0;
}
