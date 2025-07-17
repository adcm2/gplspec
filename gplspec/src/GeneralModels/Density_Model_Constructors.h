#ifndef DENSITY_EARTH_MODEL_CONSTRUCTORS_3D_H
#define DENSITY_EARTH_MODEL_CONSTRUCTORS_3D_H

#include <random>
namespace GeneralEarthModels {
// template <class model>
//    requires PlanetaryModel::BasicSphericalDensityModel<model, int, double>
// Density3D::Density3D(const model &inp_model, const int npoly, const int lMax,
//                      double max_radial_step, double maxrad)
//     : Density3D::Density3D(inp_model, PlanetaryModel::TomographyZeroModel(),
//                            npoly, lMax, max_radial_step, maxrad){};
// trying factory function initialisation
Density3D
Density3D::SphericalHomogeneousPlanet(
    double radius, double density,
    double length_norm = EarthModels::EarthConstants<double>().LengthNorm(),
    double time_norm = EarthModels::EarthConstants<double>().TimeNorm(),
    double mass_norm = EarthModels::EarthConstants<double>().MassNorm(),
    double maxradialstep = 0.01, double ballrad = 1.2) {

   // SimpleModels::spherical_model _sphere1D(1.0, 1.0);
   SimpleModels::spherical_model _sphere1D =
       SimpleModels::spherical_model::HomogeneousSphere(
           radius, density, length_norm, time_norm, mass_norm);
   Density3D testmodel = Density3D::OneDimensionalPlanetFromModel(
       _sphere1D, 5, 2, maxradialstep, ballrad);
   return testmodel;
};

// trying factory function initialisation
Density3D
Density3D::SphericalLayeredPlanet(
    std::vector<double> radius, std::vector<double> density,
    double length_norm = EarthModels::EarthConstants<double>().LengthNorm(),
    double time_norm = EarthModels::EarthConstants<double>().TimeNorm(),
    double mass_norm = EarthModels::EarthConstants<double>().MassNorm(),
    double maxradialstep = 0.01, double ballrad = 1.2) {

   // SimpleModels::spherical_model _sphere1D(1.0, 1.0);
   SimpleModels::spherical_model _sphere1D =
       SimpleModels::spherical_model::HomogeneousLayers(
           radius, density, length_norm, time_norm, mass_norm);
   Density3D testmodel = Density3D::OneDimensionalPlanetFromModel(
       _sphere1D, 5, 2, maxradialstep, ballrad);
   return testmodel;
};

Density3D
Density3D::OneDimensionalPlanetFromFile(const std::string &pathto1D,
                                        const int npoly, const int lMax,
                                        double max_radial_step, double maxrad) {
   return Density3D(TestTools::EarthModel(pathto1D),
                    PlanetaryModel::TomographyZeroModel(), npoly, lMax,
                    max_radial_step, maxrad);
};

Density3D
Density3D::SphericalThreeDimensionalPlanetFromFile(
    const std::string &pathto1D, const std::string &pathtotomo, const int npoly,
    const int lMax, double max_radial_step, double maxrad) {
   class inp_map {
    public:
      inp_map() {};
      auto RadialMapping(int i) const {
         auto lambdamap = [](double r, double theta, double phi) {
            return 0.0;
         };
         return lambdamap;
      }
   };
   inp_map map1;
   return Density3D(TestTools::EarthModel(pathto1D), Tomography(pathtotomo),
                    npoly, lMax, max_radial_step, maxrad);
};

template <class model>
   requires PlanetaryModel::BasicSphericalDensityModel<model, int, double>
Density3D
Density3D::OneDimensionalPlanetFromModel(const model &inp_model,
                                         const int npoly, const int lMax,
                                         double max_radial_step,
                                         double maxrad) {
   return Density3D(inp_model, PlanetaryModel::TomographyZeroModel(), npoly,
                    lMax, max_radial_step, maxrad);
};

template <class model, class tomomodel>
   requires PlanetaryModel::BasicSphericalDensityModel<model, int, double> &&
            PlanetaryModel::TomographyModel<tomomodel>
Density3D
Density3D::SphericalThreeDimensionalPlanetFromModel(
    const model &inp_model, const tomomodel &tomo_model, const int npoly,
    const int lMax, double max_radial_step, double maxrad) {
   return Density3D(inp_model, tomo_model, npoly, lMax, max_radial_step,
                    maxrad);
};

template <class model, class tomomodel, class mapclass>
   requires PlanetaryModel::BasicSphericalDensityModel<model, int, double> &&
            PlanetaryModel::TomographyModel<tomomodel> &&
            PlanetaryModel::RadialMappingClass<mapclass>
Density3D
Density3D::PhysicalPlanetModelWithProvidedMapping(
    const model &inp_model, const tomomodel &tomo_model,
    const mapclass &inp_map, const int npoly, const int lMax,
    double max_radial_step, double maxrad) {
   return Density3D(inp_model, tomo_model, inp_map, npoly, lMax,
                    max_radial_step, maxrad);
};

// template <class model, class tomomodel>
//    requires PlanetaryModel::BasicSphericalDensityModel<model, int, double> &&
//             PlanetaryModel::TomographyModel<tomomodel>
// Density3D
// Density3D::AsphericalTest(const model &inp_model, const tomomodel
// &tomo_model,
//                           const int npoly, const int lMax,
//                           double max_radial_step, double maxrad) {
//    // testing model add
//    class simplemap {
//     public:
//       simplemap() {

//       };
//       auto Mapping(int i) {
//          auto lambdamap = [](double r, double theta, double phi) {
//             return 0.02 * r;
//          };
//          return lambdamap;
//       }
//    };

//    return Density3D(inp_model, tomo_model, npoly, lMax, max_radial_step,
//    maxrad,
//                     5);
// };

// aspherical constructor
//  full constructor
// template <class model, class tomomodel, class mapclass>
//    requires PlanetaryModel::BasicSphericalDensityModel<model, int, double> &&
//             PlanetaryModel::TomographyModel<tomomodel> &&
//             PlanetaryModel::MappingClass<mapclass>
// Density3D::Density3D(const model &inp_model, const tomomodel &tomo_model,
//                      const mapclass &inp_map, const int npoly, const int
//                      lMax, double max_radial_step, double maxrad, int
//                      randint)
//     : Density3D(inp_model, tomo_model, npoly, lMax, max_radial_step,
//     maxrad){};

// simple
//  constructor for homogeneous physical body with single layer. Path to file
//  contains the relative file path of a data file with the spherical harmonic
//  coefficients in the form of that for Phobos. maxradialstep is now relative
//  to the planet's referential radius, ie 0.01 means minimum 100 layers
Density3D::Density3D(double physicaldensity, std::string pathtofile,
                     const int npoly, const int lMax, double ilength_norm,
                     double itime_norm, double imass_norm, double maxradialstep,
                     double ballrad, bool randdensity)
    : Density3D(physicaldensity, pathtofile, lMax, npoly, lMax, ilength_norm,
                itime_norm, imass_norm, maxradialstep, ballrad, randdensity) {};

Density3D::Density3D(double physicaldensity, std::string pathtofile,
                     const int lModMax, const int npoly, const int lMax,
                     double ilength_norm, double itime_norm, double imass_norm,
                     double maxradialstep, double ballrad, bool randdensity)
    : length_norm(ilength_norm), mass_norm(imass_norm), time_norm(itime_norm),
      density_norm(imass_norm / std::pow(ilength_norm, 3.0)),
      velocity_norm(ilength_norm / itime_norm),
      acceleration_norm(ilength_norm / std::pow(itime_norm, 2.0)),
      force_norm(imass_norm * ilength_norm / std::pow(itime_norm, 2.0)),
      stress_norm(imass_norm / (std::pow(itime_norm, 2.0) * ilength_norm)),
      inertia_norm(imass_norm * std::pow(ilength_norm, 2.0)),
      gravitational_constant_norm(std::pow(ilength_norm, 3.0) /
                                  (imass_norm * std::pow(itime_norm, 2.0))),
      gravitational_constant(6.6743015 * std::pow(10.0, -11.0) * imass_norm *
                             std::pow(itime_norm, 2.0) /
                             std::pow(ilength_norm, 3.0)),
      potential_norm(std::pow(ilength_norm / itime_norm, 2.0)) {

   _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoly + 1);
   _grid = Grid(lMax, 2);
   // polynomial order
   _poly_ord = _q.N() - 1;

   // now this constructor deals with reading in Phobos data in particular
   // find length of file
   std::size_t num_lines = TestTools::filelength(pathtofile);
   std::size_t model_lmax = (std::sqrt(1 + 8 * num_lines) - 3) / 2;
   std::vector<int> vec_l(num_lines), vec_m(num_lines);
   std::vector<double> vec_A(num_lines), vec_SA(num_lines), vec_B(num_lines),
       vec_SB(num_lines);

   std::fstream modelfile;
   modelfile.open(pathtofile, std::ios::in);
   // getting information out of file
   if (modelfile.is_open()) {
      // test first line
      // getline(modelfile, modeltitle);
      for (int idx = 0; idx < num_lines; ++idx) {
         modelfile >> vec_l[idx] >> vec_m[idx] >> vec_A[idx] >> vec_SA[idx] >>
             vec_B[idx] >> vec_SB[idx];
         // move to next line
         modelfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }

      modelfile.close();
   } else {
      std::cout << "Couldn't open!" << "\n";
      assert("Model not found!");
   }
   // appending to vec_A if model lmax smaller than lmax:
   auto coefficientnumber =
       GSHTrans::GSHIndices<GSHTrans::NonNegative>(lMax, lMax, 0).size();
   auto coefficientnumberall =
       GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
   std::vector<std::complex<double>> vec_lm_rad(coefficientnumber, 0.0);
   // std::vector<std::complex<double>> vec_lm_rad_full(coefficientnumberall,
   // 0.0);
   int maxmodell = std::min(lModMax, lMax);
   {
      std::complex<double> imag1(0.0, 1.0);
      std::size_t idxoverall = 0;
      std::size_t maxl = 0, tmp_maxl = 0;
      if (lMax > model_lmax) {
         tmp_maxl = model_lmax;
      } else {
         tmp_maxl = lMax;
      }
      if (lModMax > tmp_maxl) {
         maxl = tmp_maxl;
      } else {
         maxl = lModMax;
      }
      // maxl = std::min(lModMax, tmp_maxl);
      for (int l = 0; l < maxl + 1; ++l) {
         for (int m = 0; m < l + 1; ++m) {
            double multfact = 0.5;
            // multfact = 1.0;
            if (m == 0) {
               multfact = 1.0;
            }
            vec_lm_rad[idxoverall] =
                multfact * (vec_A[idxoverall] - imag1 * vec_B[idxoverall]);
            ++idxoverall;
         }
      }
   }
   {
      // int idxt = 0;
      // for (int idxl = 0; idxl < 3; ++idxl) {
      //    for (int idxm = 0; idxm < idxl + 1; ++idxm) {
      //       std::cout << idxl << " " << idxm << " " << vec_lm_rad[idxt] << "
      //       "
      //                 << vec_A[idxt] << "\n";
      //       ++idxt;
      //    }
      // }
   }

   // physical radius that we will use is the 00 component
   double pi_db = 3.1415926535897932;
   double physicalaverageradius = vec_A[0] / std::sqrt(4.0 * pi_db);
   double referentialradius = physicalaverageradius / length_norm;
   double maxstep = referentialradius * std::min(maxradialstep, 1.0);
   // std::cout << "Ref. rad: " << referentialradius << "\n";

   // std::cout << vec_A[0] << " " << vec_A[1] << " " << vec_A[2] << "\n";
   // std::cout << vec_lm_rad[0] << " " << vec_lm_rad[1] << " " << vec_lm_rad[2]
   //           << "\n";

   // find radial mesh
   node_data = Radial_Tools::RadialMesh(
       referentialradius, ballrad * referentialradius, maxstep, _q);
   // std::cout << referentialradius << " " << ballrad * referentialradius << "
   // "
   //           << maxstep << "\n";

   _num_layers = node_data.NumberOfElements();
   // std::cout << "# layers: " << _num_layers << "\n";
   // std::cout << "Made node_data\n";
   // finding _mat_gaussderiv
   _mat_gaussderiv.resize(_poly_ord + 1, _poly_ord + 1);
   {
      auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
                                                    _q.Points().end());
      for (int idxi = 0; idxi < _poly_ord + 1; ++idxi) {
         for (int idxj = 0; idxj < _poly_ord + 1; ++idxj) {
            _mat_gaussderiv(idxi, idxj) = pleg.Derivative(idxi, _q.X(idxj));
         }
      }
   }
   // std::cout << "Made gauss_derivative\n";
   //////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////
   // default h
   // need to actually fill out h
   // then need to find a
   auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
   // std::cout << "Got spatialsize\n";
   _vec_h = vvvecdb(_num_layers, vvecdb(_q.N(), vecdb(spatialsize, 0.0)));
   std::vector<double> vec_outerradius(spatialsize, 0.0);
   // std::vector<std::complex<double>> vec_outerrad(spatialsize, 0.0);
   // perform SH transform to get physical outer radius
   // std::cout << "Pre-transformation\n";

   _grid.InverseTransformation(lMax, 0, vec_lm_rad, vec_outerradius);
   // _grid.InverseTransformation(lMax, 0, vec_lm_rad_full, vec_outerrad);
   std::vector<std::complex<double>> vec_lm_check(coefficientnumberall, 0.0);

   // for (int idx = 0; idx < spatialsize; ++idx) {
   //    vec_outerradius[idx] = vec_outerrad[idx].real();
   // }
   // max value
   double maxval =
       *std::max_element(vec_outerradius.begin(), vec_outerradius.end());
   double minval =
       *std::min_element(vec_outerradius.begin(), vec_outerradius.end());
   std::cout << "Max radius: " << maxval << ", min element: " << minval << "\n";
   // for (auto &idx : vec_outerradius) {
   //    idx *= 1.0 / length_norm;
   // }
   std::vector<double> vec_houter(spatialsize, 0.0);
   for (int idx = 0; idx < spatialsize; ++idx) {
      vec_houter[idx] = vec_outerradius[idx] / length_norm - referentialradius;
   }
   // _grid.ForwardTransformation(lMax,0,vec_houter,vec_lm_check);
   // std::cout << vec_lm_check[0]/sqrt(3.1415926535*4.0)<<"\n";
   // std::cout << vec_outerradius[0] << " " << vec_outerradius[100] << "\n";

   // fill out h from mapping
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      int laynum = node_data.LayerNumber(idxelem);

      // looping through nodes
      for (int idxnode = 0; idxnode < npoly + 1; ++idxnode) {
         auto radr = node_data.NodeRadius(idxelem, idxnode);
         auto multfact = 1.0;

         // check if within planet
         if (radr > node_data.PlanetRadius()) {
            multfact = (node_data.OuterRadius() - radr) /
                       (node_data.OuterRadius() - node_data.PlanetRadius());
         } else {
            multfact = radr / node_data.PlanetRadius();
         }

         // fill out h
         int idxspatial = 0;
         for (auto it : _grid.CoLatitudes()) {
            for (auto ip : _grid.Longitudes()) {
               _vec_h[idxelem][idxnode][idxspatial] =
                   vec_houter[idxspatial] * multfact;
               ++idxspatial;
            }
         }
      }
   }

   // finding hlm
   //  length of coefficients for YLM
   // auto coefficientnumber =
   //     GSHTrans::GSHIndices<GSHTrans::NonNegative>(lMax, lMax, 0).size();

   // get vec_hlm
   auto vec_hlm = std::vector<std::vector<std::vector<std::complex<double>>>>(
       _num_layers, std::vector<std::vector<std::complex<double>>>(
                        npoly + 1, std::vector<std::complex<double>>(
                                       coefficientnumberall, 0.0)));
   auto indexreal = [](int l, int m) {
      if (m < 0) {
         return (l * (l + 1)) / 2 - m;
      } else {
         return (l * (l + 1)) / 2 + m;
      }
   };
   // std::cout << "Hello pre hlm\n";
   // fill out h from mapping
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      int laynum = node_data.LayerNumber(idxelem);

      // looping through nodes
      for (int idxnode = 0; idxnode < npoly + 1; ++idxnode) {
         std::vector<std::complex<double>> tmp_h(coefficientnumber);
         _grid.ForwardTransformation(lMax, 0, _vec_h[idxelem][idxnode], tmp_h);

         // looping over l,m values
         // std::size_t mycheckidx =
         //     (idxnode + nnode * idxelem) * std::pow(lMax + 1, 2);
         int idx = 0;
         for (int idxl = 0; idxl < lMax + 1; ++idxl) {
            // deal with -ve m:
            for (int idxm = -idxl; idxm < 0; ++idxm) {
               int idxlm = indexreal(idxl, idxm);
               vec_hlm[idxelem][idxnode][idx] =
                   std::pow(-1.0, idxm) * std::conj(tmp_h[idxlm]);
               ++idx;
            }

            //+ve m:
            for (int idxm = 0; idxm < idxl + 1; ++idxm) {
               int idxlm = indexreal(idxl, idxm);
               vec_hlm[idxelem][idxnode][idx] = tmp_h[idxlm];
               ++idx;
            }
         }
      }
   }
   // // std::cout << "Hello pre j\n";
   // finding j
   _vec_j = vvvecdb(_num_layers, vvecdb(_q.N(), vecdb(spatialsize, 0.0)));
   {
      // declaring typenames
      using veccomp = std::vector<double>;
      using vecvech = std::vector<veccomp>;
      auto intsize = _grid.NumberOfLongitudes() * _grid.NumberOfCoLatitudes();
      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
         // max for idxnode loop
         int idxnodemax = npoly + 1;
         int idxnodemin = 0;
         if (idxelem == 0) {
            idxnodemin = 1;
         }

         // components of derivative, ie \nabla h:
         vecvech vec_h0(npoly + 1, veccomp(intsize, 0.0));
         double inv2 = 2.0 / node_data.ElementWidth(idxelem);

         // finding h0 component
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            // fill out h
            int idxspatial = 0;
            for (auto it : _grid.CoLatitudes()) {
               for (auto ip : _grid.Longitudes()) {
                  for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                     vec_h0[idxnode][idxspatial] +=
                         _vec_h[idxelem][idxn][idxspatial] *
                         _mat_gaussderiv(idxn, idxnode) * inv2;
                  }
                  ++idxspatial;
               }
            }
         }

         // std::cout << "Hello 5 \n";

         // find a
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            int idxspatial = 0;
            for (auto it : _grid.CoLatitudes()) {
               for (auto ip : _grid.Longitudes()) {
                  // auto idxuse = idxelem * (npoly + 1) + idxnode;
                  double hdivr;
                  if (idxelem == 0 && idxnode == 0) {
                     hdivr = 0.0;
                  } else {
                     hdivr = _vec_h[idxelem][idxnode][idxspatial] /
                             node_data.NodeRadius(idxelem, idxnode);
                  }
                  _vec_j[idxelem][idxnode][idxspatial] =
                      (1.0 + hdivr) * (1.0 + hdivr) *
                      (1.0 + vec_h0[idxnode][idxspatial]);
                  ++idxspatial;
               }
            }
         }
      }
   }

   // finding F^{-1}:
   {
      Eigen::Matrix3cd mat_F0;
      mat_F0 << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
      _vec_invF =
          vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize, mat_F0)));
      // // various sizes and definitions
      // auto lMax = grid.MaxDegree();
      // auto nMax = grid.MaxUpperIndex();
      auto size =
          GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();   // dof
      auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
      auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
      auto intsize = _grid.NumberOfLongitudes() * _grid.NumberOfCoLatitudes();
      // auto nelem = vec_elemwidth.size();
      // std::size_t nelem =
      // int matlen = nelem * npoly + 1;   // size of matrix
      // auto intsize =
      // grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();

      // declaring typenames
      using veccomp = std::vector<std::complex<double>>;
      using vecvech = std::vector<veccomp>;
      using MATRIX3 = Eigen::Matrix<std::complex<double>, 3, 3>;

      // declaring variables used in calculation
      // vecvech vec_hlm(matlen, veccomp(_size0, 0.0));
      // auto vec_hlm = vecgrid_to_sphdecomp(grid, vec_h, npoly, nelem);
      // auto vec_hcheck = sphdecomp_to_vecgrid(grid, vec_hlm, npoly, nelem);

      // // std::cout << std::endl;
      // auto mat_0 = MATRIX3::Zero();
      // std::vector<std::vector<MATRIX3>> vec_f(
      //     nelem * (npoly + 1), std::vector<MATRIX3>(intsize, mat_0));

      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
         // max for idxpoly loop
         int idxpolymax = npoly + 1;
         int idxpolymin = 0;
         if (idxelem == 0) {
            idxpolymin = 1;
         }

         //    // components of derivative, ie \nabla h:
         vecvech vec_hp(npoly + 1, veccomp(_sizepm, 0.0)),
             vec_hm(npoly + 1, veccomp(_sizepm, 0.0));
         double inv2 = 2.0 / node_data.ElementWidth(idxelem);

         //    // finding h0 component
         vecvech vec_h02(npoly + 1, veccomp(intsize, 0.0));

         // finding h0 component
         for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                  auto idxouter = idxelem * npoly + idxn;
                  vec_h02[idxpoly][idxrad] += _vec_h[idxelem][idxn][idxrad] *
                                              _mat_gaussderiv(idxn, idxpoly) *
                                              inv2;
               }
            }
         }

         for (int idxpoly = idxpolymin; idxpoly < idxpolymax; ++idxpoly) {
            auto idxouter = idxelem * npoly + idxpoly;
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               // omega_l^0
               double Omegal0 = std::sqrt(static_cast<double>(idxl) *
                                          static_cast<double>(idxl + 1) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto idx0 = idxm + idxl * (1 + idxl);
                  auto idxref = idxm + idxl * (1 + idxl) - 1;
                  vec_hp[idxpoly][idxref] =
                      vec_hlm[idxelem][idxpoly][idx0] * Omegal0 /
                      node_data.NodeRadius(idxelem, idxpoly);
                  vec_hm[idxpoly][idxref] = vec_hp[idxpoly][idxref];
               }
            }
         }

         // transform back into spatial
         vecvech vec_nhp(npoly + 1, veccomp(intsize, 0.0)),
             vec_nhm(npoly + 1, veccomp(intsize, 0.0));

         for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
            _grid.InverseTransformation(lMax, 1, vec_hp[idxpoly],
                                        vec_nhp[idxpoly]);
            _grid.InverseTransformation(lMax, -1, vec_hm[idxpoly],
                                        vec_nhm[idxpoly]);
         }

         //    // find a
         for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               auto idxuse = idxelem * (npoly + 1) + idxpoly;
               auto nhm = vec_nhm[idxpoly][idxrad];
               auto nhp = vec_nhp[idxpoly][idxrad];
               // auto nh0 = vec_nh0[idxpoly][idxrad];
               auto nh0 = vec_h02[idxpoly][idxrad];
               std::complex<double> hdivr;

               if (idxuse == 0) {
                  hdivr = 0.0;
               } else {
                  hdivr = _vec_h[idxelem][idxpoly][idxrad] /
                          node_data.NodeRadius(idxelem, idxpoly);
               }

               auto tmp1 = 1.0 / (1.0 + hdivr);
               auto tmp2 = tmp1 / (1.0 + nh0);

               if (std::isnan(std::abs(tmp1))) {
                  std::cout << "idxe: " << idxelem << ", idxp: " << idxpoly
                            << std::endl;
               }

               _vec_invF[idxelem][idxpoly][idxrad](0, 2) = -tmp1;
               _vec_invF[idxelem][idxpoly][idxrad](1, 0) = -nhm * tmp2;
               _vec_invF[idxelem][idxpoly][idxrad](1, 1) =
                   tmp1 - (nh0 - hdivr) * tmp2;
               _vec_invF[idxelem][idxpoly][idxrad](1, 2) = -nhp * tmp2;
               _vec_invF[idxelem][idxpoly][idxrad](2, 0) = -tmp1;

               // if (idxelem == 0 && idxpoly == 0 && idxrad == 0) {
               //    std::cout << "Hello\n\n";
               //    std::cout << tmp1 << " " << tmp2 << " " << nhm << " "
               // << nhp
               //              << " " << hdivr << std::endl;
               //    std::cout << vec_f[0][0] << std::endl;
               //    std::cout << "Hello 2\n\n";
               // }
            }
         }
      }
      // // std::cout << vec_f[0][0] << std::endl;
      // return vec_f;
   }
   // // std::cout << "Hello post j\n";
   {

      Eigen::Matrix3cd mat_a0;
      mat_a0 << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
      _vec_a =
          vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize, mat_a0)));
      // auto mat_0 = MATRIX3::Zero();
      // std::vector<std::vector<MATRIX3>> vec_a(
      //     nelem * (npoly + 1), std::vector<MATRIX3>(intsize, mat_0));
      // declaring typenames
      using veccomp = std::vector<std::complex<double>>;
      using vecvech = std::vector<veccomp>;
      auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
      auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
      auto intsize = _grid.NumberOfLongitudes() * _grid.NumberOfCoLatitudes();

      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
         // max for idxnode loop
         int idxnodemax = npoly + 1;
         int idxnodemin = 0;
         if (idxelem == 0) {
            idxnodemin = 1;
         }

         // components of derivative, ie \nabla h:
         vecvech vec_hp(npoly + 1, veccomp(_sizepm, 0.0)),
             vec_hm(npoly + 1, veccomp(_sizepm, 0.0));
         double inv2 = 2.0 / node_data.ElementWidth(idxelem);

         vecvech vec_h02(npoly + 1, veccomp(intsize, 0.0));

         // finding h0 component
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                  vec_h02[idxnode][idxrad] += _vec_h[idxelem][idxn][idxrad] *
                                              _mat_gaussderiv(idxn, idxnode) *
                                              inv2;
               }
            }
         }

         // need to get vec_hlm, ie spherical harmonic decomposition of h
         for (int idxnode = idxnodemin; idxnode < idxnodemax; ++idxnode) {
            auto idxouter = idxelem * npoly + idxnode;
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               // omega_l^0
               double Omegal0 = std::sqrt(static_cast<double>(idxl) *
                                          static_cast<double>(idxl + 1) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto idx0 = idxm + idxl * (1 + idxl);
                  auto idxref = idx0 - 1;
                  vec_hp[idxnode][idxref] =
                      vec_hlm[idxelem][idxnode][idx0] * Omegal0 /
                      node_data.NodeRadius(idxelem, idxnode);
                  vec_hm[idxnode][idxref] = vec_hp[idxnode][idxref];
               }
            }
         }

         // transform back into spatial
         vecvech vec_nhp(npoly + 1, veccomp(intsize, 0.0)),
             vec_nhm(npoly + 1, veccomp(intsize, 0.0));

         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            // grid.InverseTransformation(lMax, 0, vec_h0[idxnode],
            // vec_nh0[idxnode]);
            _grid.InverseTransformation(lMax, 1, vec_hp[idxnode],
                                        vec_nhp[idxnode]);
            _grid.InverseTransformation(lMax, -1, vec_hm[idxnode],
                                        vec_nhm[idxnode]);
         }

         // find a
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               auto nhm = vec_nhm[idxnode][idxrad];
               auto nhp = vec_nhp[idxnode][idxrad];
               auto nh0 = vec_h02[idxnode][idxrad];
               std::complex<double> hdivr;
               if (idxelem + idxnode == 0) {
                  hdivr = 0.0;
               } else {
                  hdivr = _vec_h[idxelem][idxnode][idxrad] /
                          node_data.NodeRadius(idxelem, idxnode);
               }

               _vec_a[idxelem][idxnode][idxrad](0, 1) = -nhm;
               _vec_a[idxelem][idxnode][idxrad](1, 0) = -nhm;
               _vec_a[idxelem][idxnode][idxrad](1, 2) = -nhp;
               _vec_a[idxelem][idxnode][idxrad](2, 1) = -nhp;
               _vec_a[idxelem][idxnode][idxrad](0, 2) = -(1.0 + nh0);
               _vec_a[idxelem][idxnode][idxrad](2, 0) = -(1.0 + nh0);
               _vec_a[idxelem][idxnode][idxrad](1, 1) =
                   1.0 - nh0 + 2.0 * hdivr +
                   ((nh0 - hdivr) * (nh0 - hdivr) - 2.0 * nhp * nhm) /
                       (1.0 + nh0);
            }
         }
      }
      // _vec_a =
      //     vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize,
      //     mat_a0)));
   }

   // // {
   // //    Eigen::Matrix3cd mat_a0;
   // //    mat_a0 << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
   // //    _vec_a =
   // //        vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize,
   // //        mat_a0)));
   // // }

   // //////////////////////////////////////////////////////////////////////
   // //////////////////////////////////////////////////////////////////////
   // //////////////////////////////////////////////////////////////////////

   // finding the spectral element grid
   _spectral_info = SpectralElementTools::MatrixWeakForm(node_data, _q);

   // std::cout << "Hello\n";
   // density
   // find all model information
   _vec_density.resize(_num_layers);

   double multval = 1.0;

   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      // temporary vectors
      vvecdb tmp_density;   // to contain all values within the element
      int laynum = node_data.LayerNumber(idxelem);

      // loop through polynomial indices
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         double density_1D = 0.0;   // initialise density
         if (node_data.LayerNumber(idxelem) == 0) {
            density_1D = physicaldensity;
         }

         // fill out all values at particular radius
         vecdb tmp_vec_level_density(spatialsize, density_1D);

         // loop over all gridpoints
         int idxspatial = 0;
         for (auto idxt : _grid.CoLatitudes()) {
            for (auto idxp : _grid.Longitudes()) {
               if (randdensity) {
                  double cr = node_data.NodeRadius(idxelem, idxnode);
                  multval = 1.0 + std::sin(4.0 * idxp) * std::sin(3.0 * idxt) *
                                      cr * cr * (1 - cr);
               }
               tmp_vec_level_density[idxspatial] *=
                   _vec_j[idxelem][idxnode][idxspatial] * multval;
               ++idxspatial;
            }
         }
         // }

         tmp_density.push_back(tmp_vec_level_density);
      }

      // assign to the idxelem element
      _vec_density[idxelem] = tmp_density;
   }
};

// simple
//  constructor for homogeneous physical planet with mapping/layers
// importantly one must ensure the mapping is one-to-one.
// the mapping must also be spherical on the referential radius
template <class mapclass>
   requires PlanetaryModel::RadialMappingClass<mapclass>
Density3D::Density3D(double physicalradius, double physicaldensity,
                     double referentialradius, const mapclass &inp_map,
                     const int npoly, const int lMax, double ilength_norm,
                     double itime_norm, double imass_norm, double maxradialstep,
                     double ballrad)
    : length_norm(ilength_norm), mass_norm(imass_norm), time_norm(itime_norm),
      density_norm(imass_norm / std::pow(ilength_norm, 3.0)),
      velocity_norm(ilength_norm / itime_norm),
      acceleration_norm(ilength_norm / std::pow(itime_norm, 2.0)),
      force_norm(imass_norm * ilength_norm / std::pow(itime_norm, 2.0)),
      stress_norm(imass_norm / (std::pow(itime_norm, 2.0) * ilength_norm)),
      inertia_norm(imass_norm * std::pow(ilength_norm, 2.0)),
      gravitational_constant_norm(std::pow(ilength_norm, 3.0) /
                                  (imass_norm * std::pow(itime_norm, 2.0))),
      gravitational_constant(6.6743015 * std::pow(10.0, -11.0) * imass_norm *
                             std::pow(itime_norm, 2.0) /
                             std::pow(ilength_norm, 3.0)),
      potential_norm(std::pow(ilength_norm / itime_norm, 2.0)) {
   _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoly + 1);

   _grid = Grid(lMax, 2);
   // polynomial order
   _poly_ord = _q.N() - 1;

   // find radial mesh
   node_data = Radial_Tools::RadialMesh(
       referentialradius, ballrad * physicalradius, maxradialstep, _q);

   _num_layers = node_data.NumberOfElements();
   std::cout << "The number of elements is: " << _num_layers << "\n";

   // finding _mat_gaussderiv
   _mat_gaussderiv.resize(_poly_ord + 1, _poly_ord + 1);
   {
      auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
                                                    _q.Points().end());
      for (int idxi = 0; idxi < _poly_ord + 1; ++idxi) {
         for (int idxj = 0; idxj < _poly_ord + 1; ++idxj) {
            _mat_gaussderiv(idxi, idxj) = pleg.Derivative(idxi, _q.X(idxj));
         }
      }
   }

   //////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////
   // default h
   // need to actually fill out h
   // then need to find a
   auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
   _vec_h = vvvecdb(_num_layers, vvecdb(_q.N(), vecdb(spatialsize, 0.0)));

   // fill out h from mapping
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      int laynum = node_data.LayerNumber(idxelem);
      std::cout << "Layer number: " << laynum << "\n";
      // looping through nodes
      for (int idxnode = 0; idxnode < npoly + 1; ++idxnode) {
         auto radr = node_data.NodeRadius(idxelem, idxnode);
         auto multfact = 1.0;
         auto raduse = radr;

         // check if within planet
         if (radr > node_data.PlanetRadius()) {
            raduse = node_data.PlanetRadius();
            multfact = (node_data.OuterRadius() - radr) /
                       (node_data.OuterRadius() - node_data.PlanetRadius());
         }

         // fill out h
         int idxspatial = 0;
         for (auto it : _grid.CoLatitudes()) {
            for (auto ip : _grid.Longitudes()) {
               _vec_h[idxelem][idxnode][idxspatial] =
                   inp_map.RadialMapping(laynum)(raduse, it, ip) * multfact;
               ++idxspatial;
            }
         }
      }
   }

   // finding hlm
   //  length of coefficients for YLM
   auto coefficientnumber =
       GSHTrans::GSHIndices<GSHTrans::NonNegative>(lMax, lMax, 0).size();
   auto coefficientnumberall =
       GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
   // get vec_hlm
   auto vec_hlm = std::vector<std::vector<std::vector<std::complex<double>>>>(
       _num_layers, std::vector<std::vector<std::complex<double>>>(
                        npoly + 1, std::vector<std::complex<double>>(
                                       coefficientnumberall, 0.0)));
   auto indexreal = [](int l, int m) {
      if (m < 0) {
         return (l * (l + 1)) / 2 - m;
      } else {
         return (l * (l + 1)) / 2 + m;
      }
   };
   // std::cout << "Hello pre hlm\n";
   // fill out h from mapping
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      int laynum = node_data.LayerNumber(idxelem);

      // looping through nodes
      for (int idxnode = 0; idxnode < npoly + 1; ++idxnode) {
         std::vector<std::complex<double>> tmp_h(coefficientnumber);
         _grid.ForwardTransformation(lMax, 0, _vec_h[idxelem][idxnode], tmp_h);

         // looping over l,m values
         // std::size_t mycheckidx =
         //     (idxnode + nnode * idxelem) * std::pow(lMax + 1, 2);
         int idx = 0;
         for (int idxl = 0; idxl < lMax + 1; ++idxl) {
            // deal with -ve m:
            for (int idxm = -idxl; idxm < 0; ++idxm) {
               int idxlm = indexreal(idxl, idxm);
               vec_hlm[idxelem][idxnode][idx] =
                   std::pow(-1.0, idxm) * std::conj(tmp_h[idxlm]);
               ++idx;
            }

            //+ve m:
            for (int idxm = 0; idxm < idxl + 1; ++idxm) {
               int idxlm = indexreal(idxl, idxm);
               vec_hlm[idxelem][idxnode][idx] = tmp_h[idxlm];
               ++idx;
            }
         }
      }
   }
   // std::cout << "Hello pre j\n";
   // finding j
   _vec_j = vvvecdb(_num_layers, vvecdb(_q.N(), vecdb(spatialsize, 0.0)));
   {
      // declaring typenames
      using veccomp = std::vector<double>;
      using vecvech = std::vector<veccomp>;
      auto intsize = _grid.NumberOfLongitudes() * _grid.NumberOfCoLatitudes();
      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
         // max for idxnode loop
         int idxnodemax = npoly + 1;
         int idxnodemin = 0;
         if (idxelem == 0) {
            idxnodemin = 1;
         }

         // components of derivative, ie \nabla h:
         vecvech vec_h0(npoly + 1, veccomp(intsize, 0.0));
         double inv2 = 2.0 / node_data.ElementWidth(idxelem);

         // finding h0 component
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            // fill out h
            int idxspatial = 0;
            for (auto it : _grid.CoLatitudes()) {
               for (auto ip : _grid.Longitudes()) {
                  for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                     vec_h0[idxnode][idxspatial] +=
                         _vec_h[idxelem][idxn][idxspatial] *
                         _mat_gaussderiv(idxn, idxnode) * inv2;
                  }
                  ++idxspatial;
               }
            }
         }

         // std::cout << "Hello 5 \n";

         // find a
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            int idxspatial = 0;
            for (auto it : _grid.CoLatitudes()) {
               for (auto ip : _grid.Longitudes()) {
                  // auto idxuse = idxelem * (npoly + 1) + idxnode;
                  double hdivr;
                  if (idxelem == 0 && idxnode == 0) {
                     hdivr = 0.0;
                  } else {
                     hdivr = _vec_h[idxelem][idxnode][idxspatial] /
                             node_data.NodeRadius(idxelem, idxnode);
                  }
                  _vec_j[idxelem][idxnode][idxspatial] =
                      (1.0 + hdivr) * (1.0 + hdivr) *
                      (1.0 + vec_h0[idxnode][idxspatial]);
                  ++idxspatial;
               }
            }
         }
      }
   }

   // finding F^{-1}:
   {
      Eigen::Matrix3cd mat_F0;
      mat_F0 << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
      _vec_invF =
          vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize, mat_F0)));
      // // various sizes and definitions
      // auto lMax = grid.MaxDegree();
      // auto nMax = grid.MaxUpperIndex();
      auto size =
          GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();   // dof
      auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
      auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
      auto intsize = _grid.NumberOfLongitudes() * _grid.NumberOfCoLatitudes();
      // auto nelem = vec_elemwidth.size();
      // std::size_t nelem =
      // int matlen = nelem * npoly + 1;   // size of matrix
      // auto intsize =
      // grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();

      // declaring typenames
      using veccomp = std::vector<std::complex<double>>;
      using vecvech = std::vector<veccomp>;
      using MATRIX3 = Eigen::Matrix<std::complex<double>, 3, 3>;

      // declaring variables used in calculation
      // vecvech vec_hlm(matlen, veccomp(_size0, 0.0));
      // auto vec_hlm = vecgrid_to_sphdecomp(grid, vec_h, npoly, nelem);
      // auto vec_hcheck = sphdecomp_to_vecgrid(grid, vec_hlm, npoly, nelem);

      // // std::cout << std::endl;
      // auto mat_0 = MATRIX3::Zero();
      // std::vector<std::vector<MATRIX3>> vec_f(
      //     nelem * (npoly + 1), std::vector<MATRIX3>(intsize, mat_0));

      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
         // max for idxpoly loop
         int idxpolymax = npoly + 1;
         int idxpolymin = 0;
         if (idxelem == 0) {
            idxpolymin = 1;
         }

         //    // components of derivative, ie \nabla h:
         vecvech vec_hp(npoly + 1, veccomp(_sizepm, 0.0)),
             vec_hm(npoly + 1, veccomp(_sizepm, 0.0));
         double inv2 = 2.0 / node_data.ElementWidth(idxelem);

         //    // finding h0 component
         vecvech vec_h02(npoly + 1, veccomp(intsize, 0.0));

         // finding h0 component
         for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                  auto idxouter = idxelem * npoly + idxn;
                  vec_h02[idxpoly][idxrad] += _vec_h[idxelem][idxn][idxrad] *
                                              _mat_gaussderiv(idxn, idxpoly) *
                                              inv2;
               }
            }
         }

         for (int idxpoly = idxpolymin; idxpoly < idxpolymax; ++idxpoly) {
            auto idxouter = idxelem * npoly + idxpoly;
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               // omega_l^0
               double Omegal0 = std::sqrt(static_cast<double>(idxl) *
                                          static_cast<double>(idxl + 1) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto idx0 = idxm + idxl * (1 + idxl);
                  auto idxref = idxm + idxl * (1 + idxl) - 1;
                  vec_hp[idxpoly][idxref] =
                      vec_hlm[idxelem][idxpoly][idx0] * Omegal0 /
                      node_data.NodeRadius(idxelem, idxpoly);
                  vec_hm[idxpoly][idxref] = vec_hp[idxpoly][idxref];
               }
            }
         }

         // transform back into spatial
         vecvech vec_nhp(npoly + 1, veccomp(intsize, 0.0)),
             vec_nhm(npoly + 1, veccomp(intsize, 0.0));

         for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
            _grid.InverseTransformation(lMax, 1, vec_hp[idxpoly],
                                        vec_nhp[idxpoly]);
            _grid.InverseTransformation(lMax, -1, vec_hm[idxpoly],
                                        vec_nhm[idxpoly]);
         }

         //    // find a
         for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               auto idxuse = idxelem * (npoly + 1) + idxpoly;
               auto nhm = vec_nhm[idxpoly][idxrad];
               auto nhp = vec_nhp[idxpoly][idxrad];
               // auto nh0 = vec_nh0[idxpoly][idxrad];
               auto nh0 = vec_h02[idxpoly][idxrad];
               std::complex<double> hdivr;

               if (idxuse == 0) {
                  hdivr = 0.0;
               } else {
                  hdivr = _vec_h[idxelem][idxpoly][idxrad] /
                          node_data.NodeRadius(idxelem, idxpoly);
               }

               auto tmp1 = 1.0 / (1.0 + hdivr);
               auto tmp2 = tmp1 / (1.0 + nh0);

               if (std::isnan(std::abs(tmp1))) {
                  std::cout << "idxe: " << idxelem << ", idxp: " << idxpoly
                            << std::endl;
               }

               _vec_invF[idxelem][idxpoly][idxrad](0, 2) = -tmp1;
               _vec_invF[idxelem][idxpoly][idxrad](1, 0) = -nhm * tmp2;
               _vec_invF[idxelem][idxpoly][idxrad](1, 1) =
                   tmp1 - (nh0 - hdivr) * tmp2;
               _vec_invF[idxelem][idxpoly][idxrad](1, 2) = -nhp * tmp2;
               _vec_invF[idxelem][idxpoly][idxrad](2, 0) = -tmp1;

               // if (idxelem == 0 && idxpoly == 0 && idxrad == 0) {
               //    std::cout << "Hello\n\n";
               //    std::cout << tmp1 << " " << tmp2 << " " << nhm << " "
               // << nhp
               //              << " " << hdivr << std::endl;
               //    std::cout << vec_f[0][0] << std::endl;
               //    std::cout << "Hello 2\n\n";
               // }
            }
         }
      }
      // // std::cout << vec_f[0][0] << std::endl;
      // return vec_f;
   }
   // std::cout << "Hello post j\n";
   {

      Eigen::Matrix3cd mat_a0;
      mat_a0 << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
      _vec_a =
          vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize, mat_a0)));
      // auto mat_0 = MATRIX3::Zero();
      // std::vector<std::vector<MATRIX3>> vec_a(
      //     nelem * (npoly + 1), std::vector<MATRIX3>(intsize, mat_0));
      // declaring typenames
      using veccomp = std::vector<std::complex<double>>;
      using vecvech = std::vector<veccomp>;
      auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
      auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
      auto intsize = _grid.NumberOfLongitudes() * _grid.NumberOfCoLatitudes();

      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
         // max for idxnode loop
         int idxnodemax = npoly + 1;
         int idxnodemin = 0;
         if (idxelem == 0) {
            idxnodemin = 1;
         }

         // components of derivative, ie \nabla h:
         vecvech vec_hp(npoly + 1, veccomp(_sizepm, 0.0)),
             vec_hm(npoly + 1, veccomp(_sizepm, 0.0));
         double inv2 = 2.0 / node_data.ElementWidth(idxelem);

         vecvech vec_h02(npoly + 1, veccomp(intsize, 0.0));

         // finding h0 component
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                  vec_h02[idxnode][idxrad] += _vec_h[idxelem][idxn][idxrad] *
                                              _mat_gaussderiv(idxn, idxnode) *
                                              inv2;
               }
            }
         }

         // need to get vec_hlm, ie spherical harmonic decomposition of h
         for (int idxnode = idxnodemin; idxnode < idxnodemax; ++idxnode) {
            auto idxouter = idxelem * npoly + idxnode;
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               // omega_l^0
               double Omegal0 = std::sqrt(static_cast<double>(idxl) *
                                          static_cast<double>(idxl + 1) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto idx0 = idxm + idxl * (1 + idxl);
                  auto idxref = idx0 - 1;
                  vec_hp[idxnode][idxref] =
                      vec_hlm[idxelem][idxnode][idx0] * Omegal0 /
                      node_data.NodeRadius(idxelem, idxnode);
                  vec_hm[idxnode][idxref] = vec_hp[idxnode][idxref];
               }
            }
         }

         // transform back into spatial
         vecvech vec_nhp(npoly + 1, veccomp(intsize, 0.0)),
             vec_nhm(npoly + 1, veccomp(intsize, 0.0));

         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            // grid.InverseTransformation(lMax, 0, vec_h0[idxnode],
            // vec_nh0[idxnode]);
            _grid.InverseTransformation(lMax, 1, vec_hp[idxnode],
                                        vec_nhp[idxnode]);
            _grid.InverseTransformation(lMax, -1, vec_hm[idxnode],
                                        vec_nhm[idxnode]);
         }

         // find a
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               auto nhm = vec_nhm[idxnode][idxrad];
               auto nhp = vec_nhp[idxnode][idxrad];
               auto nh0 = vec_h02[idxnode][idxrad];
               std::complex<double> hdivr;
               if (idxelem + idxnode == 0) {
                  hdivr = 0.0;
               } else {
                  hdivr = _vec_h[idxelem][idxnode][idxrad] /
                          node_data.NodeRadius(idxelem, idxnode);
               }

               _vec_a[idxelem][idxnode][idxrad](0, 1) = -nhm;
               _vec_a[idxelem][idxnode][idxrad](1, 0) = -nhm;
               _vec_a[idxelem][idxnode][idxrad](1, 2) = -nhp;
               _vec_a[idxelem][idxnode][idxrad](2, 1) = -nhp;
               _vec_a[idxelem][idxnode][idxrad](0, 2) = -(1.0 + nh0);
               _vec_a[idxelem][idxnode][idxrad](2, 0) = -(1.0 + nh0);
               _vec_a[idxelem][idxnode][idxrad](1, 1) =
                   1.0 - nh0 + 2.0 * hdivr +
                   ((nh0 - hdivr) * (nh0 - hdivr) - 2.0 * nhp * nhm) /
                       (1.0 + nh0);
            }
         }
      }
      // _vec_a =
      //     vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize,
      //     mat_a0)));
   }

   // {
   //    Eigen::Matrix3cd mat_a0;
   //    mat_a0 << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
   //    _vec_a =
   //        vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize,
   //        mat_a0)));
   // }

   //////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////

   // finding the spectral element grid
   _spectral_info = SpectralElementTools::MatrixWeakForm(node_data, _q);

   // std::cout << "Hello\n";
   // density
   // find all model information
   _vec_density.resize(_num_layers);

   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      // temporary vectors
      vvecdb tmp_density;   // to contain all values within the element
      int laynum = node_data.LayerNumber(idxelem);

      // loop through polynomial indices
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         double density_1D = 0.0;   // initialise density
         if (node_data.LayerNumber(idxelem) == 0) {
            density_1D = physicaldensity;
         }

         // fill out all values at particular radius
         vecdb tmp_vec_level_density(spatialsize, density_1D);

         // loop over all gridpoints
         int idxspatial = 0;
         for (auto idxt : _grid.CoLatitudes()) {
            for (auto idxp : _grid.Longitudes()) {
               tmp_vec_level_density[idxspatial] *=
                   _vec_j[idxelem][idxnode][idxspatial];
               ++idxspatial;
            }
         }
         // }

         tmp_density.push_back(tmp_vec_level_density);
      }

      // assign to the idxelem element
      _vec_density[idxelem] = tmp_density;
   }
};

// constructs a model with referential density given by the 1d input and the
// tomography model, with a mapping given by the map class
template <class model, class tomomodel, class mapclass>
   requires PlanetaryModel::BasicSphericalDensityModel<model, int, double> &&
                PlanetaryModel::TomographyModel<tomomodel> &&
                PlanetaryModel::RadialMappingClass<mapclass>
Density3D::Density3D(const model &inp_model, const tomomodel &tomo_model,
                     const mapclass &inp_map, const int npoly, const int lMax,
                     double max_radial_step, double maxrad)
    : length_norm(inp_model.LengthNorm()), mass_norm(inp_model.MassNorm()),
      time_norm(inp_model.TimeNorm()),
      density_norm(inp_model.MassNorm() /
                   std::pow(inp_model.LengthNorm(), 3.0)),
      velocity_norm(inp_model.LengthNorm() / inp_model.TimeNorm()),
      acceleration_norm(inp_model.LengthNorm() /
                        std::pow(inp_model.TimeNorm(), 2.0)),
      force_norm(inp_model.MassNorm() * inp_model.LengthNorm() /
                 std::pow(inp_model.TimeNorm(), 2.0)),
      stress_norm(inp_model.MassNorm() / (std::pow(inp_model.TimeNorm(), 2.0) *
                                          inp_model.LengthNorm())),
      inertia_norm(inp_model.MassNorm() *
                   std::pow(inp_model.LengthNorm(), 2.0)),
      gravitational_constant_norm(
          std::pow(inp_model.LengthNorm(), 3.0) /
          (inp_model.MassNorm() * std::pow(inp_model.TimeNorm(), 2.0))),
      gravitational_constant(6.6743015 * std::pow(10.0, -11.0) *
                             inp_model.MassNorm() *
                             std::pow(inp_model.TimeNorm(), 2.0) /
                             std::pow(inp_model.LengthNorm(), 3.0)),
      potential_norm(
          std::pow(inp_model.LengthNorm() / inp_model.TimeNorm(), 2.0)) {

   _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoly + 1);

   _grid = Grid(lMax, 2);
   // polynomial order
   _poly_ord = _q.N() - 1;

   // find radial mesh
   node_data = Radial_Tools::RadialMesh(inp_model, _q, max_radial_step, maxrad);
   _num_layers = node_data.NumberOfElements();

   // finding _mat_gaussderiv
   _mat_gaussderiv.resize(_poly_ord + 1, _poly_ord + 1);
   {
      auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
                                                    _q.Points().end());
      for (int idxi = 0; idxi < _poly_ord + 1; ++idxi) {
         for (int idxj = 0; idxj < _poly_ord + 1; ++idxj) {
            _mat_gaussderiv(idxi, idxj) = pleg.Derivative(idxi, _q.X(idxj));
         }
      }
   }

   //////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////
   // default h
   // need to actually fill out h
   // then need to find a
   auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
   _vec_h = vvvecdb(_num_layers, vvecdb(_q.N(), vecdb(spatialsize, 0.0)));
   _vec_j = vvvecdb(_num_layers, vvecdb(_q.N(), vecdb(spatialsize, 1.0)));
   // fill out h from mapping
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      int laynum = node_data.LayerNumber(idxelem);

      // looping through nodes
      for (int idxnode = 0; idxnode < npoly + 1; ++idxnode) {
         auto radr = node_data.NodeRadius(idxelem, idxnode);
         auto multfact = 1.0;
         auto raduse = radr;

         // check if within planet
         if (radr > node_data.PlanetRadius()) {
            raduse = node_data.PlanetRadius();
            multfact = (node_data.OuterRadius() - radr) /
                       (node_data.OuterRadius() - node_data.PlanetRadius());
         }

         // fill out h
         int idxspatial = 0;
         for (auto it : _grid.CoLatitudes()) {
            for (auto ip : _grid.Longitudes()) {
               _vec_h[idxelem][idxnode][idxspatial++] =
                   inp_map.RadialMapping(laynum)(raduse, it, ip) * multfact;
            }
         }
      }
   }

   // finding hlm
   //  length of coefficients for YLM
   auto coefficientnumber =
       GSHTrans::GSHIndices<GSHTrans::NonNegative>(lMax, lMax, 0).size();
   auto coefficientnumberall =
       GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
   // get vec_hlm
   auto vec_hlm = std::vector<std::vector<std::vector<std::complex<double>>>>(
       _num_layers, std::vector<std::vector<std::complex<double>>>(
                        npoly + 1, std::vector<std::complex<double>>(
                                       coefficientnumberall, 0.0)));
   auto indexreal = [](int l, int m) {
      if (m < 0) {
         return (l * (l + 1)) / 2 - m;
      } else {
         return (l * (l + 1)) / 2 + m;
      }
   };

   // std::cout << "Finding hlm\n";
   // fill out h from mapping
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      int laynum = node_data.LayerNumber(idxelem);

      // looping through nodes
      for (int idxnode = 0; idxnode < npoly + 1; ++idxnode) {
         std::vector<std::complex<double>> tmp_h(coefficientnumber);
         _grid.ForwardTransformation(lMax, 0, _vec_h[idxelem][idxnode], tmp_h);

         // looping over l,m values
         // std::size_t mycheckidx =
         //     (idxnode + nnode * idxelem) * std::pow(lMax + 1, 2);
         int idx = 0;
         for (int idxl = 0; idxl < lMax + 1; ++idxl) {
            // deal with -ve m:
            for (int idxm = -idxl; idxm < 0; ++idxm) {
               int idxlm = indexreal(idxl, idxm);
               vec_hlm[idxelem][idxnode][idx] =
                   std::pow(-1.0, idxm) * std::conj(tmp_h[idxlm]);
               ++idx;
            }

            //+ve m:
            for (int idxm = 0; idxm < idxl + 1; ++idxm) {
               int idxlm = indexreal(idxl, idxm);
               vec_hlm[idxelem][idxnode][idx] = tmp_h[idxlm];
               ++idx;
            }
         }
      }
   }
   // std::cout << "Finding j\n";

   // finding j

   {
      // declaring typenames
      using veccomp = std::vector<double>;
      using vecvech = std::vector<veccomp>;
      auto intsize = _grid.NumberOfLongitudes() * _grid.NumberOfCoLatitudes();
      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
         // max for idxnode loop
         int idxnodemax = npoly + 1;
         int idxnodemin = 0;
         if (idxelem == 0) {
            idxnodemin = 1;
         }

         // components of derivative, ie \nabla h:
         vecvech vec_h0(npoly + 1, veccomp(intsize, 0.0));
         double inv2 = 2.0 / node_data.ElementWidth(idxelem);

         // finding h0 component
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            // fill out h
            int idxspatial = 0;
            for (auto it : _grid.CoLatitudes()) {
               for (auto ip : _grid.Longitudes()) {
                  for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                     vec_h0[idxnode][idxspatial] +=
                         _vec_h[idxelem][idxnode][idxspatial] *
                         _mat_gaussderiv(idxn, idxnode) * inv2;
                  }
                  ++idxspatial;
               }
            }
         }

         // find a
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            int idxspatial = 0;
            for (auto it : _grid.CoLatitudes()) {
               for (auto ip : _grid.Longitudes()) {
                  // auto idxuse = idxelem * (npoly + 1) + idxnode;
                  // std::cout << "FIRST\n";
                  double hdivr;
                  if (idxelem == 0 && idxnode == 0) {
                     hdivr = 0.0;
                  } else {
                     // std::cout << idxelem << " " << idxnode << " " <<
                     // idxspatial
                     //           << "\n";
                     hdivr = _vec_h[idxelem][idxnode][idxspatial] /
                             node_data.NodeRadius(idxelem, idxnode);
                  }
                  // std::cout << idxelem << " " << idxnode << " " << idxspatial
                  //           << "\n";
                  // std::cout << node_data.NodeRadius(idxelem, idxnode) <<
                  // "\n"; std::cout << vec_h0[idxnode][idxspatial] << "\n";
                  _vec_j[idxelem][idxnode][idxspatial] =
                      (1.0 + hdivr) * (1.0 + hdivr) *
                      (1.0 + vec_h0[idxnode][idxspatial]);
                  // std::cout << "Segmentational dissipation\n";
                  ++idxspatial;
               }
            }
         }
      }
   }

   // finding F^{-1}:
   {
      Eigen::Matrix3cd mat_F0;
      mat_F0 << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
      _vec_invF =
          vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize, mat_F0)));
      // // various sizes and definitions
      // auto lMax = grid.MaxDegree();
      // auto nMax = grid.MaxUpperIndex();
      auto size =
          GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();   // dof
      auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
      auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
      auto intsize = _grid.NumberOfLongitudes() * _grid.NumberOfCoLatitudes();
      // auto nelem = vec_elemwidth.size();
      // std::size_t nelem =
      // int matlen = nelem * npoly + 1;   // size of matrix
      // auto intsize =
      // grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();

      // declaring typenames
      using veccomp = std::vector<std::complex<double>>;
      using vecvech = std::vector<veccomp>;
      using MATRIX3 = Eigen::Matrix<std::complex<double>, 3, 3>;

      // declaring variables used in calculation
      // vecvech vec_hlm(matlen, veccomp(_size0, 0.0));
      // auto vec_hlm = vecgrid_to_sphdecomp(grid, vec_h, npoly, nelem);
      // auto vec_hcheck = sphdecomp_to_vecgrid(grid, vec_hlm, npoly, nelem);

      // // std::cout << std::endl;
      // auto mat_0 = MATRIX3::Zero();
      // std::vector<std::vector<MATRIX3>> vec_f(
      //     nelem * (npoly + 1), std::vector<MATRIX3>(intsize, mat_0));

      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
         // max for idxpoly loop
         int idxpolymax = npoly + 1;
         int idxpolymin = 0;
         if (idxelem == 0) {
            idxpolymin = 1;
         }

         //    // components of derivative, ie \nabla h:
         vecvech vec_hp(npoly + 1, veccomp(_sizepm, 0.0)),
             vec_hm(npoly + 1, veccomp(_sizepm, 0.0));
         double inv2 = 2.0 / node_data.ElementWidth(idxelem);

         //    // finding h0 component
         vecvech vec_h02(npoly + 1, veccomp(intsize, 0.0));

         // finding h0 component
         for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                  auto idxouter = idxelem * npoly + idxn;
                  vec_h02[idxpoly][idxrad] += _vec_h[idxelem][idxn][idxrad] *
                                              _mat_gaussderiv(idxn, idxpoly) *
                                              inv2;
               }
            }
         }

         for (int idxpoly = idxpolymin; idxpoly < idxpolymax; ++idxpoly) {
            auto idxouter = idxelem * npoly + idxpoly;
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               // omega_l^0
               double Omegal0 = std::sqrt(static_cast<double>(idxl) *
                                          static_cast<double>(idxl + 1) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto idx0 = idxm + idxl * (1 + idxl);
                  auto idxref = idxm + idxl * (1 + idxl) - 1;
                  vec_hp[idxpoly][idxref] =
                      vec_hlm[idxelem][idxpoly][idx0] * Omegal0 /
                      node_data.NodeRadius(idxelem, idxpoly);
                  vec_hm[idxpoly][idxref] = vec_hp[idxpoly][idxref];
               }
            }
         }

         // transform back into spatial
         vecvech vec_nhp(npoly + 1, veccomp(intsize, 0.0)),
             vec_nhm(npoly + 1, veccomp(intsize, 0.0));

         for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
            _grid.InverseTransformation(lMax, 1, vec_hp[idxpoly],
                                        vec_nhp[idxpoly]);
            _grid.InverseTransformation(lMax, -1, vec_hm[idxpoly],
                                        vec_nhm[idxpoly]);
         }

         //    // find a
         for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               auto idxuse = idxelem * (npoly + 1) + idxpoly;
               auto nhm = vec_nhm[idxpoly][idxrad];
               auto nhp = vec_nhp[idxpoly][idxrad];
               // auto nh0 = vec_nh0[idxpoly][idxrad];
               auto nh0 = vec_h02[idxpoly][idxrad];
               std::complex<double> hdivr;

               if (idxuse == 0) {
                  hdivr = 0.0;
               } else {
                  hdivr = _vec_h[idxelem][idxpoly][idxrad] /
                          node_data.NodeRadius(idxelem, idxpoly);
               }

               auto tmp1 = 1.0 / (1.0 + hdivr);
               auto tmp2 = tmp1 / (1.0 + nh0);

               if (std::isnan(std::abs(tmp1))) {
                  std::cout << "idxe: " << idxelem << ", idxp: " << idxpoly
                            << std::endl;
               }

               _vec_invF[idxelem][idxpoly][idxrad](0, 2) = -tmp1;
               _vec_invF[idxelem][idxpoly][idxrad](1, 0) = -nhm * tmp2;
               _vec_invF[idxelem][idxpoly][idxrad](1, 1) =
                   tmp1 - (nh0 - hdivr) * tmp2;
               _vec_invF[idxelem][idxpoly][idxrad](1, 2) = -nhp * tmp2;
               _vec_invF[idxelem][idxpoly][idxrad](2, 0) = -tmp1;

               // if (idxelem == 0 && idxpoly == 0 && idxrad == 0) {
               //    std::cout << "Hello\n\n";
               //    std::cout << tmp1 << " " << tmp2 << " " << nhm << " "
               // << nhp
               //              << " " << hdivr << std::endl;
               //    std::cout << vec_f[0][0] << std::endl;
               //    std::cout << "Hello 2\n\n";
               // }
            }
         }
      }
      // // std::cout << vec_f[0][0] << std::endl;
      // return vec_f;
   }
   std::cout << "Finding a\n";
   {

      Eigen::Matrix3cd mat_a0;
      mat_a0 << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
      _vec_a =
          vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize, mat_a0)));
      // auto mat_0 = MATRIX3::Zero();
      // std::vector<std::vector<MATRIX3>> vec_a(
      //     nelem * (npoly + 1), std::vector<MATRIX3>(intsize, mat_0));
      // declaring typenames
      using veccomp = std::vector<std::complex<double>>;
      using vecvech = std::vector<veccomp>;
      auto _size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
      auto _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
      auto intsize = _grid.NumberOfLongitudes() * _grid.NumberOfCoLatitudes();

      for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
         // max for idxnode loop
         int idxnodemax = npoly + 1;
         int idxnodemin = 0;
         if (idxelem == 0) {
            idxnodemin = 1;
         }

         // components of derivative, ie \nabla h:
         vecvech vec_hp(npoly + 1, veccomp(_sizepm, 0.0)),
             vec_hm(npoly + 1, veccomp(_sizepm, 0.0));
         double inv2 = 2.0 / node_data.ElementWidth(idxelem);

         vecvech vec_h02(npoly + 1, veccomp(intsize, 0.0));

         // finding h0 component
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               for (int idxn = 0; idxn < npoly + 1; ++idxn) {
                  vec_h02[idxnode][idxrad] += _vec_h[idxelem][idxn][idxrad] *
                                              _mat_gaussderiv(idxn, idxnode) *
                                              inv2;
               }
            }
         }

         // need to get vec_hlm, ie spherical harmonic decomposition of h
         for (int idxnode = idxnodemin; idxnode < idxnodemax; ++idxnode) {
            auto idxouter = idxelem * npoly + idxnode;
            for (int idxl = 1; idxl < lMax + 1; ++idxl) {
               // omega_l^0
               double Omegal0 = std::sqrt(static_cast<double>(idxl) *
                                          static_cast<double>(idxl + 1) / 2.0);
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  auto idx0 = idxm + idxl * (1 + idxl);
                  auto idxref = idx0 - 1;
                  vec_hp[idxnode][idxref] =
                      vec_hlm[idxelem][idxnode][idx0] * Omegal0 /
                      node_data.NodeRadius(idxelem, idxnode);
                  vec_hm[idxnode][idxref] = vec_hp[idxnode][idxref];
               }
            }
         }

         // transform back into spatial
         vecvech vec_nhp(npoly + 1, veccomp(intsize, 0.0)),
             vec_nhm(npoly + 1, veccomp(intsize, 0.0));

         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            // grid.InverseTransformation(lMax, 0, vec_h0[idxnode],
            // vec_nh0[idxnode]);
            _grid.InverseTransformation(lMax, 1, vec_hp[idxnode],
                                        vec_nhp[idxnode]);
            _grid.InverseTransformation(lMax, -1, vec_hm[idxnode],
                                        vec_nhm[idxnode]);
         }

         // find a
         for (int idxnode = 0; idxnode < idxnodemax; ++idxnode) {
            for (int idxrad = 0; idxrad < intsize; ++idxrad) {
               auto nhm = vec_nhm[idxnode][idxrad];
               auto nhp = vec_nhp[idxnode][idxrad];
               auto nh0 = vec_h02[idxnode][idxrad];
               std::complex<double> hdivr;
               if (idxelem + idxnode == 0) {
                  hdivr = 0.0;
               } else {
                  hdivr = _vec_h[idxelem][idxnode][idxrad] /
                          node_data.NodeRadius(idxelem, idxnode);
               }

               _vec_a[idxelem][idxnode][idxrad](0, 1) = -nhm;
               _vec_a[idxelem][idxnode][idxrad](1, 0) = -nhm;
               _vec_a[idxelem][idxnode][idxrad](1, 2) = -nhp;
               _vec_a[idxelem][idxnode][idxrad](2, 1) = -nhp;
               _vec_a[idxelem][idxnode][idxrad](0, 2) = -(1.0 + nh0);
               _vec_a[idxelem][idxnode][idxrad](2, 0) = -(1.0 + nh0);
               _vec_a[idxelem][idxnode][idxrad](1, 1) =
                   1.0 - nh0 + 2.0 * hdivr +
                   ((nh0 - hdivr) * (nh0 - hdivr) - 2.0 * nhp * nhm) /
                       (1.0 + nh0);
            }
         }
      }
      for (int idxrad = 0; idxrad < intsize; ++idxrad) {

         _vec_a[0][0][idxrad](0, 0) = 0.0;
         _vec_a[0][0][idxrad](0, 1) = 0.0;
         _vec_a[0][0][idxrad](0, 2) = -1.0;
         _vec_a[0][0][idxrad](1, 0) = 0.0;
         _vec_a[0][0][idxrad](1, 1) = 1.0;
         _vec_a[0][0][idxrad](1, 2) = 0.0;
         _vec_a[0][0][idxrad](2, 0) = -1.0;
         _vec_a[0][0][idxrad](2, 1) = 0.0;
         _vec_a[0][0][idxrad](2, 2) = 0.0;
      }
   }

   // {
   //    Eigen::Matrix3cd mat_a0;
   //    mat_a0 << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
   //    _vec_a =
   //        vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize,
   //        mat_a0)));
   // }

   //////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////

   // finding the spectral element grid
   _spectral_info = SpectralElementTools::MatrixWeakForm(node_data, _q);

   // lambda to convert between radius and depth (as in depth needed for
   // tomography model)
   auto rad_to_depth = [](double radius, double outerradius,
                          double lengthnorm) {
      double depth = outerradius - radius;
      depth *= lengthnorm / 1000.0;
      return depth;
   };
   auto colatitude_to_latitude = [](double colatitude) {
      double pi_db = 3.1415926535897932384626433;
      double latitude = pi_db / 2.0 - colatitude;
      latitude *= 180.0 / pi_db;
      return latitude;
   };
   auto longitude_to_degrees = [](double longitude) {
      double multfact = 180.0 / 3.1415926535897932384626433;
      return multfact * longitude;
   };

   // density
   // find all model information
   _vec_density.resize(_num_layers);

   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      // temporary vectors
      vvecdb tmp_density;   // to contain all values within the element
      int laynum = node_data.LayerNumber(idxelem);

      // loop through polynomial indices
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         double rad_current =
             node_data.NodeRadius(idxelem, idxnode);   // current radius
         // if (!isrefmodel) {
         //    rad_current += _vec_h[idxelem][idxnode][0];
         // }
         double density_1D = 0.0;   // initialise density

         // set density iff inside input model
         if (laynum < inp_model.NumberOfLayers()) {
            density_1D = inp_model.Density(laynum)(rad_current);
         }

         // fill out all values at particular radius
         vecdb tmp_vec_level_density(spatialsize, density_1D);
         auto depth = rad_to_depth(rad_current, node_data.OuterRadius(),
                                   inp_model.LengthNorm());

         // check within tomography model
         if (depth > tomo_model.GetDepths()[0] &&
             depth < tomo_model.GetDepths().back()) {

            // loop over all gridpoints
            int idxspatial = 0;
            for (auto idxt : _grid.CoLatitudes()) {
               auto latitude = colatitude_to_latitude(idxt);
               for (auto idxp : _grid.Longitudes()) {
                  auto longitude = longitude_to_degrees(idxp);
                  tmp_vec_level_density[idxspatial] *=
                      (1.0 + 0.005 * tomo_model.GetValueAt(depth, longitude,
                                                           latitude));
                  ++idxspatial;
               }
            }
         }

         tmp_density.push_back(tmp_vec_level_density);
      }

      // assign to the idxelem element
      _vec_density[idxelem] = tmp_density;
   }
};

// with just npoly and lMax
//  with 1D and 3D input models that satisfies requirements
template <class model, class tomomodel>
   requires PlanetaryModel::BasicSphericalDensityModel<model, int, double> &&
                PlanetaryModel::TomographyModel<tomomodel>
Density3D::Density3D(const model &inp_model, const tomomodel &tomo_model,
                     const int npoly, const int lMax, double max_radial_step,
                     double maxrad)
    : length_norm(inp_model.LengthNorm()), mass_norm(inp_model.MassNorm()),
      time_norm(inp_model.TimeNorm()),
      density_norm(inp_model.MassNorm() /
                   std::pow(inp_model.LengthNorm(), 3.0)),
      velocity_norm(inp_model.LengthNorm() / inp_model.TimeNorm()),
      acceleration_norm(inp_model.LengthNorm() /
                        std::pow(inp_model.TimeNorm(), 2.0)),
      force_norm(inp_model.MassNorm() * inp_model.LengthNorm() /
                 std::pow(inp_model.TimeNorm(), 2.0)),
      stress_norm(inp_model.MassNorm() / (std::pow(inp_model.TimeNorm(), 2.0) *
                                          inp_model.LengthNorm())),
      inertia_norm(inp_model.MassNorm() *
                   std::pow(inp_model.LengthNorm(), 2.0)),
      gravitational_constant_norm(
          std::pow(inp_model.LengthNorm(), 3.0) /
          (inp_model.MassNorm() * std::pow(inp_model.TimeNorm(), 2.0))),
      gravitational_constant(6.6743015 * std::pow(10.0, -11.0) *
                             inp_model.MassNorm() *
                             std::pow(inp_model.TimeNorm(), 2.0) /
                             std::pow(inp_model.LengthNorm(), 3.0)),
      potential_norm(
          std::pow(inp_model.LengthNorm() / inp_model.TimeNorm(), 2.0)) {

   _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoly + 1);

   // std::cout << "Hello\n";
   // if (lMax > 1)
   _grid = Grid(lMax, 2);

   // polynomial order
   _poly_ord = _q.N() - 1;
   // std::cout << "Hello 1.1\n";

   // find radial mesh
   // std::cout << "Hello 1.2\n";
   node_data = Radial_Tools::RadialMesh(inp_model, _q, max_radial_step, maxrad);
   _num_layers = node_data.NumberOfElements();
   // std::cout << "Hello 2\n";
   // default h
   auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
   _vec_h = vvvecdb(_num_layers, vvecdb(_q.N(), vecdb(spatialsize, 0.0)));
   {
      _vec_j = std::vector<std::vector<std::vector<double>>>(
          _num_layers, std::vector<std::vector<double>>(
                           _q.N(), std::vector<double>(spatialsize, 1.0)));
      Eigen::Matrix3cd mat_a0;
      mat_a0 << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
      _vec_a =
          vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize, mat_a0)));
      _vec_invF =
          vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize, mat_a0)));
   }
   // std::cout << "Hello 3\n";
   // finding _mat_gaussderiv
   _mat_gaussderiv.resize(_poly_ord + 1, _poly_ord + 1);
   {
      auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
                                                    _q.Points().end());
      for (int idxi = 0; idxi < _poly_ord + 1; ++idxi) {
         for (int idxj = 0; idxj < _poly_ord + 1; ++idxj) {
            _mat_gaussderiv(idxi, idxj) = pleg.Derivative(idxi, _q.X(idxj));
         }
      }
   }
   // std::cout << "Hello 4\n";
   // finding the spectral element grid
   _spectral_info = SpectralElementTools::MatrixWeakForm(node_data, _q);

   // lambda to convert between radius and depth (as in depth needed for
   // tomography model)
   auto rad_to_depth = [](double radius, double outerradius,
                          double lengthnorm) {
      double depth = outerradius - radius;
      depth *= lengthnorm / 1000.0;
      return depth;
   };
   auto colatitude_to_latitude = [](double colatitude) {
      double pi_db = 3.1415926535897932384626433;
      double latitude = pi_db / 2.0 - colatitude;
      latitude *= 180.0 / pi_db;
      return latitude;
   };
   auto longitude_to_degrees = [](double longitude) {
      double multfact = 180.0 / 3.1415926535897932384626433;
      return multfact * longitude;
   };

   // density
   // find all model information
   _vec_density.resize(_num_layers);

   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      // temporary vectors
      vvecdb tmp_density;   // to contain all values within the element
      int laynum = node_data.LayerNumber(idxelem);

      // loop through polynomial indices
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         double rad_current =
             node_data.NodeRadius(idxelem, idxnode);   // current radius
         double density_1D = 0.0;                      // initialise density

         // set density iff inside input model
         if (laynum < inp_model.NumberOfLayers()) {
            density_1D = inp_model.Density(laynum)(rad_current);
         }

         // fill out all values at particular radius
         vecdb tmp_vec_level_density(spatialsize, density_1D);
         auto depth = rad_to_depth(rad_current, node_data.PlanetRadius(),
                                   inp_model.LengthNorm());

         // check within tomography model
         if (depth > tomo_model.GetDepths()[0] &&
             depth < tomo_model.GetDepths().back()) {
            // std::cout << rad_current << "\n";
            // loop over all gridpoints
            int idxspatial = 0;
            for (auto idxt : _grid.CoLatitudes()) {
               auto latitude = colatitude_to_latitude(idxt);
               for (auto idxp : _grid.Longitudes()) {
                  auto longitude = longitude_to_degrees(idxp);
                  tmp_vec_level_density[idxspatial] *=
                      (1.0 + 0.005 * tomo_model.GetValueAt(depth, longitude,
                                                           latitude));
                  // tmp_vec_level_density[idxspatial] *= 1.0;
                  ++idxspatial;
               }
            }
         }

         tmp_density.push_back(tmp_vec_level_density);
      }

      // assign to the idxelem element
      _vec_density[idxelem] = tmp_density;
   }
};

// with path to file for 1D and tomography models
// Density3D::Density3D(const std::string &pathto1D, const int npoly,
//                      const int lMax, double max_radial_step, double
//                      maxrad)
//     : Density3D::Density3D(TestTools::EarthModel(pathto1D),
//                            PlanetaryModel::TomographyZeroModel(),
//                            npoly, lMax, max_radial_step, maxrad) {};

// with path to file for 1D and tomography models
// Density3D::Density3D(const std::string &pathto1D, const std::string
// &pathtotomo,
//                      const int npoly, const int lMax, double
//                      max_radial_step, double maxrad)
//     : Density3D::Density3D(TestTools::EarthModel(pathto1D),
//                            Tomography(pathtotomo), npoly, lMax,
//                            max_radial_step, maxrad) {};

}   // namespace GeneralEarthModels
#endif