#ifndef TEST_TOOLS_H
#define TEST_TOOLS_H

#include <Eigen/Core>
#include <Eigen/Dense>
// #include <Eigen/CholmodSupport>
#include <Eigen/Sparse>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>
// #include <TomographyModels/All>
#include "SphericalGeometryPreconditioner.h"
#include "Shared_Utilities.h"
#include "Model_Input.h"
#include <FFTWpp/Ranges>
#include <GSHTrans/All>
#include <TomographyModels/All>
#include <algorithm>
#include <complex>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>

#include <chrono>
#include <random>
#include <vector>

namespace TestTools {

enum class CoefficientOrdering {
   RadialClumped,
   SphericalDecomp
};

using MATRIX3cd = Eigen::Matrix<std::complex<double>, 3, 3>;
using EARTHMATRIX3 = std::vector<std::vector<MATRIX3cd>>;
using EARTHVEC = std::vector<std::vector<std::vector<double>>>;
using RADIUSVEC = std::vector<std::vector<std::complex<double>>>;
using EARTHVECX = std::vector<RADIUSVEC>;



}   // namespace TestTools

// if (idxinner > 0) {
//             if (radius == layered_radii[laynum][idxinner - 1]) {
//                layered_radii.push_back({radius});
//                layered_rho.push_back({rho});
//                layered_vpv.push_back({vpv});
//                layered_vsv.push_back({vsv});
//                layered_qkappa.push_back({qkappa});
//                layered_qshear.push_back({qshear});
//                layered_vph.push_back({vph});
//                layered_vsh.push_back({vsh});
//                layered_eta.push_back({eta});
//                idxinner = 0;
//                ++laynum;
//             } else {
//                layered_radii[laynum].push_back(radius);
//                layered_rho[laynum].push_back(rho);
//                layered_vpv[laynum].push_back(vpv);
//                layered_vsv[laynum].push_back(vsv);
//                layered_qkappa[laynum].push_back(qkappa);
//                layered_qshear[laynum].push_back(qshear);
//                layered_vph[laynum].push_back(vph);
//                layered_vsh[laynum].push_back(vsh);
//                layered_eta[laynum].push_back(eta);
//             }
//          } else {
//             layered_radii[laynum].push_back(radius);
//             layered_rho[laynum].push_back(rho);
//             layered_vpv[laynum].push_back(vpv);
//             layered_vsv[laynum].push_back(vsv);
//             layered_qkappa[laynum].push_back(qkappa);
//             layered_qshear[laynum].push_back(qshear);
//             layered_vph[laynum].push_back(vph);
//             layered_vsh[laynum].push_back(vsh);
//             layered_eta[laynum].push_back(eta);
//          }

#endif
