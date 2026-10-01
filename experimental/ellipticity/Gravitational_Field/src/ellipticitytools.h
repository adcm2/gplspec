#ifndef TESTELLIP_TOOLS_H
#define TESTELLIP_TOOLS_H

#include "SphericalGeometryPreconditioner.h"
#include "Spherical_Integrator.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <FFTWpp/Ranges>
#include <GSHTrans/All>
#include <GaussQuad/All>
#include <Gravitational_Field/Test>
#include <Interpolation/All>
#include <PlanetaryModel/All>
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

namespace EllipticityTools {

using namespace GravityFunctions;

using MATRIX3cd = Eigen::Matrix<std::complex<double>, 3, 3>;
using EARTHMATRIX3 = std::vector<std::vector<MATRIX3cd>>;
using EARTHVEC = std::vector<std::vector<std::vector<double>>>;
using RADIUSVEC = std::vector<std::vector<std::complex<double>>>;
using EARTHVECX = std::vector<RADIUSVEC>;

class FDEllipticitySolver {};

}   // namespace EllipticityTools

#endif