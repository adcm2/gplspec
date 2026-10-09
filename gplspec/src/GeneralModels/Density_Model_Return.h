#ifndef DENSITY_EARTH_MODEL_RETURN_3D_H
#define DENSITY_EARTH_MODEL_RETURN_3D_H

#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <vector>
#include "../Radial_Tools.h"
#include "../Shared_Utilities.h"
#include "Earth_Density_Models_3D.h"
#include <GSHTrans/All>

namespace GPLSpec::detail {

struct RotationOutputAngles {
   double alpha;
   double beta;
   double gamma;
};

inline RotationOutputAngles
RotationOutputAnglesFor(const std::vector<double> &vec_p1,
                        const std::vector<double> &vec_p2) {
   double theta1 = vec_p1[0];
   double phi1 = vec_p1[1];
   double theta2 = vec_p2[0];
   double phi2 = vec_p2[1];
   double cosphi = std::cos(theta2) * std::cos(theta1) +
                   std::sin(theta2) * std::sin(theta1) * std::cos(phi2 - phi1);
   double sinphi = std::sqrt(1 - cosphi * cosphi);
   auto tmp1 = std::sin(theta2) * std::cos(theta1) * std::cos(phi2) -
               std::cos(theta2) * std::sin(theta1) * std::cos(phi1);
   auto tmp2 = std::cos(theta2) * std::sin(theta1) * std::sin(phi1) -
               std::sin(theta2) * std::cos(theta1) * std::sin(phi2);
   auto tmp3 = std::sin(theta2) * std::sin(theta1) * std::sin(phi2 - phi1);
   auto tmp4 = -std::cos(theta1) * cosphi + std::cos(theta2);
   auto tmp5 = std::cos(theta1) * sinphi;

   if (std::abs(tmp1) < std::numeric_limits<double>::epsilon()) {
      tmp1 = 0.0;
   }
   if (std::abs(tmp2) < std::numeric_limits<double>::epsilon()) {
      tmp2 = 0.0;
   }
   if (std::abs(tmp4) < std::numeric_limits<double>::epsilon()) {
      tmp4 = 0.0;
   }
   if (std::abs(tmp5) < std::numeric_limits<double>::epsilon()) {
      tmp5 = 0.0;
   }

   double alpha = std::atan2(tmp1, tmp2);
   double beta = std::acos(tmp3 / sinphi);
   double gamma = std::atan2(tmp4, tmp5);
   return {alpha, beta, gamma};
}

inline std::vector<Eigen::MatrixXcd>
RotationOutputMatrices(int lMax, double alpha, double beta, double gamma) {
   auto vec_wig = std::vector<Eigen::MatrixXcd>(lMax + 1);
   auto wigtemp =
       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                        GSHTrans::Single, GSHTrans::ColumnMajor>(
           lMax, lMax, lMax, beta);
   for (int l = 0; l < lMax + 1; ++l) {
      Eigen::MatrixXcd mat_tmp = Eigen::MatrixXcd::Zero(2 * l + 1, 2 * l + 1);
      auto multval = GPLSpec::detail::RotationHarmonicNormalization(l);
      for (int m = -l; m < l + 1; ++m) {
         auto dl = wigtemp[m];
         for (int mp = -l; mp < l + 1; ++mp) {
            std::complex<double> i1(0.0, 1.0);
            auto tmpmult =
                multval * exp(i1 * (static_cast<double>(m) * gamma +
                                    static_cast<double>(mp) * alpha));
            mat_tmp(m + l, mp + l) = dl[l, mp] * tmpmult;
         }
      }
      vec_wig[l] = mat_tmp;
   }
   return vec_wig;
}

}   // namespace GPLSpec::detail

namespace GeneralEarthModels {

// info on nodes
auto
Density3D::Num_Elements() const {
   return _num_layers;
};
auto
Density3D::Poly_Order() const {
   return _poly_ord;
};
const Radial_Tools::RadialMesh &
Density3D::Node_Information() const {
   return node_data;
};

// interpolation polynomial derivative
auto
Density3D::GaussDerivative() const {
   return _mat_gaussderiv;
};
auto
Density3D::GaussDerivative(int idxi, int idxj) const {
   return _mat_gaussderiv(idxi, idxj);
};
auto
Density3D::SpectralElementInformation() {
   return _spectral_info;
};
const Density3D::Grid &
Density3D::GSH_Grid() const {
   return _grid;
};

// density
auto
Density3D::Density() const {
   return _vec_density;
};
//    auto DensityAtRadialNode(int, int) const;
//    auto Density_Point(int, int, int) const;

// h
auto
Density3D::Mapping() const {
   return _vec_h;
};
//    auto MappingAtRadialNode(int, int) const;
//    auto Mapping_Point(int, int, int) const;

auto
Density3D::LaplaceTensor() const {
   return _vec_a;
};

const std::vector<std::vector<std::vector<Eigen::Matrix3cd>>> &
Density3D::LaplaceTensorRef() const {
   return _vec_a;
};
//    auto LaplaceTensorAtRadialNode(int, int) const;
//    auto LaplaceTensor_Point(int, int, int) const;

// lower triangle
void
Density3D::SetLowerTriangle() {
   _spectral_info.set_left_lower();
};
//    void SetLowerTriangle() { _spectral_info.set_left_lower(); };

// norms
double const
Density3D::LengthNorm() const {
   return length_norm;
};
double const
Density3D::MassNorm() const {
   return mass_norm;
};
double const
Density3D::TimeNorm() const {
   return time_norm;
}
double const
Density3D::DensityNorm() const {
   return density_norm;
};
double const
Density3D::InertiaNorm() const {
   return inertia_norm;
};
double const
Density3D::VelocityNorm() const {
   return velocity_norm;
};
double const
Density3D::AccelerationNorm() const {
   return acceleration_norm;
};
double const
Density3D::ForceNorm() const {
   return force_norm;
};
double const
Density3D::StressNorm() const {
   return stress_norm;
};
double const
Density3D::GravitationalConstantNorm() const {
   return gravitational_constant_norm;
};
double const
Density3D::GravitationalConstant() const {
   return gravitational_constant;
};
double const
Density3D::PotentialNorm() const {
   return potential_norm;
};

// density
auto
Density3D::DensityAtRadialNode(int idxelem, int idxnode) const {
   return _vec_density[idxelem][idxnode];
};
auto
Density3D::Density_Point(int idxelem, int idxnode, int idxspatial) const {
   return _vec_density[idxelem][idxnode][idxspatial];
};

// mapping
auto
Density3D::MappingAtRadialNode(int idxelem, int idxnode) const {
   return _vec_h[idxelem][idxnode];
};
auto
Density3D::Mapping_Point(int idxelem, int idxnode, int idxspatial) const {
   return _vec_h[idxelem][idxnode][idxspatial];
};

auto
Density3D::Jacobian() const {
   return _vec_j;
};
auto
Density3D::JacobianAtRadialNode(int idxelem, int idxnode) const {
   return _vec_j[idxelem][idxnode];
};
auto
Density3D::Jacobian_Point(int idxelem, int idxnode, int idxspatial) const {
   return _vec_j[idxelem][idxnode][idxspatial];
};

// F^{-1}
auto
Density3D::InverseF() const {
   return _vec_invF;
};
const std::vector<std::vector<std::vector<Eigen::Matrix3cd>>> &
Density3D::InverseFRef() const {
   return _vec_invF;
};
auto
Density3D::InverseFAtRadialNode(int idxelem, int idxnode) const {
   return _vec_invF[idxelem][idxnode];
};
auto
Density3D::InverseF_Point(int idxelem, int idxnode, int idxspatial) const {
   return _vec_invF[idxelem][idxnode][idxspatial];
};

// return Laplace tensor a = J C^{-1}
auto
Density3D::LaplaceTensorAtRadialNode(int idxelem, int idxnode) const {
   return _vec_a[idxelem][idxnode];
};
auto
Density3D::LaplaceTensor_Point(int idxelem, int idxnode, int idxspatial) const {
   return _vec_a[idxelem][idxnode][idxspatial];
};

// return physical radius at a point
auto
Density3D::PhysicalRadius_Point(int idxelem, int idxnode,
                                int idxspatial) const {
   double refrad = node_data.NodeRadius(idxelem, idxnode);
   return refrad + _vec_h[idxelem][idxnode][idxspatial];
};
auto
Density3D::PhysicalRadius_Line(int idxspatial) const {
   std::vector<double> out_radii;
   for (int idx = 0; idx < this->Num_Elements(); ++idx) {
      double refrad = this->PhysicalRadius_Point(idx, 0, idxspatial);
      out_radii.push_back(refrad);
   }
   out_radii.push_back(node_data.OuterRadius());
   return out_radii;
};

auto
Density3D::Volume() const {

   // need to run through the model, take the spherical harmonic transform at
   // each node, then get sum
   int tot_elem = 0;
   for (int idxelem = 0; idxelem < this->Node_Information().NumberOfElements();
        ++idxelem) {
      if (Node_Information().LayerNumber(idxelem) ==
          Node_Information().NumberOfLayers() - 1) {
         tot_elem = idxelem;
         break;
      }
   }
   // size of vector to contain all coefficients in spherical harmonic transform
   auto coefficientnumberall = GSHTrans::GSHIndices<GSHTrans::All>(
                                   _grid.MaxDegree(), _grid.MaxDegree(), 0)
                                   .Size();
   auto coefficientnumberpos = GSHTrans::GSHIndices<GSHTrans::NonNegative>(
                                   _grid.MaxDegree(), _grid.MaxDegree(), 0)
                                   .Size();
   double volret = 0.0;
   for (int idxelem = 0; idxelem < tot_elem; ++idxelem) {
      double voltemp = 0.0;
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         std::vector<std::complex<double>> vec_jsph(coefficientnumberpos, 0.0);
         _grid.ForwardTransformation(_grid.MaxDegree(), 0,
                                     _vec_j[idxelem][idxnode], vec_jsph);
         voltemp += std::pow(node_data.NodeRadius(idxelem, idxnode), 2.0) *
                    vec_jsph[0].real() * _q.W(idxnode);
      }
      voltemp *= node_data.ElementWidth(idxelem) / 2.0;
      volret += voltemp;
   }
   return volret * sqrt(4.0 * 3.1415926535) * std::pow(this->LengthNorm(), 3.0);
};

auto
Density3D::Mass() const {

   // need to run through the model, take the spherical harmonic transform at
   // each node, then get sum
   int tot_elem = 0;
   for (int idxelem = 0; idxelem < this->Node_Information().NumberOfElements();
        ++idxelem) {
      if (Node_Information().LayerNumber(idxelem) ==
          Node_Information().NumberOfLayers() - 1) {
         tot_elem = idxelem;
         break;
      }
   }
   // size of vector to contain all coefficients in spherical harmonic transform
   int lMax = _grid.MaxDegree();
   auto coefficientnumberpos =
       GSHTrans::GSHIndices<GSHTrans::NonNegative>(lMax, lMax, 0).Size();
   double massret = 0.0;
   for (int idxelem = 0; idxelem < tot_elem; ++idxelem) {
      double masstemp = 0.0;
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         std::vector<std::complex<double>> vec_jsph(coefficientnumberpos, 0.0);
         _grid.ForwardTransformation(lMax, 0, _vec_density[idxelem][idxnode],
                                     vec_jsph);
         masstemp += std::pow(node_data.NodeRadius(idxelem, idxnode), 2.0) *
                     vec_jsph[0].real() * _q.W(idxnode);
      }
      masstemp *= node_data.ElementWidth(idxelem) / 2.0;
      massret += masstemp;
   }
   // GSHTrans uses the full-precision orthonormal Y00 normalization.
   return massret * std::sqrt(4.0 * std::numbers::pi_v<double>) *
          this->MassNorm();
};

auto
Density3D::SingleEigenVectorToGeneralFormat(
    const Eigen::VectorXcd &vec_pot) const {

   // e.g. for potential take in vector that comes out of solution and convert
   // to vector of vectors (of vectors) of lm coefficients:
   std::vector<std::vector<std::vector<std::complex<double>>>> vec_output(
       _num_layers,
       std::vector<std::vector<std::complex<double>>>(
           _poly_ord + 1, std::vector<std::complex<double>>(
                              std::pow(_grid.MaxDegree() + 1, 2), 0.0)));

   // number of l,m values
   std::size_t numlm = std::pow(_grid.MaxDegree() + 1, 2);

   // fill out
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
      std::size_t idxpot = idxelem * _poly_ord * numlm;
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         std::size_t idxlm = 0;
         for (int idxl = 0; idxl < _grid.MaxDegree() + 1; ++idxl) {
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               vec_output[idxelem][idxnode][idxlm] = vec_pot(idxpot);
               ++idxlm;
               ++idxpot;
            }
         }
      }
   }

   return vec_output;
};


}   // namespace GeneralEarthModels
#endif
