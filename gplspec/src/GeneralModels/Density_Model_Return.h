#ifndef DENSITY_EARTH_MODEL_RETURN_3D_H
#define DENSITY_EARTH_MODEL_RETURN_3D_H

#include <string>
#include <fstream>

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
auto
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
Density3D::q() const {
   return _q;
};
auto
Density3D::SpectralElementInformation() {
   return _spectral_info;
};
auto
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
                                   .size();
   auto coefficientnumberpos = GSHTrans::GSHIndices<GSHTrans::NonNegative>(
                                   _grid.MaxDegree(), _grid.MaxDegree(), 0)
                                   .size();
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
       GSHTrans::GSHIndices<GSHTrans::NonNegative>(lMax, lMax, 0).size();
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
   return massret * sqrt(4.0 * 3.1415926535) * this->MassNorm();
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

auto
Density3D::PowerSTD(
    const std::vector<std::vector<std::vector<std::complex<double>>>> &vec_pot)
    const {

   // e.g. for potential take in vector that comes out of solution and convert
   // to vector of vectors (of vectors) of lm coefficients:
   std::vector<std::vector<std::vector<double>>> vec_output(
       _num_layers,
       std::vector<std::vector<double>>(
           _poly_ord + 1, std::vector<double>(_grid.MaxDegree() + 1, 0.0)));
   // number of l,m values
   // std::size_t numlm = std::pow(_grid.MaxDegree() + 1, 2);

   // fill out
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         std::size_t idxlm = 0;
         double df = std::norm(vec_pot[idxelem][idxnode][0]);
         for (int idxl = 0; idxl < _grid.MaxDegree() + 1; ++idxl) {
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               std::complex<double> tmp = vec_pot[idxelem][idxnode][idxlm++];
               vec_output[idxelem][idxnode][idxl] += std::norm(tmp) / df;
            }
         }
      }
   }

   return vec_output;
};

// rotation to output a specific slice
auto
Density3D::RotateSliceToEquator(
    std::vector<double> &vec_p1, std::vector<double> &vec_p2,
    const std::vector<std::vector<std::vector<std::complex<double>>>>
        &vec_potential) {
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
   auto tmp4 = std::cos(theta1) * cosphi - std::cos(theta2);
   auto tmp5 = std::cos(theta1) * sinphi;

   if (tmp1 < std::numeric_limits<double>::epsilon()) {
      tmp1 = 0.0;
   }
   if (tmp2 < std::numeric_limits<double>::epsilon()) {
      tmp2 = 0.0;
   }
   if (tmp4 < std::numeric_limits<double>::epsilon()) {
      tmp4 = 0.0;
   }
   if (tmp5 < std::numeric_limits<double>::epsilon()) {
      tmp5 = 0.0;
   }

   double alpha = std::atan2(tmp1, tmp2);
   double beta = std::acos(tmp3 / sinphi);
   double gamma = std::atan2(tmp4, tmp5);

   // declare vector of matrices
   // int lmax_v = 2;
   // double theta_rot = std::numbers::pi_v<double> / 2.0;
   auto vec_wig = std::vector<Eigen::MatrixXcd>(_grid.MaxDegree() + 1);
   for (int l = 0; l < _grid.MaxDegree() + 1; ++l) {
      // temporary
      Eigen::MatrixXcd mat_tmp = Eigen::MatrixXcd::Zero(2 * l + 1, 2 * l + 1);
      int rowidx = 0;
      auto multval = std::sqrt((4.0 * 3.1415926535) / (2 * l + 1));

      // fill out matrix
      for (int m = -l; m < l + 1; ++m) {
         auto wigtemp = GSHTrans::Wigner(l, l, m, beta);
         auto dl = wigtemp(l);
         int colidx = 0;
         for (int mp = -l; mp < l + 1; ++mp) {
            std::complex<double> i1(0.0, 1.0);
            auto tmpmult =
                multval * exp(i1 * (static_cast<double>(m) * gamma +
                                    static_cast<double>(mp) * alpha));
            mat_tmp(rowidx, colidx) = dl(mp) * tmpmult;
            ++colidx;
         }
         ++rowidx;
      }
      vec_wig[l] = mat_tmp;
      // std::cout << "\n" << mat_tmp << "\n";
   }
   int nelem = this->Num_Elements();
   // double normfactor = this->PotentialNorm();
   std::vector<std::vector<std::vector<std::complex<double>>>> vec_output(
       _num_layers,
       std::vector<std::vector<std::complex<double>>>(
           _poly_ord + 1, std::vector<std::complex<double>>(
                              std::pow(_grid.MaxDegree() + 1, 2), 0.0)));
   for (int i = 0; i < nelem; ++i) {
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         int idxlm = 0;
         for (int idxl = 0; idxl < this->GSH_Grid().MaxDegree() + 1; ++idxl) {
            int rowidx = 0;
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               int idxlmp = idxl * idxl;
               int colidx = 0;
               for (int idxmp = -idxl; idxmp < idxl + 1; ++idxmp) {
                  vec_output[i][idxnode][idxlm] +=
                      vec_wig[idxl](rowidx, colidx) *
                      vec_potential[i][idxnode][idxlmp];
                  ++idxlmp;
                  ++colidx;
               }
               ++rowidx;
               ++idxlm;
            }
         }
      }
   };
   return vec_output;
};

void
Density3D::ReferentialOutputAtElement(
    const std::string str_pathtofolder,
    const std::vector<std::vector<std::vector<std::complex<double>>>>
        &vec_fullinformation,
    bool sensitivitykernel = false) const {

   assert(vec_fullinformation.size() == this->Num_Elements());

   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   int nelem = this->Num_Elements();
   double normfactor = this->PotentialNorm();
   std::string pathtofile;
   if (sensitivitykernel) {
      // std::cout << "\n\n\n\n\n HELLO THERE\n\n\n\n\n";
      normfactor *= 1.0 / this->MassNorm();
      pathtofile = str_pathtofolder + "/SensitivityKernel.out";
   } else {
      pathtofile = str_pathtofolder + "/MatrixSolution.out";
   }

   auto file2 = std::ofstream(pathtofile);

   for (int i = 0; i < nelem; ++i) {
      double currentrad;
      if (i == nelem) {
         currentrad = this->Node_Information().ElementUpperRadius(i) *
                      this->LengthNorm();
      } else {
         currentrad = this->Node_Information().ElementLowerRadius(i) *
                      this->LengthNorm();
      }

      //
      // first output
      //
      file2 << std::setprecision(16) << currentrad;

      //
      // working through lm
      //
      int idxlm = 0;
      for (int idxl = 0; idxl < this->GSH_Grid().MaxDegree() + 1; ++idxl) {
         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
            file2 << ";" << vec_fullinformation[i][0][idxlm].real() * normfactor
                  << ";"
                  << vec_fullinformation[i][0][idxlm].imag() * normfactor;
            ++idxlm;
         }
      }
      file2 << std::endl;
   };

   {
      double currentrad =
          this->Node_Information().ElementUpperRadius(nelem - 1) *
          this->LengthNorm();
      file2 << std::setprecision(16) << currentrad;
      int idxlm = 0;
      for (int idxl = 0; idxl < this->GSH_Grid().MaxDegree() + 1; ++idxl) {
         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
            file2 << ";"
                  << vec_fullinformation[nelem - 1][this->Poly_Order()][idxlm]
                             .real() *
                         normfactor
                  << ";"
                  << vec_fullinformation[nelem - 1][this->Poly_Order()][idxlm]
                             .imag() *
                         normfactor;
            ++idxlm;
         }
      }
      file2 << std::endl;
   }

   // close
   file2.close();
};

// physical radius
void
Density3D::ReferentialOutputSlice(
    const std::string str_pathtofolder,
    const std::vector<std::vector<std::vector<std::complex<double>>>>
        &vec_fullinformation,
    bool sensitivitykernel = false) const {

   assert(vec_fullinformation.size() == this->Num_Elements());

   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   // std::string pathtofile = str_pathtofolder + "/MatrixSolution.out";

   // int nelem = this->Num_Elements();
   int nelem = this->node_data.NumberOfElements();
   double normfactor = this->PotentialNorm();
   std::string pathtofile;
   if (sensitivitykernel) {
      // std::cout << "\n\n\n\n\n HELLO THERE\n\n\n\n\n";
      normfactor *= 1.0 / this->MassNorm();
      pathtofile = str_pathtofolder + "/SensitivityKernelSlice.out";
   } else {
      pathtofile = str_pathtofolder + "/MatrixSolutionSlice.out";
   }
   auto file2 = std::ofstream(pathtofile);

   int ncolat = _grid.NumberOfCoLatitudes();
   bool evenn = false;
   if (ncolat % 2 == 0) {
      evenn = true;
   }
   int idxlow = 0;
   if (evenn) {
      idxlow = ncolat / 2 - 1;
   } else {
      idxlow = std::floor(static_cast<double>(ncolat) / 2.0);
   }
   for (int i = 0; i < nelem + 1; ++i) {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize);
      int nlayer = i;
      int nnode = 0;
      if (i == nelem) {
         nlayer = nelem - 1;
         nnode = this->Poly_Order();
      }
      _grid.InverseTransformation(lmax, 0, vec_fullinformation[nlayer][nnode],
                                  vectmp);

      // output
      int idxspatial = 0;
      int idxcolat = 0;
      std::vector<double> vec_reals(_grid.NumberOfLongitudes(), 0.0),
          vec_imags(_grid.NumberOfLongitudes(), 0.0);
      for (auto it : _grid.CoLatitudes()) {
         int idxlong = 0;
         for (auto ip : _grid.Longitudes()) {
            // file2 << std::setprecision(16)
            //       << (node_data.ElementLowerRadius(i)) * this->LengthNorm()
            //       << ";" << it << ";" << ip << ";"
            //       << vectmp[idxspatial].real() * normfactor << ";"
            //       << vectmp[idxspatial].imag() * normfactor;
            if (evenn) {
               if (idxcolat == idxlow || idxcolat == idxlow + 1) {
                  vec_reals[idxlong] +=
                      0.5 * vectmp[idxspatial].real() * normfactor;
                  vec_imags[idxlong] +=
                      0.5 * vectmp[idxspatial].imag() * normfactor;
               }
            } else {
               if (idxcolat == idxlow) {
                  vec_reals[idxlong] += vectmp[idxspatial].real() * normfactor;
                  vec_imags[idxlong] += vectmp[idxspatial].imag() * normfactor;
               }
            }
            // if (idxspatial < spatialsize - 1) {
            //    file2 << ";";
            // }
            ++idxspatial;
            ++idxlong;
         }
         ++idxcolat;
      }
      idxcolat = 0;
      {
         int idxlong = 0;
         double currentrad = 0.0;
         if (i < nelem) {
            currentrad = node_data.ElementLowerRadius(nlayer);

         } else {
            currentrad = node_data.ElementUpperRadius(nlayer);
         }
         for (auto ip : _grid.Longitudes()) {
            file2 << std::setprecision(16) << currentrad * this->LengthNorm()
                  << ";" << ip << ";" << vec_reals[idxlong] << ";"
                  << vec_imags[idxlong];
            ++idxspatial;
            ++idxlong;
            if (idxlong < _grid.NumberOfLongitudes()) {
               file2 << ";";
            }
         }

         file2 << std::endl;
      };
   };

   // close
   file2.close();
};

// physical radius
void
Density3D::PhysicalOutputSlice(
    const std::string str_pathtofolder,
    const std::vector<std::vector<std::vector<std::complex<double>>>>
        &vec_fullinformation,
    bool sensitivitykernel = false) const {

   assert(vec_fullinformation.size() == this->Num_Elements());

   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   // std::string pathtofile = str_pathtofolder + "/MatrixSolution.out";

   // int nelem = this->Num_Elements();
   int nelem = this->node_data.NumberOfElements();
   double normfactor = this->PotentialNorm();
   std::string pathtofile;
   if (sensitivitykernel) {
      // std::cout << "\n\n\n\n\n HELLO THERE\n\n\n\n\n";
      normfactor *= 1.0 / this->MassNorm();
      pathtofile = str_pathtofolder + "/SensitivityKernelSlice.out";
   } else {
      pathtofile = str_pathtofolder + "/MatrixSolutionSlice.out";
   }
   auto file2 = std::ofstream(pathtofile);

   int ncolat = _grid.NumberOfCoLatitudes();
   bool evenn = false;
   if (ncolat % 2 == 0) {
      evenn = true;
   }
   int idxlow = 0;
   if (evenn) {
      idxlow = ncolat / 2 - 1;
   } else {
      idxlow = std::floor(static_cast<double>(ncolat) / 2.0);
   }
   for (int i = 0; i < nelem + 1; ++i) {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize);
      int nlayer = i;
      int nnode = 0;
      if (i == nelem) {
         nlayer = nelem - 1;
         nnode = this->Poly_Order();
      }
      _grid.InverseTransformation(lmax, 0, vec_fullinformation[nlayer][nnode],
                                  vectmp);
      double rad_ref = 0.0;
      if (i < nelem) {
         rad_ref = node_data.ElementLowerRadius(nlayer);

      } else {
         rad_ref = node_data.ElementUpperRadius(nlayer);
      }

      // output
      int idxspatial = 0;
      int idxcolat = 0;
      std::vector<double> vec_reals(_grid.NumberOfLongitudes(), 0.0),
          vec_imags(_grid.NumberOfLongitudes(), 0.0),
          vec_radii(_grid.NumberOfLongitudes(), 0.0);
      for (auto it : _grid.CoLatitudes()) {
         int idxlong = 0;
         for (auto ip : _grid.Longitudes()) {
            // file2 << std::setprecision(16)
            //       << (node_data.ElementLowerRadius(i)) * this->LengthNorm()
            //       << ";" << it << ";" << ip << ";"
            //       << vectmp[idxspatial].real() * normfactor << ";"
            //       << vectmp[idxspatial].imag() * normfactor;

            if (evenn) {
               if (idxcolat == idxlow || idxcolat == idxlow + 1) {
                  vec_reals[idxlong] +=
                      0.5 * vectmp[idxspatial].real() * normfactor;
                  vec_imags[idxlong] +=
                      0.5 * vectmp[idxspatial].imag() * normfactor;

                  vec_radii[idxlong] +=
                      0.5 * (rad_ref + _vec_h[nlayer][nnode][idxspatial]) *
                      this->LengthNorm();
               }
            } else {
               if (idxcolat == idxlow) {
                  vec_reals[idxlong] = vectmp[idxspatial].real() * normfactor;
                  vec_imags[idxlong] = vectmp[idxspatial].imag() * normfactor;
                  vec_radii[idxlong] =
                      (rad_ref + _vec_h[nlayer][nnode][idxspatial]) *
                      this->LengthNorm();
               }
            }
            // if (idxspatial < spatialsize - 1) {
            //    file2 << ";";
            // }
            ++idxspatial;
            ++idxlong;
         }
         ++idxcolat;
      }
      idxcolat = 0;
      {
         int idxlong = 0;

         for (auto ip : _grid.Longitudes()) {

            // double physrad = rad_ref *;
            file2 << std::setprecision(16) << vec_radii[idxlong] << ";" << ip
                  << ";" << vec_reals[idxlong] << ";" << vec_imags[idxlong];
            ++idxspatial;
            ++idxlong;
            if (idxlong < _grid.NumberOfLongitudes()) {
               file2 << ";";
            }
         }

         file2 << std::endl;
      };
   };

   // close
   file2.close();
};

// physical radius
void
Density3D::PhysicalOutputAtElement(
    const std::string str_pathtofolder,
    const std::vector<std::vector<std::vector<std::complex<double>>>>
        &vec_fullinformation,
    bool sensitivitykernel = false) const {

   assert(vec_fullinformation.size() == this->Num_Elements());

   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   // std::string pathtofile = str_pathtofolder + "/MatrixSolution.out";

   // int nelem = this->Num_Elements();
   int nelem = this->node_data.NumberOfElements();
   double normfactor = this->PotentialNorm();
   std::string pathtofile;
   if (sensitivitykernel) {
      // std::cout << "\n\n\n\n\n HELLO THERE\n\n\n\n\n";
      normfactor *= 1.0 / this->MassNorm();
      pathtofile = str_pathtofolder + "/SensitivityKernelSpatial.out";
   } else {
      pathtofile = str_pathtofolder + "/MatrixSolution.out";
   }
   // std::cout << pathtofile << "\n\n\n";
   auto file2 = std::ofstream(pathtofile);

   for (int i = 0; i < nelem; ++i) {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize);
      _grid.InverseTransformation(lmax, 0, vec_fullinformation[i][0], vectmp);
      // if (i == 10) {
      //    for (int idxt = 0; idxt < 30; ++idxt) {
      //       std::cout << vec_fullinformation[i][0][idxt] << " "
      //                 << vectmp[idxt] * normfactor << "\n";
      //    }
      // }

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            file2 << std::setprecision(16)
                  << (node_data.ElementLowerRadius(i) +
                      _vec_h[i][0][idxspatial]) *
                         this->LengthNorm()
                  << ";" << it << ";" << ip << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   };

   {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize);
      _grid.InverseTransformation(
          lmax, 0, vec_fullinformation[nelem - 1][this->Poly_Order()], vectmp);

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            file2 << std::setprecision(16)
                  << (node_data.ElementUpperRadius(nelem - 1) +
                      _vec_h[nelem - 1][this->Poly_Order()][idxspatial]) *
                         this->LengthNorm()
                  << ";" << it << ";" << ip << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   }

   // close
   file2.close();
};

// physical radius
void
Density3D::ReferentialOutputRotated(
    const std::string str_pathtofolder, std::vector<double> &vec_p1,
    std::vector<double> &vec_p2,
    const std::vector<std::vector<std::vector<std::complex<double>>>>
        &vec_fullinformation) const {

   // get rotation information
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

   std::cout << "\ntmp1: " << tmp1 << ", tmp2: " << tmp2 << ", alpha: " << alpha
             << "\n";

   std::cout << "\n\n\n" << alpha << " " << beta << " " << gamma << "\n\n\n";
   // declare vector of matrices
   // int lmax_v = 2;
   // double theta_rot = std::numbers::pi_v<double> / 2.0;
   auto vec_wig = std::vector<Eigen::MatrixXcd>(_grid.MaxDegree() + 1);
   for (int l = 0; l < _grid.MaxDegree() + 1; ++l) {
      // temporary
      Eigen::MatrixXcd mat_tmp = Eigen::MatrixXcd::Zero(2 * l + 1, 2 * l + 1);
      int rowidx = 0;
      auto multval = std::sqrt((4.0 * 3.1415926535) / (2 * l + 1));

      // fill out matrix
      for (int m = -l; m < l + 1; ++m) {
         auto wigtemp = GSHTrans::Wigner(l, l, m, beta);
         auto dl = wigtemp(l);
         int colidx = 0;
         for (int mp = -l; mp < l + 1; ++mp) {
            std::complex<double> i1(0.0, 1.0);
            auto tmpmult =
                multval * exp(i1 * (static_cast<double>(m) * gamma +
                                    static_cast<double>(mp) * alpha));
            mat_tmp(rowidx, colidx) = dl(mp) * tmpmult;
            ++colidx;
         }
         ++rowidx;
      }
      vec_wig[l] = mat_tmp;
      // std::cout << "\n" << mat_tmp << "\n";
   }
   int nelem = this->Num_Elements();

   assert(vec_fullinformation.size() == this->Num_Elements());

   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   // std::string pathtofile = str_pathtofolder + "/MatrixSolutionRotated.out";
   std::string pathtofile = str_pathtofolder;
   auto file2 = std::ofstream(pathtofile);
   // int nelem = this->Num_Elements();
   // int nelem = this->node_data.NumberOfElements();
   double normfactor = this->PotentialNorm();
   // get vec_hlm
   auto coefficientnumber = GSHTrans::GSHIndices<GSHTrans::NonNegative>(
                                _grid.MaxDegree(), _grid.MaxDegree(), 0)
                                .size();
   auto coefficientnumberall = GSHTrans::GSHIndices<GSHTrans::All>(
                                   _grid.MaxDegree(), _grid.MaxDegree(), 0)
                                   .size();
   auto lMax = _grid.MaxDegree();
   // int npoly = _q.N() - 1;
   auto vec_hlm = std::vector<std::vector<std::vector<std::complex<double>>>>(
       _num_layers, std::vector<std::vector<std::complex<double>>>(
                        _poly_ord + 1, std::vector<std::complex<double>>(
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
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
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

   // perform rotations
   //  double normfactor = this->PotentialNorm();
   std::vector<std::vector<std::vector<std::complex<double>>>> vec_output(
       _num_layers,
       std::vector<std::vector<std::complex<double>>>(
           _poly_ord + 1, std::vector<std::complex<double>>(
                              std::pow(_grid.MaxDegree() + 1, 2), 0.0)));
   std::vector<std::vector<std::vector<std::complex<double>>>> vec_roth(
       _num_layers,
       std::vector<std::vector<std::complex<double>>>(
           _poly_ord + 1, std::vector<std::complex<double>>(
                              std::pow(_grid.MaxDegree() + 1, 2), 0.0)));
   for (int i = 0; i < nelem; ++i) {
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         int idxlm = 0;
         for (int idxl = 0; idxl < this->GSH_Grid().MaxDegree() + 1; ++idxl) {
            int rowidx = 0;
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               int idxlmp = idxl * idxl;
               int colidx = 0;
               for (int idxmp = -idxl; idxmp < idxl + 1; ++idxmp) {
                  vec_output[i][idxnode][idxlm] +=
                      vec_wig[idxl](rowidx, colidx) *
                      vec_fullinformation[i][idxnode][idxlmp];
                  vec_roth[i][idxnode][idxlm] += vec_wig[idxl](rowidx, colidx) *
                                                 vec_hlm[i][idxnode][idxlmp];
                  ++idxlmp;
                  ++colidx;
               }
               ++rowidx;
               ++idxlm;
            }
         }
      }
   };

   for (int i = 0; i < nelem; ++i) {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize),
          vechrot(spatialsize);
      _grid.InverseTransformation(lmax, 0, vec_output[i][0], vectmp);
      _grid.InverseTransformation(lmax, 0, vec_roth[i][0], vechrot);

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            file2 << std::setprecision(16)
                  << (node_data.ElementLowerRadius(i)) * this->LengthNorm()
                  << ";" << it << ";" << ip << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   };

   {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize),
          vechrot(spatialsize);
      _grid.InverseTransformation(
          lmax, 0, vec_output[nelem - 1][this->Poly_Order()], vectmp);
      _grid.InverseTransformation(
          lmax, 0, vec_roth[nelem - 1][this->Poly_Order()], vechrot);

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            file2 << std::setprecision(16)
                  << (node_data.ElementUpperRadius(nelem - 1)) *
                         this->LengthNorm()
                  << ";" << it << ";" << ip << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   }

   // close
   file2.close();
};

// physical radius
void
Density3D::ModelDensityOutputRotated(const std::string str_pathtofolder,
                                     std::vector<double> &vec_p1,
                                     std::vector<double> &vec_p2,
                                     bool physical = false) const {

   // get rotation information
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

   std::cout << "\ntmp1: " << tmp1 << ", tmp2: " << tmp2 << ", alpha: " << alpha
             << "\n";

   std::cout << "\n\n\n" << alpha << " " << beta << " " << gamma << "\n\n\n";
   // declare vector of matrices
   // int lmax_v = 2;
   // double theta_rot = std::numbers::pi_v<double> / 2.0;
   auto vec_wig = std::vector<Eigen::MatrixXcd>(_grid.MaxDegree() + 1);
   for (int l = 0; l < _grid.MaxDegree() + 1; ++l) {
      // temporary
      Eigen::MatrixXcd mat_tmp = Eigen::MatrixXcd::Zero(2 * l + 1, 2 * l + 1);
      int rowidx = 0;
      auto multval = std::sqrt((4.0 * 3.1415926535) / (2 * l + 1));

      // fill out matrix
      for (int m = -l; m < l + 1; ++m) {
         auto wigtemp = GSHTrans::Wigner(l, l, m, beta);
         auto dl = wigtemp(l);
         int colidx = 0;
         for (int mp = -l; mp < l + 1; ++mp) {
            std::complex<double> i1(0.0, 1.0);
            auto tmpmult =
                multval * exp(i1 * (static_cast<double>(m) * gamma +
                                    static_cast<double>(mp) * alpha));
            mat_tmp(rowidx, colidx) = dl(mp) * tmpmult;
            ++colidx;
         }
         ++rowidx;
      }
      vec_wig[l] = mat_tmp;
      // std::cout << "\n" << mat_tmp << "\n";
   }
   int nelem = this->Num_Elements();

   // assert(vec_fullinformation.size() == this->Num_Elements());

   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   // std::string pathtofile = str_pathtofolder + "/MatrixSolutionRotated.out";
   std::string pathtofile = str_pathtofolder;
   auto file2 = std::ofstream(pathtofile);
   // int nelem = this->Num_Elements();
   // int nelem = this->node_data.NumberOfElements();
   double normfactor = this->PotentialNorm();
   // get vec_hlm
   auto coefficientnumber = GSHTrans::GSHIndices<GSHTrans::NonNegative>(
                                _grid.MaxDegree(), _grid.MaxDegree(), 0)
                                .size();
   auto coefficientnumberall = GSHTrans::GSHIndices<GSHTrans::All>(
                                   _grid.MaxDegree(), _grid.MaxDegree(), 0)
                                   .size();
   auto lMax = _grid.MaxDegree();
   // int npoly = _q.N() - 1;
   auto vec_hlm = std::vector<std::vector<std::vector<std::complex<double>>>>(
       _num_layers, std::vector<std::vector<std::complex<double>>>(
                        _poly_ord + 1, std::vector<std::complex<double>>(
                                           coefficientnumberall, 0.0)));
   auto vec_rholm = vec_hlm;
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
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         std::vector<std::complex<double>> tmp_h(coefficientnumber),
             tmp_rho(coefficientnumber);
         auto spatialsize =
             _grid.Longitudes().size() * _grid.CoLatitudes().size();
         std::vector<double> vec_denstmp(spatialsize, 0.0);
         if (physical) {
            for (int idx = 0; idx < spatialsize; ++idx) {
               vec_denstmp[idx] = _vec_density[idxelem][idxnode][idx] /
                                  _vec_j[idxelem][idxnode][idx];
            }
         } else {
            // for (int idx = 0; idx < spatialsize; ++idx) {
            vec_denstmp = _vec_density[idxelem][idxnode];
            // }
         }

         _grid.ForwardTransformation(lMax, 0, _vec_h[idxelem][idxnode], tmp_h);
         _grid.ForwardTransformation(lMax, 0, vec_denstmp, tmp_rho);

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
               vec_rholm[idxelem][idxnode][idx] =
                   std::pow(-1.0, idxm) * std::conj(tmp_rho[idxlm]);
               ++idx;
            }

            //+ve m:
            for (int idxm = 0; idxm < idxl + 1; ++idxm) {
               int idxlm = indexreal(idxl, idxm);
               vec_hlm[idxelem][idxnode][idx] = tmp_h[idxlm];
               vec_rholm[idxelem][idxnode][idx] = tmp_rho[idxlm];
               ++idx;
            }
         }
      }
   }

   // need density in terms of spherical harmonics

   // perform rotations
   //  double normfactor = this->PotentialNorm();
   std::vector<std::vector<std::vector<std::complex<double>>>> vec_output(
       _num_layers,
       std::vector<std::vector<std::complex<double>>>(
           _poly_ord + 1, std::vector<std::complex<double>>(
                              std::pow(_grid.MaxDegree() + 1, 2), 0.0)));
   std::vector<std::vector<std::vector<std::complex<double>>>> vec_roth(
       _num_layers,
       std::vector<std::vector<std::complex<double>>>(
           _poly_ord + 1, std::vector<std::complex<double>>(
                              std::pow(_grid.MaxDegree() + 1, 2), 0.0)));
   for (int i = 0; i < nelem; ++i) {
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         int idxlm = 0;
         for (int idxl = 0; idxl < this->GSH_Grid().MaxDegree() + 1; ++idxl) {
            int rowidx = 0;
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               int idxlmp = idxl * idxl;
               int colidx = 0;
               for (int idxmp = -idxl; idxmp < idxl + 1; ++idxmp) {
                  vec_output[i][idxnode][idxlm] +=
                      vec_wig[idxl](rowidx, colidx) *
                      vec_rholm[i][idxnode][idxlmp];
                  vec_roth[i][idxnode][idxlm] += vec_wig[idxl](rowidx, colidx) *
                                                 vec_hlm[i][idxnode][idxlmp];
                  ++idxlmp;
                  ++colidx;
               }
               ++rowidx;
               ++idxlm;
            }
         }
      }
   };

   for (int i = 0; i < nelem; ++i) {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize),
          vechrot(spatialsize);
      _grid.InverseTransformation(lmax, 0, vec_output[i][0], vectmp);
      _grid.InverseTransformation(lmax, 0, vec_roth[i][0], vechrot);

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            double radout = 0.0;
            if (physical) {
               radout = (node_data.ElementLowerRadius(i) +
                         vechrot[idxspatial].real()) *
                        this->LengthNorm();
            } else {
               radout = (node_data.ElementLowerRadius(i)) * this->LengthNorm();
            }
            file2 << std::setprecision(16) << radout << ";" << it << ";" << ip
                  << ";" << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   };

   {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize),
          vechrot(spatialsize);
      _grid.InverseTransformation(
          lmax, 0, vec_output[nelem - 1][this->Poly_Order()], vectmp);
      _grid.InverseTransformation(
          lmax, 0, vec_roth[nelem - 1][this->Poly_Order()], vechrot);

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            double radout = 0.0;
            if (physical) {
               radout = (node_data.ElementUpperRadius(nelem - 1) +
                         vechrot[idxspatial].real()) *
                        this->LengthNorm();
            } else {
               radout = (node_data.ElementUpperRadius(nelem - 1)) *
                        this->LengthNorm();
            }
            file2 << std::setprecision(16) << radout << ";" << it << ";" << ip
                  << ";" << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   }

   // close
   file2.close();
};

// physical radius
void
Density3D::PhysicalOutputRotated(
    const std::string str_pathtofolder, std::vector<double> &vec_p1,
    std::vector<double> &vec_p2,
    const std::vector<std::vector<std::vector<std::complex<double>>>>
        &vec_fullinformation) const {

   // get rotation information
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

   std::cout << "\ntmp1: " << tmp1 << ", tmp2: " << tmp2 << ", alpha: " << alpha
             << "\n";

   std::cout << "\n\n\n" << alpha << " " << beta << " " << gamma << "\n\n\n";
   // declare vector of matrices
   // int lmax_v = 2;
   // double theta_rot = std::numbers::pi_v<double> / 2.0;
   auto vec_wig = std::vector<Eigen::MatrixXcd>(_grid.MaxDegree() + 1);
   for (int l = 0; l < _grid.MaxDegree() + 1; ++l) {
      // temporary
      Eigen::MatrixXcd mat_tmp = Eigen::MatrixXcd::Zero(2 * l + 1, 2 * l + 1);
      int rowidx = 0;
      auto multval = std::sqrt((4.0 * 3.1415926535) / (2 * l + 1));

      // fill out matrix
      for (int m = -l; m < l + 1; ++m) {
         auto wigtemp = GSHTrans::Wigner(l, l, m, beta);
         auto dl = wigtemp(l);
         int colidx = 0;
         for (int mp = -l; mp < l + 1; ++mp) {
            std::complex<double> i1(0.0, 1.0);
            auto tmpmult =
                multval * exp(i1 * (static_cast<double>(m) * gamma +
                                    static_cast<double>(mp) * alpha));
            mat_tmp(rowidx, colidx) = dl(mp) * tmpmult;
            ++colidx;
         }
         ++rowidx;
      }
      vec_wig[l] = mat_tmp;
      // std::cout << "\n" << mat_tmp << "\n";
   }
   int nelem = this->Num_Elements();

   assert(vec_fullinformation.size() == this->Num_Elements());

   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   // std::string pathtofile = str_pathtofolder + "/MatrixSolutionRotated.out";
   std::string pathtofile = str_pathtofolder;
   auto file2 = std::ofstream(pathtofile);
   // int nelem = this->Num_Elements();
   // int nelem = this->node_data.NumberOfElements();
   double normfactor = this->PotentialNorm();
   // get vec_hlm
   auto coefficientnumber = GSHTrans::GSHIndices<GSHTrans::NonNegative>(
                                _grid.MaxDegree(), _grid.MaxDegree(), 0)
                                .size();
   auto coefficientnumberall = GSHTrans::GSHIndices<GSHTrans::All>(
                                   _grid.MaxDegree(), _grid.MaxDegree(), 0)
                                   .size();
   auto lMax = _grid.MaxDegree();
   // int npoly = _q.N() - 1;
   auto vec_hlm = std::vector<std::vector<std::vector<std::complex<double>>>>(
       _num_layers, std::vector<std::vector<std::complex<double>>>(
                        _poly_ord + 1, std::vector<std::complex<double>>(
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
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
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

   // perform rotations
   //  double normfactor = this->PotentialNorm();
   std::vector<std::vector<std::vector<std::complex<double>>>> vec_output(
       _num_layers,
       std::vector<std::vector<std::complex<double>>>(
           _poly_ord + 1, std::vector<std::complex<double>>(
                              std::pow(_grid.MaxDegree() + 1, 2), 0.0)));
   std::vector<std::vector<std::vector<std::complex<double>>>> vec_roth(
       _num_layers,
       std::vector<std::vector<std::complex<double>>>(
           _poly_ord + 1, std::vector<std::complex<double>>(
                              std::pow(_grid.MaxDegree() + 1, 2), 0.0)));
   for (int i = 0; i < nelem; ++i) {
      for (int idxnode = 0; idxnode < _poly_ord + 1; ++idxnode) {
         int idxlm = 0;
         for (int idxl = 0; idxl < this->GSH_Grid().MaxDegree() + 1; ++idxl) {
            int rowidx = 0;
            for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               int idxlmp = idxl * idxl;
               int colidx = 0;
               for (int idxmp = -idxl; idxmp < idxl + 1; ++idxmp) {
                  vec_output[i][idxnode][idxlm] +=
                      vec_wig[idxl](rowidx, colidx) *
                      vec_fullinformation[i][idxnode][idxlmp];
                  vec_roth[i][idxnode][idxlm] += vec_wig[idxl](rowidx, colidx) *
                                                 vec_hlm[i][idxnode][idxlmp];
                  ++idxlmp;
                  ++colidx;
               }
               ++rowidx;
               ++idxlm;
            }
         }
      }
   };

   for (int i = 0; i < nelem; ++i) {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize),
          vechrot(spatialsize);
      _grid.InverseTransformation(lmax, 0, vec_output[i][0], vectmp);
      _grid.InverseTransformation(lmax, 0, vec_roth[i][0], vechrot);

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            file2 << std::setprecision(16)
                  << (node_data.ElementLowerRadius(i) +
                      vechrot[idxspatial].real()) *
                         this->LengthNorm()
                  << ";" << it << ";" << ip << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   };

   {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize),
          vechrot(spatialsize);
      _grid.InverseTransformation(
          lmax, 0, vec_output[nelem - 1][this->Poly_Order()], vectmp);
      _grid.InverseTransformation(
          lmax, 0, vec_roth[nelem - 1][this->Poly_Order()], vechrot);

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            file2 << std::setprecision(16)
                  << (node_data.ElementUpperRadius(nelem - 1) +
                      vechrot[idxspatial].real()) *
                         this->LengthNorm()
                  << ";" << it << ";" << ip << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   }

   // close
   file2.close();
};

void
Density3D::CartesianOutputAtElement(
    const std::string str_pathtofolder,
    const std::vector<std::vector<std::vector<std::complex<double>>>>
        &vec_fullinformation) const {

   assert(vec_fullinformation.size() == this->Num_Elements());

   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   std::string pathtofile = str_pathtofolder + "/MatrixSolution.out";
   auto file2 = std::ofstream(pathtofile);
   // int nelem = this->Num_Elements();
   int nelem = this->node_data.NumberOfElements();
   double normfactor = this->PotentialNorm();

   for (int i = 0; i < nelem; ++i) {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize);
      _grid.InverseTransformation(lmax, 0, vec_fullinformation[i][0], vectmp);

      // output
      int idxspatial = 0;
      file2 << this->Node_Information().LayerNumber(i) << ";";

      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            double physrad =
                (node_data.ElementLowerRadius(i) + _vec_h[i][0][idxspatial]) *
                this->LengthNorm();
            file2 << std::setprecision(16)
                  << physrad * std::sin(it) * std::cos(ip) << ";"
                  << physrad * std::sin(it) * std::sin(ip) << ";"
                  << physrad * std::cos(it) << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   };

   {

      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize);
      _grid.InverseTransformation(
          lmax, 0, vec_fullinformation[nelem - 1][this->Poly_Order()], vectmp);

      // output
      int idxspatial = 0;
      file2 << this->Node_Information().LayerNumber(nelem - 1) << ";";
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            double physrad =
                (node_data.ElementUpperRadius(nelem - 1) +
                 _vec_h[nelem - 1][this->Poly_Order()][idxspatial]) *
                this->LengthNorm();
            file2 << std::setprecision(16)
                  << physrad * std::sin(it) * std::cos(ip) << ";"
                  << physrad * std::sin(it) * std::sin(ip) << ";"
                  << physrad * std::cos(it) << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   }

   // close
   file2.close();
};
// // physical radius
// void
// Density3D::SphericalHarmonicOutputAtPhysicalRadius(
//     const std::string str_pathtofolder,
//     const std::vector<std::vector<std::vector<std::complex<double>>>>
//         &vec_fullinformation) const {

//    assert(vec_fullinformation.size() == this->Num_Elements());

//    // outputting result
//    // std::string pathtofile = "./work/cleanbench1.out";
//    std::string pathtofile = str_pathtofolder + "/MatrixSolution.out";
//    auto file2 = std::ofstream(pathtofile);
//    int nelem = this->Num_Elements();
//    double normfactor = this->PotentialNorm();

//    for (int i = 0; i < nelem; ++i) {

//       // transform to spatial
//       int lmax = _grid.MaxDegree();
//       // auto spatialsize = _grid.Longitudes().size() *
//       _grid.CoLatitudes().size();
//       // std::vector<std::complex<double>> vectmp(spatialsize);
//       // _grid.InverseTransformation(lmax, 0, vec_fullinformation[i][0],
//       vectmp);

//       // output
//       int idxspatial = 0;
//       for (auto it : _grid.CoLatitudes()) {
//          for (auto ip : _grid.Longitudes()) {
//             file2 << std::setprecision(16)
//                   << (node_data.ElementLowerRadius(i) +
//                       _vec_h[i][0][idxspatial]) *
//                          this->LengthNorm()
//                   << ";" << it << ";" << ip << ";"
//                   << vectmp[idxspatial].real() * normfactor << ";"
//                   << vectmp[idxspatial].imag() * normfactor;
//             if (idxspatial < spatialsize - 1) {
//                file2 << ";";
//             }
//             ++idxspatial;
//          }
//       }

//       file2 << std::endl;
//    };

//    {

//       // transform to spatial
//       int lmax = _grid.MaxDegree();
//       auto spatialsize = _grid.Longitudes().size() *
//       _grid.CoLatitudes().size(); std::vector<std::complex<double>>
//       vectmp(spatialsize); _grid.InverseTransformation(
//           lmax, 0, vec_fullinformation[nelem - 1][this->Poly_Order()],
//           vectmp);

//       // output
//       int idxspatial = 0;
//       for (auto it : _grid.CoLatitudes()) {
//          for (auto ip : _grid.Longitudes()) {
//             file2 << std::setprecision(16)
//                   << (node_data.ElementUpperRadius(nelem - 1) +
//                       _vec_h[nelem - 1][this->Poly_Order()][idxspatial]) *
//                          this->LengthNorm()
//                   << ";" << it << ";" << ip << ";"
//                   << vectmp[idxspatial].real() * normfactor << ";"
//                   << vectmp[idxspatial].imag() * normfactor;
//             if (idxspatial < spatialsize - 1) {
//                file2 << ";";
//             }
//             ++idxspatial;
//          }
//       }

//       file2 << std::endl;
//    }

//    // close
//    file2.close();
// };
void
Density3D::OutputAtElement(const std::string str_pathtofolder,
                           const std::vector<std::vector<std::complex<double>>>
                               &vec_fullinformation) const {
   assert(vec_fullinformation.size() == this->Num_Elements() + 1);
   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   std::string pathtofile = str_pathtofolder + "/IntegralSolution.out";
   auto file2 = std::ofstream(pathtofile);
   int nelem = this->Num_Elements();
   double normfactor = this->PotentialNorm();

   for (int i = 0; i < nelem; ++i) {
      double currentrad;
      if (i == nelem) {
         currentrad = this->Node_Information().ElementUpperRadius(i) *
                      this->LengthNorm();
      } else {
         currentrad = this->Node_Information().ElementLowerRadius(i) *
                      this->LengthNorm();
      }

      //
      // first output
      //
      file2 << std::setprecision(16) << currentrad;

      //
      // working through lm
      //
      int idxlm = 0;
      for (int idxl = 0; idxl < this->GSH_Grid().MaxDegree() + 1; ++idxl) {
         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
            file2 << ";" << vec_fullinformation[i][idxlm].real() * normfactor
                  << ";" << vec_fullinformation[i][idxlm].imag() * normfactor;
            ++idxlm;
         }
      }
      file2 << std::endl;
   };

   {
      double currentrad =
          this->Node_Information().ElementUpperRadius(nelem - 1) *
          this->LengthNorm();
      file2 << std::setprecision(16) << currentrad;
      int idxlm = 0;
      for (int idxl = 0; idxl < this->GSH_Grid().MaxDegree() + 1; ++idxl) {
         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
            file2 << ";"
                  << vec_fullinformation[nelem][idxlm].real() * normfactor
                  << ";"
                  << vec_fullinformation[nelem][idxlm].imag() * normfactor;
            ++idxlm;
         }
      }
      file2 << std::endl;
   }

   // close
   file2.close();
};

void
Density3D::PhysicalOutputAtElement(
    const std::string str_pathtofolder,
    const std::vector<std::vector<std::complex<double>>> &vec_fullinformation)
    const {
   assert(vec_fullinformation.size() == this->Num_Elements() + 1);
   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   std::string pathtofile = str_pathtofolder + "/IntegralSolution.out";
   auto file2 = std::ofstream(pathtofile);
   int nelem = this->Num_Elements();
   double normfactor = this->PotentialNorm();

   for (int i = 0; i < nelem; ++i) {
      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize);
      _grid.InverseTransformation(lmax, 0, vec_fullinformation[i], vectmp);

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            file2 << std::setprecision(16)
                  << (node_data.ElementLowerRadius(i) +
                      _vec_h[i][0][idxspatial]) *
                         this->LengthNorm()
                  << ";" << it << ";" << ip << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   };

   {
      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize);
      _grid.InverseTransformation(lmax, 0, vec_fullinformation.back(), vectmp);

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            file2 << std::setprecision(16)
                  << (node_data.ElementUpperRadius(nelem - 1) +
                      _vec_h[nelem - 1][this->Poly_Order()][idxspatial]) *
                         this->LengthNorm()
                  << ";" << it << ";" << ip << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   }

   // close
   file2.close();
};

void
Density3D::CartesianOutputAtElement(
    const std::string str_pathtofolder,
    const std::vector<std::vector<std::complex<double>>> &vec_fullinformation)
    const {
   assert(vec_fullinformation.size() == this->Num_Elements() + 1);
   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   std::string pathtofile = str_pathtofolder + "/IntegralSolution.out";
   auto file2 = std::ofstream(pathtofile);
   int nelem = this->Num_Elements();
   double normfactor = this->PotentialNorm();

   for (int i = 0; i < nelem; ++i) {
      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize);
      _grid.InverseTransformation(lmax, 0, vec_fullinformation[i], vectmp);

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            double physrad =
                (node_data.ElementLowerRadius(i) + _vec_h[i][0][idxspatial]) *
                this->LengthNorm();
            file2 << std::setprecision(16)
                  << physrad * std::sin(it) * std::cos(ip) << ";"
                  << physrad * std::sin(it) * std::sin(ip) << ";"
                  << physrad * std::cos(it) << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   };

   {
      // transform to spatial
      int lmax = _grid.MaxDegree();
      auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
      std::vector<std::complex<double>> vectmp(spatialsize);
      _grid.InverseTransformation(lmax, 0, vec_fullinformation.back(), vectmp);

      // output
      int idxspatial = 0;
      for (auto it : _grid.CoLatitudes()) {
         for (auto ip : _grid.Longitudes()) {
            double physrad =
                (node_data.ElementUpperRadius(nelem - 1) +
                 _vec_h[nelem - 1][this->Poly_Order()][idxspatial]) *
                this->LengthNorm();
            file2 << std::setprecision(16)
                  << physrad * std::sin(it) * std::cos(ip) << ";"
                  << physrad * std::sin(it) * std::sin(ip) << ";"
                  << physrad * std::cos(it) << ";"
                  << vectmp[idxspatial].real() * normfactor << ";"
                  << vectmp[idxspatial].imag() * normfactor;
            if (idxspatial < spatialsize - 1) {
               file2 << ";";
            }
            ++idxspatial;
         }
      }

      file2 << std::endl;
   }

   // close
   file2.close();
};

void
Density3D::OutputAtElement(const std::string str_pathtofolder,
                           const std::vector<double> &vec_exactsol) const {

   assert(vec_exactsol.size() == this->Num_Elements() + 1);
   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   std::string pathtofile = str_pathtofolder + "/ExactSolution.out";
   auto file2 = std::ofstream(pathtofile);
   int nelem = this->Num_Elements();
   double normfactor = this->PotentialNorm();

   for (int i = 0; i < nelem; ++i) {
      double currentrad;
      if (i == nelem) {
         currentrad = this->Node_Information().ElementUpperRadius(i) *
                      this->LengthNorm();
      } else {
         currentrad = this->Node_Information().ElementLowerRadius(i) *
                      this->LengthNorm();
      }
      //
      // first output
      //
      file2 << std::setprecision(16) << currentrad << ";"
            << vec_exactsol[i] * normfactor << std::endl;
   };

   {
      double currentrad =
          this->Node_Information().ElementUpperRadius(nelem - 1) *
          this->LengthNorm();
      file2 << std::setprecision(16) << currentrad << ";"
            << vec_exactsol.back() * normfactor << std::endl;
   }

   // close
   file2.close();
};
void
Density3D::OutputAtElement(const std::string str_pathtofolder,
                           const std::vector<double> &vec_radii,
                           const std::vector<double> &vec_exactsol) const {
   assert(vec_radii.size() == vec_exactsol.size());
   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   std::string pathtofile = str_pathtofolder + "/ExactSolution.out";
   auto file2 = std::ofstream(pathtofile);
   int nelem = vec_radii.size();
   double normfactor = this->PotentialNorm();

   for (int i = 0; i < nelem; ++i) {
      file2 << std::setprecision(16) << vec_radii[i] * this->LengthNorm() << ";"
            << vec_exactsol[i] * normfactor << std::endl;
   };

   // close
   file2.close();
};

void
Density3D::OutputAtElement(
    const std::string str_pathtofolder,
    const std::vector<std::vector<double>> &vec_radii,
    const std::vector<std::vector<double>> &vec_exactsol) const {
   assert(vec_radii.size() == vec_exactsol.size());
   assert(vec_radii[0].size() == vec_exactsol[0].size());

   // outputting result
   // std::string pathtofile = "./work/cleanbench1.out";
   std::string pathtofile = str_pathtofolder + "/ExactSolution.out";
   auto file2 = std::ofstream(pathtofile);
   int nspatial = vec_radii.size();
   int nelem = vec_radii[0].size();
   double normfactor = this->PotentialNorm();

   for (int i = 0; i < nelem; ++i) {
      for (int idx = 0; idx < nspatial; ++idx) {
         file2 << std::setprecision(16)
               << vec_radii[idx][i] * this->LengthNorm() << ";"
               << vec_exactsol[idx][i] * normfactor;
         if (idx < nspatial - 1) {
            file2 << ";";
         }
      }
      file2 << std::endl;
   };

   // close
   file2.close();
};
}   // namespace GeneralEarthModels
#endif