#ifndef PSEUDOSPECTRAL_OPERATOR3D_GUARD_H
#define PSEUDOSPECTRAL_OPERATOR3D_GUARD_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <FFTWpp/Ranges>
#include <GSHTrans/All>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <iostream>
// #include <Eigen/IterativeLinearSolvers>
// #include <unsupported/Eigen/IterativeSolvers>
// #include "Timer_Class.h"
// #include "Earth_Density_Models_3D.h"
// #include "Earth_General_Models_1D.h"
// #include "Radial_Tools.h"
#include "Spherical_Integrator.h"

template <typename MRScalar> class MatrixReplacement3D;
// template <typename MRScalar> class MatrixReplacement3D;
using Eigen::SparseMatrix;

namespace Eigen {
namespace internal {
// MatrixReplacement3D looks-like a SparseMatrix, so let's inherits its traits:
template <typename MRScalar>
struct traits<MatrixReplacement3D<MRScalar>>
    : public Eigen::internal::traits<Eigen::SparseMatrix<MRScalar>> {};
}   // namespace internal
}   // namespace Eigen

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
template <typename MRScalar>
class MatrixReplacement3D
    : public Eigen::EigenBase<MatrixReplacement3D<MRScalar>> {
 public:
   // Required typedefs, constants, and method:
   using Scalar = MRScalar;
   using RealScalar = double;
   using StorageIndex = int;
   using Index = Eigen::EigenBase<MatrixReplacement3D<MRScalar>>::Index;
   // using namespace GSHTrans;
   // using Real = double;
   using MRange = GSHTrans::All;
   using NRange = GSHTrans::All;
   using Grid = GSHTrans::GaussLegendreGrid<RealScalar, MRange, NRange>;
   using Quadrature = GaussQuad::Quadrature1D<RealScalar>;
   using MATRIX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
   using MATRIX3 = Eigen::Matrix<Scalar, 3, 3>;
   using vvveceig = std::vector<std::vector<std::vector<Eigen::Matrix3cd>>>;

   // typedef double Scalar;
   // typedef double RealScalar;
   // typedef int StorageIndex;
   enum {
      ColsAtCompileTime = Eigen::Dynamic,
      MaxColsAtCompileTime = Eigen::Dynamic,
      IsRowMajor = false
   };

   // Index rows() const { return this->my_matrix().rows(); }
   // Index cols() const { return this->my_matrix().cols(); }
   Index rows() const {
      return this->_indexclass.matrixdimension(this->nelem());
   }
   Index cols() const {
      return this->_indexclass.matrixdimension(this->nelem());
   }

   template <typename Rhs>
   Eigen::Product<MatrixReplacement3D<Scalar>, Rhs, Eigen::AliasFreeProduct>
   operator*(const Eigen::MatrixBase<Rhs> &x) const {
      return Eigen::Product<MatrixReplacement3D<Scalar>, Rhs,
                            Eigen::AliasFreeProduct>(*this, x.derived());
   }

   // constructor with 3D class
   MatrixReplacement3D(GeneralEarthModels::Density3D &inp_model)
       : _inp_model(&inp_model),
         _indexclass(0, inp_model.GSH_Grid().MaxDegree(),
                     inp_model.Poly_Order()),
         _vec_a(inp_model.LaplaceTensorRef()) {};

   ////////////////////////////////////////////////////////////////////
   //////////////////////  CONSTRUCTOR WITH MAP ///////////////////////
   ////////////////////////////////////////////////////////////////////

   // will need constructor, storage for vec_a and maybe a return for the vector
   // depending upon a boolean that changes depending upon the constructor
   // constructor with 3D class

   MatrixReplacement3D(
       GeneralEarthModels::Density3D &inp_model,
       std::vector<std::vector<std::vector<Eigen::Matrix3cd>>> &vec_a)
       : _inp_model(&inp_model),
         _indexclass(0, inp_model.GSH_Grid().MaxDegree(),
                     inp_model.Poly_Order()),
         _separate_map{true}, _vec_a{vec_a} {};

   MatrixReplacement3D(GeneralEarthModels::Density3D &inp_model,
                       const GeneralEarthModels::MappingPerturbation &inp_map)
       : _inp_model(&inp_model),
         _indexclass(0, inp_model.GSH_Grid().MaxDegree(),
                     inp_model.Poly_Order()),
         _separate_map{true}, _vec_a{inp_map.ref_da()} {};

   //  return model
   auto Model() const { return _inp_model; }
   const std::vector<std::vector<std::vector<Eigen::Matrix3cd>>> &
   matrix_A() const {
      return _vec_a;
   }

   int nelem() const { return _inp_model->Num_Elements(); }
   int npoly() const { return _inp_model->Poly_Order(); }

   void IncludeBoundary() const { _inc_boundterm = true; }
   void NoBoundary() const { _inc_boundterm = false; }
   bool BoundaryStatus() const { return _inc_boundterm; }

   std::size_t myidxfunc(int idxelem, int idxpoly, int idxl, int idxm) const {
      return static_cast<int>(
          std::pow(idxl, 2) + idxl + idxm +
          (idxpoly + npoly() * idxelem) *
              std::pow(_inp_model->GSH_Grid().MaxDegree() + 1, 2));
   };
   std::size_t myidxpm(int idxelem, int idxpoly, int idxl, int idxm) const {
      return static_cast<int>(
          std::pow(idxl, 2) + idxl + idxm - 1 +
          (idxpoly + npoly() * idxelem) *
              (std::pow(_inp_model->GSH_Grid().MaxDegree() + 1, 2) - 1));
   };
   std::size_t idxsmall(int idxl, int idxm) const {
      return static_cast<int>(idxl * idxl + idxl + idxm);
   };
   std::size_t idxpmsmall(int idxl, int idxm) const {
      return static_cast<int>(idxl * idxl + idxl + idxm - 1);
   };

 private:
   // reference to Density3D model
   const GeneralEarthModels::Density3D *_inp_model;
   const SpectralElementTools::MatrixIndices _indexclass;
   const std::vector<std::vector<std::vector<Eigen::Matrix3cd>>> &_vec_a;
   // Grid _grid;
   mutable bool _inc_boundterm = true;
   mutable bool _separate_map = false;
};

// Implementation of MatrixReplacement3D * Eigen::DenseVector though a
// specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {

template <typename MRScalar, typename Rhs>
struct generic_product_impl<MatrixReplacement3D<MRScalar>, Rhs, SparseShape,
                            DenseShape,
                            GemvProduct>   // GEMV stands for matrix-vector
    : generic_product_impl_base<
          MatrixReplacement3D<MRScalar>, Rhs,
          generic_product_impl<MatrixReplacement3D<MRScalar>, Rhs>> {
   // typedef typename Product<MatrixReplacement3D<MRScalar>, Rhs>::Scalar
   // Scalar;
   using Scalar = Product<MatrixReplacement3D<MRScalar>, Rhs>::Scalar;
   using Real = Product<MatrixReplacement3D<MRScalar>, Rhs>::RealScalar;
   using Complex = std::complex<Real>;

   // if we have a vector x, we are defining x + alpha * lhs * rhs
   template <typename Dest>
   static void scaleAndAddTo(Dest &dst,
                             const MatrixReplacement3D<MRScalar> &lhs,
                             const Rhs &rhs, const Scalar &alpha) {
      // This method implements "dst += alpha * lhs * rhs" inplace

      // using statements and declarations
      using MATRIX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
      using MATRIX3 = Eigen::Matrix<Scalar, 3, 3>;
      using VECTOR = Eigen::Vector<Complex, Eigen::Dynamic>;
      VECTOR vec_output2 = VECTOR::Zero(lhs.rows());

      int tcount = 0;
      int laynum = 0;
      int lMax = lhs.Model()->GSH_Grid().MaxDegree();
      //   std::cout << "\n\n lMax: " << lMax << "\n\n";
      int size0 = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 0).size();
      int sizepm = GSHTrans::GSHIndices<GSHTrans::All>(lMax, lMax, 1).size();
      int sizespatial = lhs.Model()->GSH_Grid().NumberOfLongitudes() *
                        lhs.Model()->GSH_Grid().NumberOfCoLatitudes();
      //   std::vector<std::vector<std::vector<Eigen::Matrix3cd>>> _vec_a =
      //       lhs.Model()->LaplaceTensor();
      // std::cout << _vec_a[0][0][0] << "\n";
      // std::cout << _vec_a[1][0][0] << "\n";
      for (int idxelem = 0; idxelem < lhs.nelem(); ++idxelem) {
         double inv2 =
             2.0 / lhs.Model()->Node_Information().ElementWidth(idxelem);

         // looping over quadrature nodes
         int idxpolymax = lhs.npoly() + 1;
         for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
            // radius, inverse and squared
            auto rval =
                lhs.Model()->Node_Information().NodeRadius(idxelem, idxpoly);
            auto invr = 1 / rval;
            auto rval2 = rval * rval;

            ///////////////////////////////////////
            // step 1: finding nabla zeta
            // std::cout << "Step 1a\n";
            std::vector<Scalar> gsph_nz0(size0, 0.0), gsph_nzp1(sizepm, 0.0),
                gsph_nzm1(sizepm, 0.0);   // declaration

            // looping over l,m values
            {

               // cache friendly trial
               //  getting first index
               std::size_t idx1 = lhs.myidxfunc(idxelem, 0, 0, 0);

               // looping over radii
               for (int idxn = 0; idxn < lhs.npoly() + 1; ++idxn) {
                  std::size_t idx2 = 0;
                  auto multfact =
                      lhs.Model()->GaussDerivative(idxn, idxpoly) * rval * inv2;

                  // looping over l and m
                  for (int idxl = 0; idxl < lMax + 1; ++idxl) {
                     for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {

                        // increment local temporary
                        gsph_nz0[idx2] += rhs(idx1) * multfact;
                        ++idx1;
                        ++idx2;
                     }
                  }
               }
            }

            {
               std::size_t idx1 = 0;
               std::size_t idx2 = lhs.myidxfunc(idxelem, idxpoly, 1, -1);
               for (int idxl = 1; idxl < lMax + 1; ++idxl) {
                  // omega_l^0
                  double Omegal0 =
                      std::sqrt(static_cast<double>(idxl) *
                                static_cast<double>(idxl + 1) / 2.0);
                  auto multfact = Omegal0;
                  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {

                     auto tmp = rhs(idx2) * multfact;
                     gsph_nzp1[idx1] += tmp;
                     gsph_nzm1[idx1] += tmp;

                     // increment
                     ++idx1;
                     ++idx2;
                  }
               }
            }

            ///////////////////////////////////////////////////////////////////////////
            // step 2: transforming (\nabla \zeta) into spatial domain
            // std::cout << "Step 2a\n";
            FFTWpp::vector<Complex> spatial_nzm1(sizespatial, 0.0),
                spatial_nzp1(sizespatial, 0.0),
                spatial_nz0(sizespatial, 0.0);   // declaration

            // GSPH transforms:
            lhs.Model()->GSH_Grid().InverseTransformation(
                lMax, 0, gsph_nz0,
                spatial_nz0);   // 0th order
            // std::cout << "Between transformations\n";
            lhs.Model()->GSH_Grid().InverseTransformation(
                lMax, 1, gsph_nzp1,
                spatial_nzp1);   //+1 order
            lhs.Model()->GSH_Grid().InverseTransformation(
                lMax, -1, gsph_nzm1,
                spatial_nzm1);   //-1 order

            /////////////////////////////////////////////////////////////
            // step 3: finding q = a\nabla \zeta
            // declaration
            // std::cout << "Step 3a\n";
            FFTWpp::vector<Scalar> spatial_qm(sizespatial, 0.0),
                spatial_qp(sizespatial, 0.0), spatial_q0(sizespatial, 0.0);

            // multiplying through by matrix a
            {
               for (int idxr = 0; idxr < sizespatial; ++idxr) {

                  //   spatial_qm[idxr] -=
                  //       _vec_a[idxelem][idxpoly][idxr](0, 0) *
                  //       spatial_nzp1[idxr];
                  spatial_qm[idxr] -=
                      lhs.matrix_A()[idxelem][idxpoly][idxr](0, 0) *
                      spatial_nzp1[idxr];
                  spatial_qm[idxr] +=
                      lhs.matrix_A()[idxelem][idxpoly][idxr](0, 1) *
                      spatial_nz0[idxr];
                  spatial_qm[idxr] -=
                      lhs.matrix_A()[idxelem][idxpoly][idxr](0, 2) *
                      spatial_nzm1[idxr];

                  spatial_q0[idxr] -=
                      lhs.matrix_A()[idxelem][idxpoly][idxr](1, 0) *
                      spatial_nzp1[idxr];
                  spatial_q0[idxr] +=
                      lhs.matrix_A()[idxelem][idxpoly][idxr](1, 1) *
                      spatial_nz0[idxr];
                  spatial_q0[idxr] -=
                      lhs.matrix_A()[idxelem][idxpoly][idxr](1, 2) *
                      spatial_nzm1[idxr];

                  spatial_qp[idxr] -=
                      lhs.matrix_A()[idxelem][idxpoly][idxr](2, 0) *
                      spatial_nzp1[idxr];
                  spatial_qp[idxr] +=
                      lhs.matrix_A()[idxelem][idxpoly][idxr](2, 1) *
                      spatial_nz0[idxr];
                  spatial_qp[idxr] -=
                      lhs.matrix_A()[idxelem][idxpoly][idxr](2, 2) *
                      spatial_nzm1[idxr];
               }
            }

            /////////////////////////////////////////////////////////////
            // step 4: transforming q from spatial to spherical harmonic
            // basis
            //  transformation for q^{-1}
            // std::cout << "Step 4a\n";
            std::vector<Complex> gsph_q0(size0, 0.0), gsph_qp1(sizepm, 0.0),
                gsph_qm1(sizepm, 0.0);   // declaration of GSPH vectors
            {
               // transformation for q^{-1}
               lhs.Model()->GSH_Grid().ForwardTransformation(
                   lMax, -1, spatial_qm, gsph_qm1);

               // transformation for q^0
               lhs.Model()->GSH_Grid().ForwardTransformation(
                   lMax, 0, spatial_q0, gsph_q0);

               // transformation for q^{+1}
               lhs.Model()->GSH_Grid().ForwardTransformation(
                   lMax, +1, spatial_qp, gsph_qp1);
            }

            /////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////
            // step 5: evaluating radial integrals to give Ax
            // std::cout << "Step 5a\n";
            // 0 component term
            {
               // multiplication factor and index
               std::size_t idx2 = lhs.myidxfunc(idxelem, 0, 0, 0);
               // auto mult1 = lhs.gaussweight(idxpoly) * rval2;
               auto mult1 = lhs.Model()->q().W(idxpoly) * rval;

               // loop over radii
               for (int idxn = 0; idxn < lhs.npoly() + 1; ++idxn) {

                  // index and multiplication factor
                  std::size_t idx3 = 0;
                  // auto multfact = lhs.polyderiv(idxn, idxpoly) * mult1;
                  auto multfact =
                      lhs.Model()->GaussDerivative(idxn, idxpoly) * mult1;
                  // loop over l and m
                  for (int idxl = 0; idxl < lMax + 1; ++idxl) {
                     for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                        // increment global
                        vec_output2(idx2) -= gsph_q0[idx3] * multfact;

                        // increment indices
                        ++idx2;
                        ++idx3;
                     }
                  }
               }
            }
            // std::cout << "Step 5b\n";
            // pm component
            {
               std::size_t idx1 = lhs.myidxfunc(idxelem, idxpoly, 1, -1);
               std::size_t idx3 = 0;
               // auto mult1 = lhs.elemwidth(idxelem) * 0.5 *
               //              lhs.gaussweight(idxpoly) * rval;
               auto mult1 =
                   lhs.Model()->Node_Information().ElementWidth(idxelem) * 0.5 *
                   lhs.Model()->q().W(idxpoly);
               for (int idxl = 1; idxl < lMax + 1; ++idxl) {
                  // multiplication factors
                  double omegal0 = std::sqrt(static_cast<double>(idxl) *
                                             (static_cast<double>(idxl + 1)) *
                                             0.5);   // omega_l^0
                  auto multfact = omegal0 * mult1;

                  // loop over m
                  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                     vec_output2(idx1) -=
                         multfact * (gsph_qm1[idx3] +
                                     gsph_qp1[idx3]);   // global increment

                     // increment indices
                     ++idx1;
                     ++idx3;
                  }
               }
            }
            // std::cout << "Step 5c\n";
            // external term
            if (lhs.BoundaryStatus()) {
               if (idxelem == lhs.nelem() - 1 && idxpoly == lhs.npoly()) {
                  for (int idxl = 0; idxl < lMax + 1; ++idxl) {

                     for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                        // index of vector
                        std::size_t idx1 =
                            lhs.myidxfunc(idxelem, idxpoly, idxl, idxm);

                        vec_output2(idx1) -=
                            (static_cast<double>(idxl) + 1.0) *
                            lhs.Model()->Node_Information().OuterRadius() *
                            rhs(idx1);
                     }
                  }
               }
            }

            if (idxelem)
               ++tcount;
         };
      };

      // final return
      dst.noalias() += alpha * vec_output2;
   }
};
}   // namespace internal
}   // namespace Eigen
#endif