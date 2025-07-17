#ifndef PSEUDOSPECTRAL_OPERATOR_GUARD_H
#define PSEUDOSPECTRAL_OPERATOR_GUARD_H

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
#include "Radial_Tools.h"
#include "Spherical_Integrator.h"

template <typename MRScalar> class MatrixReplacement;
// template <typename MRScalar> class MatrixReplacement;
using Eigen::SparseMatrix;

namespace Eigen {
namespace internal {
// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
template <typename MRScalar>
struct traits<MatrixReplacement<MRScalar>>
    : public Eigen::internal::traits<Eigen::SparseMatrix<MRScalar>> {};
}   // namespace internal
}   // namespace Eigen

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
template <typename MRScalar>
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<MRScalar>> {
 public:
   // Required typedefs, constants, and method:
   using Scalar = MRScalar;
   using RealScalar = double;
   using StorageIndex = int;
   using Index = Eigen::EigenBase<MatrixReplacement<MRScalar>>::Index;
   // using namespace GSHTrans;
   // using Real = double;
   using MRange = GSHTrans::All;
   using NRange = GSHTrans::All;
   using Grid = GSHTrans::GaussLegendreGrid<RealScalar, MRange, NRange>;
   using Quadrature = GaussQuad::Quadrature1D<RealScalar>;
   using MATRIX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
   using MATRIX3 = Eigen::Matrix<Scalar, 3, 3>;

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
   Index rows() const { return this->_matlen * this->_size0; }
   Index cols() const { return this->_matlen * this->_size0; }

   template <typename Rhs>
   Eigen::Product<MatrixReplacement<Scalar>, Rhs, Eigen::AliasFreeProduct>
   operator*(const Eigen::MatrixBase<Rhs> &x) const {
      return Eigen::Product<MatrixReplacement<Scalar>, Rhs,
                            Eigen::AliasFreeProduct>(*this, x.derived());
   }

   // Custom API:
   // MatrixReplacement() {}

   // void attachMyMatrix(const SparseMatrix<Scalar> &mat)
   // {
   //     mp_mat = &mat;
   // }
   // const SparseMatrix<Scalar> my_matrix() const { return *mp_mat; }

   // constructor with my new spherical 1D class
   MatrixReplacement(GeneralEarthModels::spherical_1D &inp_model, Grid &grid)
       : _grid{&grid} {

      // useful numbers
      auto intsize = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
      // std::cout << "intsize: " << intsize << "\n";
      int nelem = inp_model.Num_Elements();
      int npoly = inp_model.Poly_Order();
      std::vector<RealScalar> vec_noderadii =
          inp_model.Node_Information().ElementNodes();
      Quadrature tmp_q = inp_model.q();

      // name abbreviations
      // using MATRIX3 = Eigen::Matrix<std::complex<double>, 3, 3>;
      using vecM = std::vector<MATRIX3>;
      using vvecM = std::vector<vecM>;

      // defining vec_a, ie vector of all a's
      MATRIX3 mat_0;
      mat_0 << 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0;
      // std::cout << mat_0 << "\n";
      vvecM tmp_a(nelem * (npoly + 1), vecM(intsize, mat_0));

      // setting up everything
      defineNodes(vec_noderadii);
      defineGrid(grid);
      defineQuadrature(tmp_q);
      defineAsphericity(tmp_a);
      // std::cout << "Hello\n";
   };

   // constructor with 3D class
   MatrixReplacement(GeneralEarthModels::Density3D &inp_model) {
      // We have grid, q, vec_a, only need vec_noderadii
      std::vector<double> vec_noderadii;
      vec_noderadii.push_back(0.0);
      for (int idx = 0; idx < inp_model.Node_Information().NumberOfElements();
           ++idx) {
         vec_noderadii.push_back(
             inp_model.Node_Information().ElementUpperRadius(idx));
      }
      defineNodes(vec_noderadii);
      defineGrid(inp_model.GSH_Grid());
      defineQuadrature(inp_model.q());
      defineAsphericity(inp_model.LaplaceTensor());
   };

   // constructor that takes arguments
   MatrixReplacement(const std::vector<RealScalar> &vec_noderadii,
                     const Grid &grid, const Quadrature &q,
                     const std::vector<std::vector<MATRIX3>> &vec_a)
       : _grid{&grid} {
      // use constituent functions
      defineNodes(vec_noderadii);
      defineGrid(grid);
      defineQuadrature(q);
      defineAsphericity(vec_a);
   };

   // add data outside of constructor
   void addData(const std::vector<RealScalar> &vec_noderadii, const Grid &grid,
                const Quadrature &q,
                const std::vector<std::vector<MATRIX3>> &vec_a) {
      // use constituent functions
      defineNodes(vec_noderadii);
      defineGrid(grid);
      defineQuadrature(q);
      defineAsphericity(vec_a);
   }

   // attaching reference to node radii
   void defineNodes(const std::vector<RealScalar> &vec_noderadii) {
      // std::cout << "Define nodes\n";
      _vec_noderadii = vec_noderadii;
      _nelem = vec_noderadii.size() - 1;
      _matlen = _nelem * _npoly + 1;
      _vec_elem_width.resize(_nelem);
      for (int idx = 0; idx < _nelem; ++idx) {
         _vec_elem_width[idx] = _vec_noderadii[idx + 1] - _vec_noderadii[idx];
      }
   }
   const std::vector<RealScalar> noderadii() const { return _vec_noderadii; }
   RealScalar elemwidth(int i) const { return _vec_elem_width[i]; }
   RealScalar node(int i) const { return _vec_noderadii[i]; }
   int nelem() const { return _nelem; }
   int npoly() const { return _npoly; }

   // attach reference to grid
   void defineGrid(const Grid &grid) {
      // _grid = &grid;
      // std::cout << "Define grid\n";
      _lMax = grid.MaxDegree();
      _nMax = grid.MaxUpperIndex();
      _size0 = GSHTrans::GSHIndices<GSHTrans::All>(_lMax, _lMax, 0).size();
      _sizepm = GSHTrans::GSHIndices<GSHTrans::All>(_lMax, _lMax, 1).size();
      _sizespatial = grid.NumberOfLongitudes() * grid.NumberOfCoLatitudes();
   }
   const Grid gridused() const { return *_grid; }
   const Grid *gridusedp() const { return _grid; }
   int lMax() const { return _lMax; }
   std::size_t gridsize() const { return _sizespatial; }
   std::size_t size0() const { return _size0; }
   std::size_t sizepm() const { return _sizepm; }

   // reference to Gaussian quadrature
   void defineQuadrature(const Quadrature &q) {
      // std::cout << "Define quadrature\n";
      _q = q;
      _npoly = q.Points().size() - 1;
      for (int idx = 0; idx < q.Points().size(); ++idx) {
         _vec_gaussquadpoints.push_back(q.X(idx));
         _vec_gaussquadweights.push_back(q.W(idx));
      }
      _matlen = _nelem * _npoly + 1;
      _mat_gaussderiv.resize(_npoly + 1, _npoly + 1);
      {
         auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
                                                       _q.Points().end());
         for (int idxi = 0; idxi < _npoly + 1; ++idxi) {
            for (int idxj = 0; idxj < _npoly + 1; ++idxj) {
               _mat_gaussderiv(idxi, idxj) = pleg.Derivative(idxi, _q.X(idxj));
            }
         }
      }
   }
   const Quadrature q() const { return _q; }
   RealScalar polyderiv(int idxi, int idxj) const {
      return _mat_gaussderiv(idxi, idxj);
      // returns the derivative of the ith interpolating polynomial at the jth
      // quadrature point
   }
   RealScalar gausspoint(int idxi) const { return _vec_gaussquadpoints[idxi]; }
   RealScalar gaussweight(int idxi) const {
      return _vec_gaussquadweights[idxi];
   }

   // matrix a
   void defineAsphericity(const std::vector<std::vector<MATRIX3>> &vec_a) {
      _vec_a = vec_a;
      _vec_a2 = &vec_a;
      // assert(vec_a.size() == this->_nelem * (this->_npoly + 1) &&
      //  "Incorrect size of vec_a");
   }

   void defineAsphericity(
       const std::vector<std::vector<std::vector<MATRIX3>>> &vec_a) {
      // std::cout << vec_a.size() << " " << vec_a[0].size() << "\n";
      for (int idxelem = 0; idxelem < vec_a.size(); ++idxelem) {
         for (int idxnode = 0; idxnode < vec_a[idxelem].size(); ++idxnode) {
            _vec_a.push_back(vec_a[idxelem][idxnode]);
         }
      }
      _vec_a2 = &_vec_a;
      // std::cout << _vec_a.size() << " " << _vec_a[0].size() << "\n";
   }
   std::vector<std::vector<MATRIX3>> total_a() const { return _vec_a; }
   std::vector<MATRIX3> grid_a(int idxr) const { return _vec_a[idxr]; }
   MATRIX3 mat_a(int idxradius, int idxgrid) const {
      return _vec_a[idxradius][idxgrid];
   }
   const std::vector<std::vector<MATRIX3>> *total_ap() const { return _vec_a2; }

   void IncludeBoundary() const { _inc_boundterm = true; }
   void NoBoundary() const { _inc_boundterm = false; }
   bool BoundaryStatus() const { return _inc_boundterm; }

   std::size_t myidxfunc(int idxelem, int idxpoly, int idxl, int idxm) const {
      return static_cast<int>(std::pow(idxl, 2) + idxl + idxm +
                              (idxpoly + this->_npoly * idxelem) *
                                  std::pow(this->_lMax + 1, 2));
   };
   std::size_t myidxpm(int idxelem, int idxpoly, int idxl, int idxm) const {
      return static_cast<int>(std::pow(idxl, 2) + idxl + idxm - 1 +
                              (idxpoly + this->_npoly * idxelem) *
                                  (std::pow(this->_lMax + 1, 2) - 1));
   };
   std::size_t idxsmall(int idxl, int idxm) const {
      return static_cast<int>(idxl * idxl + idxl + idxm);
   };
   std::size_t idxpmsmall(int idxl, int idxm) const {
      return static_cast<int>(idxl * idxl + idxl + idxm - 1);
   };

 private:
   // const SparseMatrix<Scalar> *mp_mat;
   std::vector<RealScalar> _vec_noderadii, _vec_elem_width;
   std::vector<RealScalar> _vec_gaussquadpoints, _vec_gaussquadweights;
   const Grid *_grid;
   Quadrature _q;
   std::size_t _nelem{0}, _npoly{0}, _matlen{0};
   long int _lMax{0}, _nMax{0}, _size0{0}, _sizepm{0}, _sizespatial{0};
   std::vector<std::vector<MATRIX3>> _vec_a;
   const std::vector<std::vector<MATRIX3>> *_vec_a2;
   // Interpolation::LagrangePolynomial _pleg;
   Eigen::Matrix<RealScalar, Eigen::Dynamic, Eigen::Dynamic> _mat_gaussderiv;
   mutable bool _inc_boundterm = true;
};

// Implementation of MatrixReplacement * Eigen::DenseVector though a
// specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {

template <typename MRScalar, typename Rhs>
struct generic_product_impl<MatrixReplacement<MRScalar>, Rhs, SparseShape,
                            DenseShape,
                            GemvProduct>   // GEMV stands for matrix-vector
    : generic_product_impl_base<
          MatrixReplacement<MRScalar>, Rhs,
          generic_product_impl<MatrixReplacement<MRScalar>, Rhs>> {
   // typedef typename Product<MatrixReplacement<MRScalar>, Rhs>::Scalar
   // Scalar;
   using Scalar = Product<MatrixReplacement<MRScalar>, Rhs>::Scalar;
   using Real = Product<MatrixReplacement<MRScalar>, Rhs>::RealScalar;
   using Complex = std::complex<Real>;

   // if we have a vector x, we are defining x + alpha * lhs * rhs
   template <typename Dest>
   static void scaleAndAddTo(Dest &dst, const MatrixReplacement<MRScalar> &lhs,
                             const Rhs &rhs, const Scalar &alpha) {
      // This method implements "dst += alpha * lhs * rhs" inplace

      // using statements and declarations
      using MATRIX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
      using MATRIX3 = Eigen::Matrix<Scalar, 3, 3>;
      using VECTOR = Eigen::Vector<Complex, Eigen::Dynamic>;
      VECTOR vec_output2 = VECTOR::Zero(lhs.rows());

      int tcount = 0;
      int laynum = 0;
      for (int idxelem = 0; idxelem < lhs.nelem(); ++idxelem) {
         // std::cout << "idxelem: " << idxelem << "\n";
         double inv2 = 2.0 / lhs.elemwidth(idxelem);
         // looping over quadrature nodes
         int idxpolymax = lhs.npoly() + 1;
         for (int idxpoly = 0; idxpoly < idxpolymax; ++idxpoly) {
            // radius, inverse and squared
            auto rval = GravityFunctions::StandardIntervalMap(
                lhs.gausspoint(idxpoly), lhs.node(idxelem),
                lhs.node(idxelem + 1));
            auto invr = 1 / rval;
            auto rval2 = rval * rval;
            // if (idxelem == 106 || idxelem == 107) {
            //    std::cout << rval << " " << lhs.elemwidth(idxelem) << "\n";
            // }

            ///////////////////////////////////////
            // step 1: finding nabla zeta
            // std::cout << "Step 1a\n";
            std::vector<Scalar> gsph_nz0(lhs.size0(), 0.0),
                gsph_nzp1(lhs.sizepm(), 0.0),
                gsph_nzm1(lhs.sizepm(), 0.0);   // declaration

            // looping over l,m values
            {

               // cache friendly trial
               //  getting first index
               std::size_t idx1 = lhs.myidxfunc(idxelem, 0, 0, 0);

               // looping over radii
               for (int idxn = 0; idxn < lhs.npoly() + 1; ++idxn) {
                  std::size_t idx2 = 0;
                  auto multfact = lhs.polyderiv(idxn, idxpoly) * inv2;
                  // if (idxelem == 106 && idxpoly == 0) {
                  //    std::cout << std::setprecision(16) << multfact * rval
                  //              << "\n";
                  // }
                  // looping over l and m
                  for (int idxl = 0; idxl < lhs.lMax() + 1; ++idxl) {
                     for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {

                        // increment local temporary
                        gsph_nz0[idx2] += rhs(idx1) * multfact;

                        // if (idxm == 0 && idxl == 0 && idxelem == 106 &&
                        //     idxpoly == 0) {
                        //    std::cout << idx2 << ": " << std::setprecision(16)
                        //              << gsph_nz0[idx2] * rval << " "
                        //              << rhs(idx1) * multfact * rval << "\n";
                        // }
                        // increment indices
                        ++idx1;
                        ++idx2;
                     }
                  }
               }
               // if (idxelem == 106 && idxpoly == 0) {
               //    int idxgsph = 0;
               //    for (int idxl = 0; idxl < lhs.lMax() + 1; ++idxl) {
               //       for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
               //          if (idxm == 0) {
               //             std::cout << idxgsph << ": " <<
               //             std::setprecision(16)
               //                       << gsph_nz0[idxgsph] * rval << "\n";
               //          }
               //          ++idxgsph;
               //       }
               //    }
               // }
               // for (auto &idx : gsph_nz0) {
               //    idx *= inv2;
               // }
            }
            // std::cout << "Step 1b\n";
            // only do the pm if not at the zero radius
            // if (tcount > 0) {
            {
               std::size_t idx1 = 0;
               std::size_t idx2 = lhs.myidxfunc(idxelem, idxpoly, 1, -1);
               for (int idxl = 1; idxl < lhs.lMax() + 1; ++idxl) {
                  // omega_l^0
                  double Omegal0 =
                      std::sqrt(static_cast<double>(idxl) *
                                static_cast<double>(idxl + 1) / 2.0);
                  auto multfact = Omegal0 * invr;
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
            // }
            // if (idxelem == 106) {
            //    std::cout << std::setprecision(12) << gsph_nzm1[0] * rval << "
            //    "
            //              << gsph_nz0[0] * rval << " " << gsph_nzp1[0] * rval
            //              << "\n";
            // }
            ///////////////////////////////////////////////////////////////////////////
            // step 2: transforming (\nabla \zeta) into spatial domain
            // std::cout << "Step 2a\n";
            FFTWpp::vector<Complex> spatial_nzm1(lhs.gridsize(), 0.0),
                spatial_nzp1(lhs.gridsize(), 0.0),
                spatial_nz0(lhs.gridsize(), 0.0);   // declaration

            // std::cout << "Gridsize: " << lhs.gridsize()
            //           << ", lmax: " << lhs.lMax() << "\n";
            // GSPH transforms:
            // lhs.gridused().InverseTransformation(lhs.lMax(), 0, gsph_nz0,
            //                                      spatial_nz0);   // 0th order
            lhs.gridusedp()->InverseTransformation(lhs.lMax(), 0, gsph_nz0,
                                                   spatial_nz0);   // 0th order
            // std::cout << "Between transformations\n";
            lhs.gridusedp()->InverseTransformation(lhs.lMax(), 1, gsph_nzp1,
                                                   spatial_nzp1);   //+1 order
            lhs.gridusedp()->InverseTransformation(lhs.lMax(), -1, gsph_nzm1,
                                                   spatial_nzm1);   //-1 order
            // lhs.gridused().InverseTransformation(lhs.lMax(), 1, gsph_nzp1,
            //                                      spatial_nzp1);   //+1 order
            // lhs.gridused().InverseTransformation(lhs.lMax(), -1, gsph_nzm1,
            //                                      spatial_nzm1);   //-1 order

            // if (idxelem == 106) {
            //    // for (auto &idxtrial : spatial_nz0) {
            //    std::cout << std::setprecision(12) << spatial_nzp1[0] * rval
            //              << " " << spatial_nz0[0] * rval << " "
            //              << spatial_nzp1[0] * rval << "\n";

            //    // }
            // }

            /////////////////////////////////////////////////////////////
            // step 3: finding q = a\nabla \zeta
            // declaration
            // std::cout << "Step 3a\n";
            FFTWpp::vector<Scalar> spatial_qm(lhs.gridsize(), 0.0),
                spatial_qp(lhs.gridsize(), 0.0),
                spatial_q0(lhs.gridsize(), 0.0);

            // multiplying through by matrix a
            {
               auto matreturn = lhs.total_ap();
               auto &matuse = *matreturn;
               std::size_t idxext = idxelem * (lhs.npoly() + 1) + idxpoly;
               // std::cout << "gridsize(): " << lhs.gridsize() << "\n";
               // std::cout << "Hello\n";
               for (int idxr = 0; idxr < lhs.gridsize(); ++idxr) {
                  // std::cout << "Before declaration of mat_a: " << idxr <<
                  // "\n";
                  // MATRIX3 mat_a = matuse[idxext][idxr];
                  MATRIX3 mat_a = lhs.mat_a(idxext, idxr);
                  // std::cout << "After declaration of mat_a: " << idxr <<
                  // "\n";
                  spatial_qm[idxr] -= mat_a(0, 0) * spatial_nzp1[idxr];
                  spatial_qm[idxr] += mat_a(0, 1) * spatial_nz0[idxr];
                  spatial_qm[idxr] -= mat_a(0, 2) * spatial_nzm1[idxr];

                  spatial_q0[idxr] -= mat_a(1, 0) * spatial_nzp1[idxr];
                  spatial_q0[idxr] += mat_a(1, 1) * spatial_nz0[idxr];
                  spatial_q0[idxr] -= mat_a(1, 2) * spatial_nzm1[idxr];

                  spatial_qp[idxr] -= mat_a(2, 0) * spatial_nzp1[idxr];
                  spatial_qp[idxr] += mat_a(2, 1) * spatial_nz0[idxr];
                  spatial_qp[idxr] -= mat_a(2, 2) * spatial_nzm1[idxr];
               }
            }

            /////////////////////////////////////////////////////////////
            // step 4: transforming q from spatial to spherical harmonic
            // basis
            //  transformation for q^{-1}
            // std::cout << "Step 4a\n";
            std::vector<Complex> gsph_q0(lhs.size0(), 0.0),
                gsph_qp1(lhs.sizepm(), 0.0),
                gsph_qm1(lhs.sizepm(), 0.0);   // declaration of GSPH vectors
            {
               // transformation for q^{-1}

               lhs.gridusedp()->ForwardTransformation(lhs.lMax(), -1,
                                                      spatial_qm, gsph_qm1);
               // lhs.gridused().ForwardTransformation(lhs.lMax(), -1,
               // spatial_qm,
               //                                      gsph_qm1);

               // transformation for q^0
               // lhs.gridused().ForwardTransformation(lhs.lMax(), 0,
               // spatial_q0,
               //                                      gsph_q0);
               lhs.gridusedp()->ForwardTransformation(lhs.lMax(), 0, spatial_q0,
                                                      gsph_q0);

               // transformation for q^{+1}
               // lhs.gridused().ForwardTransformation(lhs.lMax(), +1,
               // spatial_qp,
               //                                      gsph_qp1);
               lhs.gridusedp()->ForwardTransformation(lhs.lMax(), +1,
                                                      spatial_qp, gsph_qp1);
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
               auto mult1 = lhs.gaussweight(idxpoly) * rval2;

               // loop over radii
               for (int idxn = 0; idxn < lhs.npoly() + 1; ++idxn) {

                  // index and multiplication factor
                  std::size_t idx3 = 0;
                  auto multfact = lhs.polyderiv(idxn, idxpoly) * mult1;

                  // loop over l and m
                  for (int idxl = 0; idxl < lhs.lMax() + 1; ++idxl) {
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
               auto mult1 = lhs.elemwidth(idxelem) * 0.5 *
                            lhs.gaussweight(idxpoly) * rval;
               for (int idxl = 1; idxl < lhs.lMax() + 1; ++idxl) {

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
                  for (int idxl = 0; idxl < lhs.lMax() + 1; ++idxl) {

                     for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                        // index of vector
                        std::size_t idx1 =
                            lhs.myidxfunc(idxelem, idxpoly, idxl, idxm);

                        vec_output2(idx1) -= (static_cast<double>(idxl) + 1.0) *
                                             lhs.node(lhs.nelem()) * rhs(idx1);
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