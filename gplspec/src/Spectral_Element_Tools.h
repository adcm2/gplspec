#ifndef SPECTRAL_ELEMENT_TOOLS_H
#define SPECTRAL_ELEMENT_TOOLS_H

#include "Radial_Tools.h"
#include <GaussQuad/All>
#include <Interpolation/All>

namespace SpectralElementTools {
class MatrixIndices {
 public:
   MatrixIndices() {};
   MatrixIndices(int lMin, int lMax, int npoly)
       : _lMax{lMax}, _lMin{lMin}, _npoly{npoly} {
      _linlen = (lMax + 1 + lMin) * (lMax + 1 - lMin);
      _elemlen = npoly * _linlen;   // number of points in first npoly rows
   };

   // output
   int linelength() const { return this->_linlen; };
   std::size_t fullindex(int idxelem, int idxl, int idxm, int idxi) {
      std::size_t tmpidx = this->_elemlen * idxelem;
      tmpidx += this->_linlen * idxi;
      tmpidx += idxl * idxl + idxl + idxm - this->_lMin * this->_lMin;
      return tmpidx;
   };

   // size of a matrix
   std::size_t matrixdimension(int elemnum) const {
      return this->_elemlen * elemnum + this->_linlen;
   };

 private:
   int _lMax, _lMin, _linlen, _npoly, _elemlen;
};

class LaplacianWeakForm {
 public:
   // constructors
   LaplacianWeakForm() {};

   LaplacianWeakForm(const GaussQuad::Quadrature1D<double> &q)
       : _q{q}, quadrature_declared{true} {
      auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
                                                    _q.Points().end());
      mat_gauss_deriv.reserve(_q.N());
      for (int idxk = 0; idxk < _q.N(); ++idxk) {
         std::vector<std::vector<double>> mat_tmp(_q.N(),
                                                  std::vector<double>(_q.N()));
         for (int idxi = 0; idxi < _q.N(); ++idxi) {

            std::vector<double> vec_tmp(_q.N());
            for (int idxj = 0; idxj < _q.N(); ++idxj) {
               vec_tmp[idxj] = pleg.Derivative(idxi, q.X(idxk)) *
                               pleg.Derivative(idxj, q.X(idxk));
            }
            mat_tmp[idxi] = vec_tmp;
         }
         mat_gauss_deriv.push_back(mat_tmp);
      }
   };

   LaplacianWeakForm(int l, double r0, double r1,
                     const GaussQuad::Quadrature1D<double> &q)
       : _l{l}, _r0{r0}, _r1{r1}, _q{q}, quadrature_declared{true},
         fully_declared{true}, l_avail{true}, interval_avail{true} {
      auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
                                                    _q.Points().end());
      mat_gauss_deriv.reserve(_q.N());
      for (int idxk = 0; idxk < _q.N(); ++idxk) {
         std::vector<std::vector<double>> mat_tmp(_q.N(),
                                                  std::vector<double>(_q.N()));
         for (int idxi = 0; idxi < _q.N(); ++idxi) {

            std::vector<double> vec_tmp(_q.N());
            for (int idxj = 0; idxj < _q.N(); ++idxj) {
               vec_tmp[idxj] = pleg.Derivative(idxi, q.X(idxk)) *
                               pleg.Derivative(idxj, q.X(idxk));
            }
            mat_tmp[idxi] = vec_tmp;
         }
         mat_gauss_deriv.push_back(mat_tmp);
      }
   };

   // appending information on
   void append_l(int l) {
      _l = l;
      l_avail = true;
      if (interval_avail) {
         fully_declared = true;
      };
      Omega_l0_two_square =
          static_cast<double>(l) * static_cast<double>((l + 1));
   };
   void append_interval(double r0, double r1) {
      _r0 = r0;
      _r1 = r1;
      invnodespacetwo = 2.0 / (r1 - r0);
      _aval = 0.5 * (r1 - r0);
      _bval = 0.5 * (r1 + r0);
      interval_avail = true;
      if (l_avail) {
         fully_declared = true;
      };
   };

   void declare_ortho() { is_ortho = true; }

   // outputting info about l, m, interval
   int l() {
      assert(l_avail && "l not provided");
      return _l;
   };
   //    int m() { return _m; };
   double r0() {
      assert(interval_avail && "interval not provided");
      return _r0;
   };
   double r1() {
      assert(interval_avail && "interval not provided");
      return _r1;
   };
   double weight(int n) { return _q.W(n); };

   // output matrix element
   double matrixelement(int i, int j) {
      assert(fully_declared && "Not fully declared");
      double tmp = 0.0;
      for (int idxk = 0; idxk < _q.N(); ++idxk) {
         double tmp1 = _aval * _q.X(idxk) + _bval;
         tmp += tmp1 * tmp1 * mat_gauss_deriv[idxk][i][j] * _q.W(idxk);
      }
      tmp *= invnodespacetwo;

      if (i == j) {
         tmp += Omega_l0_two_square * _q.W(i) * _aval;
      }

      return tmp;
   };

 private:
   bool quadrature_declared{false}, fully_declared{false}, l_avail{false},
       interval_avail{false}, is_ortho{false};
   int _l;
   double _r0, _r1, Omega_l0_two_square, invnodespacetwo, _aval, _bval;
   //    std::vector<double> _vec_nodes, _vec_a, _vec_b, _vec_insp;
   std::vector<std::vector<std::vector<double>>> mat_gauss_deriv;
   GaussQuad::Quadrature1D<double> _q;
};

class MatrixWeakForm {
 public:
   MatrixWeakForm() {};
   // MatrixWeakForm(const std::vector<double> &,
   //                const GaussQuad::Quadrature1D<double> &);
   MatrixWeakForm(const Radial_Tools::RadialMesh &,
                  const GaussQuad::Quadrature1D<double> &);

   // optional arguments:
   void add_dtn_map() { dirtoneum = true; };
   void no_dtn_map() { dirtoneum = false; };   // no DTN mapping
   void set_full_matrix() { leftlowertriangle = false; };
   void set_left_lower() { leftlowertriangle = true; };

   // output matrix element
   double matrixelement(int, int, int, int);

   // function to return the full spectral element matrix
   template <typename MSCALAR>
   Eigen::SparseMatrix<MSCALAR> fullmatrix(int, int);

 private:
   double outerradius;
   GaussQuad::Quadrature1D<double> _q;
   int _numlayers, _polyord;
   // std::vector<double> _vec_insp;
   std::vector<std::vector<std::vector<double>>> mat_gauss_deriv;
   Radial_Tools::RadialMesh _radial_mesh;
   bool dirtoneum{true}, leftlowertriangle{false};
};

// MatrixWeakForm::MatrixWeakForm(const std::vector<double> &vec_noderadii,
//                                const GaussQuad::Quadrature1D<double> &q)
//     : outerradius{vec_noderadii.back()} {
//    if (vec_noderadii.size() < 2) {
//       assert("Need to append a vector of nodes.");
//    };
//    _numlayers = vec_noderadii.size() - 1;
//    _polyord = q.N() - 1;
//    _q = q;

//    // finding values of derivatives at polynomial nodes
//    auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
//                                                  _q.Points().end());
//    mat_gauss_deriv.reserve(_q.N());
//    for (int idxk = 0; idxk < _q.N(); ++idxk) {
//       std::vector<std::vector<double>> mat_tmp(_q.N(),
//                                                std::vector<double>(_q.N()));
//       for (int idxi = 0; idxi < _q.N(); ++idxi) {

//          std::vector<double> vec_tmp(_q.N());
//          for (int idxj = 0; idxj < _q.N(); ++idxj) {
//             vec_tmp[idxj] = pleg.Derivative(idxi, q.X(idxk)) *
//                             pleg.Derivative(idxj, q.X(idxk));
//          }
//          mat_tmp[idxi] = vec_tmp;
//       }
//       mat_gauss_deriv.push_back(mat_tmp);
//    };

//    // initialise
//    _vec_a = std::vector<double>(vec_noderadii.size() - 1, 0.0);
//    _vec_b = _vec_a;
//    _vec_insp = _vec_a;

//    // fillout
//    for (int idx1 = 0; idx1 < vec_noderadii.size() - 1; ++idx1) {
//       _vec_a[idx1] = 0.5 * (vec_noderadii[idx1 + 1] - vec_noderadii[idx1]);
//       _vec_b[idx1] = 0.5 * (vec_noderadii[idx1 + 1] + vec_noderadii[idx1]);
//       _vec_insp[idx1] = 1.0 / _vec_a[idx1];
//    };
// };

MatrixWeakForm::MatrixWeakForm(const Radial_Tools::RadialMesh &radial_mesh,
                               const GaussQuad::Quadrature1D<double> &q)
    : outerradius{radial_mesh.OuterRadius()}, _radial_mesh{radial_mesh} {

   // fill out private variables
   if (_radial_mesh.NumberOfElements() < 2) {
      assert("Need to have a model with more than one element.");
   };
   _numlayers = _radial_mesh.NumberOfElements();
   _polyord = q.N() - 1;
   _q = q;

   // finding values of derivatives at polynomial nodes
   auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
                                                 _q.Points().end());
   mat_gauss_deriv.reserve(_q.N());
   for (int idxk = 0; idxk < _q.N(); ++idxk) {
      std::vector<std::vector<double>> mat_tmp(_q.N(),
                                               std::vector<double>(_q.N()));
      for (int idxi = 0; idxi < _q.N(); ++idxi) {

         std::vector<double> vec_tmp(_q.N());
         for (int idxj = 0; idxj < _q.N(); ++idxj) {
            vec_tmp[idxj] = pleg.Derivative(idxi, q.X(idxk)) *
                            pleg.Derivative(idxj, q.X(idxk));
         }
         mat_tmp[idxi] = vec_tmp;
      }
      mat_gauss_deriv.push_back(mat_tmp);
   };
};

double
MatrixWeakForm::matrixelement(int laynum, int l, int i, int j) {

   assert((laynum > -1 && laynum < _numlayers + 1) && "Not a layer");

   // main integral, ie r^2 f_1 f_2
   double tmp = 0.0;
   for (int idxk = 0; idxk < _polyord + 1; ++idxk) {
      double tmp1 = _radial_mesh.NodeRadius(laynum, idxk);
      tmp += tmp1 * tmp1 * mat_gauss_deriv[idxk][i][j] * _q.W(idxk);
   }
   tmp *= 2.0 / _radial_mesh.ElementWidth(laynum);

   // additional diagonal term
   if (i == j) {
      double Omega_l0_two_square =
          static_cast<double>(l) * static_cast<double>(l + 1);
      tmp += Omega_l0_two_square * _q.W(i) * 0.5 *
             _radial_mesh.ElementWidth(laynum);
   }

   return tmp;
};

template <typename MSCALAR>
Eigen::SparseMatrix<MSCALAR>
MatrixWeakForm::fullmatrix(int lmin, int lmax) {

   // declare indexing class instance
   MatrixIndices myindex(lmin, lmax, this->_polyord);

   Eigen::SparseMatrix<MSCALAR> mat_return;
   mat_return.resize(myindex.matrixdimension(this->_numlayers),
                     myindex.matrixdimension(this->_numlayers));

   // tripletlist vector
   using T = Eigen::Triplet<MSCALAR>;
   std::vector<T> tripletlist;
   tripletlist.reserve(
       (_numlayers * (((_polyord + 1) * (_polyord + 2)) / 2 - 1) + 1) *
       std::pow(lmax + 1, 2));

   // filling out matrix
   for (int idxelem = 0; idxelem < this->_numlayers; ++idxelem) {
      // over elements
      // loop of polynomial orders
      for (int idxi = 0; idxi < this->_polyord + 1; ++idxi) {
         // over internal nodes within each element
         int jmax;
         if (leftlowertriangle) {
            jmax = idxi + 1;
         } else {
            jmax = _polyord + 1;
         };
         for (int idxj = 0; idxj < jmax; ++idxj) {
            // over internal nodes within each element
            // base indices
            std::size_t xidx = myindex.fullindex(idxelem, lmin, -lmin, idxi);
            std::size_t yidx = myindex.fullindex(idxelem, lmin, -lmin, idxj);

            // loop over l's and m's
            for (int idxl = lmin; idxl < lmax + 1; ++idxl) {   // over l
               auto tmp = this->matrixelement(idxelem, idxl, idxi, idxj);

               // tmp doesn't change between different m's for the same l
               for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {   // over m
                  tripletlist.push_back(T(xidx, yidx, tmp));
                  ++xidx;
                  ++yidx;
               }
            }
         }
      }
   }

   // DTN map
   if (dirtoneum) {

      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
         std::size_t xidx =
             myindex.fullindex(_numlayers - 1, idxl, -idxl, _polyord);
         auto tmp = static_cast<MSCALAR>(idxl + 1) * outerradius;

         // tmp doesn't change between different m's for the same l
         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
            tripletlist.push_back(T(xidx, xidx, tmp));
            ++xidx;
         }
      }
   }
   mat_return.setFromTriplets(tripletlist.begin(), tripletlist.end());
   return mat_return;
};
}   // namespace SpectralElementTools
#endif