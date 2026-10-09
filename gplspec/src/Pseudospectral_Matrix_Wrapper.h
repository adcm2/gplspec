#ifndef PSEUDOSPECTRAL_OPERATOR_GUARD_H
#define PSEUDOSPECTRAL_OPERATOR_GUARD_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <GSHTrans/All>

#include <complex>
#include <cstddef>

#include "GeneralModels/Earth_General_Models_1D.h"
#include "Operator_Actions.h"

// Eigen adapter for the spherical operator action. The mathematical action and
// its boundary choice live outside Eigen's product specialization.
template <typename MRScalar> class MatrixReplacement;

namespace Eigen::internal {
template <typename MRScalar>
struct traits<MatrixReplacement<MRScalar>>
    : public traits<Eigen::SparseMatrix<MRScalar>> {};
}  // namespace Eigen::internal

template <typename MRScalar>
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<MRScalar>> {
 public:
   using Scalar = MRScalar;
   using RealScalar = double;
   using StorageIndex = int;
   using Index = Eigen::EigenBase<MatrixReplacement<MRScalar>>::Index;
   using Grid = GSHTrans::GaussLegendreGrid<double, GSHTrans::All, GSHTrans::All>;
   enum { ColsAtCompileTime = Eigen::Dynamic,
          MaxColsAtCompileTime = Eigen::Dynamic,
          IsRowMajor = false };

   // The model and grid are borrowed and must outlive this adapter.
   MatrixReplacement(GeneralEarthModels::spherical_1D &model,
                     const Grid &grid, bool include_boundary = true)
       : _model(&model), _grid(&grid), _include_boundary(include_boundary) {}
   MatrixReplacement(GeneralEarthModels::spherical_1D &, Grid &&,
                     bool = true) = delete;
   MatrixReplacement(GeneralEarthModels::spherical_1D &, const Grid &&,
                     bool = true) = delete;

   Index rows() const { return Dimension(); }
   Index cols() const { return Dimension(); }

   template <typename Rhs>
   Eigen::Product<MatrixReplacement<Scalar>, Rhs, Eigen::AliasFreeProduct>
   operator*(const Eigen::MatrixBase<Rhs> &x) const {
      return Eigen::Product<MatrixReplacement<Scalar>, Rhs,
                            Eigen::AliasFreeProduct>(*this, x.derived());
   }

   template <typename Derived>
   Eigen::VectorXcd Apply(const Eigen::MatrixBase<Derived> &input) const {
      return GPLSpec::OperatorActions::ApplySpherical1D(
          *_model, *_grid, input, _include_boundary);
   }

   MatrixReplacement WithBoundary(bool include_boundary) const {
      return MatrixReplacement(*_model, *_grid, include_boundary);
   }

 private:
   Index Dimension() const {
      const auto &mesh = _model->Node_Information();
      const auto harmonic_count = static_cast<Index>(_grid->MaxDegree() + 1) *
                                  (_grid->MaxDegree() + 1);
      return mesh.NumberOfGlobalRadialNodes() * harmonic_count;
   }

   GeneralEarthModels::spherical_1D *_model;
   const Grid *_grid;
   bool _include_boundary;
};

namespace Eigen::internal {
template <typename MRScalar, typename Rhs>
struct generic_product_impl<MatrixReplacement<MRScalar>, Rhs, SparseShape,
                            DenseShape, GemvProduct>
    : generic_product_impl_base<
          MatrixReplacement<MRScalar>, Rhs,
          generic_product_impl<MatrixReplacement<MRScalar>, Rhs>> {
   using ProductScalar =
       typename Product<MatrixReplacement<MRScalar>, Rhs>::Scalar;

   template <typename Dest>
   static void scaleAndAddTo(Dest &dst,
                             const MatrixReplacement<MRScalar> &lhs,
                             const Rhs &rhs, const ProductScalar &alpha) {
      dst.noalias() += alpha * lhs.Apply(rhs);
   }
};
}  // namespace Eigen::internal

#endif
