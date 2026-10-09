#ifndef PSEUDOSPECTRAL_OPERATOR3D_GUARD_H
#define PSEUDOSPECTRAL_OPERATOR3D_GUARD_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <GSHTrans/All>

#include <complex>
#include <cstddef>

#include "GeneralModels/Earth_Density_Models_3D.h"
#include "GeneralModels/Perturbation_Models.h"
#include "Operator_Actions.h"

// Eigen adapter for the mapped operator action.
template <typename MRScalar> class MatrixReplacement3D;

namespace Eigen::internal {
template <typename MRScalar>
struct traits<MatrixReplacement3D<MRScalar>>
    : public traits<Eigen::SparseMatrix<MRScalar>> {};
}  // namespace Eigen::internal

template <typename MRScalar>
class MatrixReplacement3D
    : public Eigen::EigenBase<MatrixReplacement3D<MRScalar>> {
 public:
   using Scalar = MRScalar;
   using RealScalar = double;
   using StorageIndex = int;
   using Index = Eigen::EigenBase<MatrixReplacement3D<MRScalar>>::Index;
   enum { ColsAtCompileTime = Eigen::Dynamic,
          MaxColsAtCompileTime = Eigen::Dynamic,
          IsRowMajor = false };

   // The model and its borrowed tensor storage must outlive this adapter.
   explicit MatrixReplacement3D(GeneralEarthModels::Density3D &model,
                                bool include_boundary = true)
       : _model(&model), _tensor_field(&model.LaplaceTensorRef()),
         _include_boundary(include_boundary) {}

   MatrixReplacement3D(
       GeneralEarthModels::Density3D &model,
       const GeneralEarthModels::MappingPerturbation &perturbation,
       bool include_boundary = true)
       : _model(&model), _tensor_field(&perturbation.ref_da()),
         _include_boundary(include_boundary) {}

   MatrixReplacement3D(GeneralEarthModels::Density3D &,
                       GeneralEarthModels::MappingPerturbation &&,
                       bool = true) = delete;
   MatrixReplacement3D(GeneralEarthModels::Density3D &,
                       const GeneralEarthModels::MappingPerturbation &&,
                       bool = true) = delete;

   Index rows() const { return Dimension(); }
   Index cols() const { return Dimension(); }

   template <typename Rhs>
   Eigen::Product<MatrixReplacement3D<Scalar>, Rhs, Eigen::AliasFreeProduct>
   operator*(const Eigen::MatrixBase<Rhs> &x) const {
      return Eigen::Product<MatrixReplacement3D<Scalar>, Rhs,
                            Eigen::AliasFreeProduct>(*this, x.derived());
   }

   template <typename Derived>
   Eigen::VectorXcd Apply(const Eigen::MatrixBase<Derived> &input) const {
      return GPLSpec::OperatorActions::ApplyMapped3D(
          *_model, *_tensor_field, input, _include_boundary);
   }

   MatrixReplacement3D WithBoundary(bool include_boundary) const {
      return MatrixReplacement3D(*_model, *_tensor_field, include_boundary);
   }

 private:
   MatrixReplacement3D(
       GeneralEarthModels::Density3D &model,
       const GPLSpec::OperatorActions::TensorField3D &tensor_field,
       bool include_boundary)
       : _model(&model), _tensor_field(&tensor_field),
         _include_boundary(include_boundary) {}

   Index Dimension() const {
      const auto &mesh = _model->Node_InformationP();
      const auto harmonic_count = static_cast<Index>(
          _model->GSH_Grid().MaxDegree() + 1) *
          (_model->GSH_Grid().MaxDegree() + 1);
      return mesh.NumberOfGlobalRadialNodes() * harmonic_count;
   }

   GeneralEarthModels::Density3D *_model;
   const GPLSpec::OperatorActions::TensorField3D *_tensor_field;
   bool _include_boundary;
};

namespace Eigen::internal {
template <typename MRScalar, typename Rhs>
struct generic_product_impl<MatrixReplacement3D<MRScalar>, Rhs, SparseShape,
                            DenseShape, GemvProduct>
    : generic_product_impl_base<
          MatrixReplacement3D<MRScalar>, Rhs,
          generic_product_impl<MatrixReplacement3D<MRScalar>, Rhs>> {
   using ProductScalar =
       typename Product<MatrixReplacement3D<MRScalar>, Rhs>::Scalar;

   template <typename Dest>
   static void scaleAndAddTo(Dest &dst,
                             const MatrixReplacement3D<MRScalar> &lhs,
                             const Rhs &rhs, const ProductScalar &alpha) {
      dst.noalias() += alpha * lhs.Apply(rhs);
   }
};
}  // namespace Eigen::internal

#endif
