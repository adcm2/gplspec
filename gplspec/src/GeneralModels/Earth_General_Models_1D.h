#ifndef GENERAL_EARTH_MODEL_1D_H
#define GENERAL_EARTH_MODEL_1D_H

#include "../Radial_Tools.h"
#include "../Spectral_Element_Tools.h"
#include "../testtools.h"
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>

// This namespace includes several different types of models. There are multiple
// constructors in general and this file contains the information for the
// spherically symmetric 1D model
namespace GeneralEarthModels {

// class spherical_test : public Radial_Tools::node_info
class spherical_1D {
 public:
   spherical_1D() {};   // default

   // with 1D input model that satisfies requirements for 1D model
   template <class model>
      requires PlanetaryModel::SphericalElasticModel<model>
   spherical_1D(const model &, const GaussQuad::Quadrature1D<double> &, double,
                double);

   // with path to file
   spherical_1D(const std::string &, const GaussQuad::Quadrature1D<double> &,
                double, double);

   /////////////////////////////////////////////////////
   /////////////////////////////////////////////////////
   //////////////////    outputs    ////////////////////
   /////////////////////////////////////////////////////
   /////////////////////////////////////////////////////

   // info on nodes
   auto Num_Elements() const { return _num_layers; };
   auto Poly_Order() const { return _poly_ord; };
   auto Vector_Nodes() const { return _vec_nodes; };
   //    auto All_Nodes() const { return _vec_allnodes; };
   auto Node_Information() { return node_data; };

   // interpolation polynomial derivative
   auto GaussDerivative() const { return _mat_gaussderiv; };
   auto q() const { return _q; };
   auto SpectralElementInformation() { return _spectral_info; };

   // lower triangle
   void SetLowerTriangle() { _spectral_info.set_left_lower(); };
   //    void SetLowerTriangle() { _spectral_info.set_left_lower(); };

   // norms
   double const LengthNorm() const { return length_norm; };
   double const MassNorm() const { return mass_norm; };
   double const TimeNorm() const { return time_norm; }
   double const DensityNorm() const { return density_norm; };
   double const InertiaNorm() const { return inertia_norm; };
   double const VelocityNorm() const { return velocity_norm; };
   double const AccelerationNorm() const { return acceleration_norm; };
   double const ForceNorm() const { return force_norm; };
   double const StressNorm() const { return stress_norm; };
   double const GravitationalConstant() const {
      return gravitational_constant;
   };
   double const GravitationalConstanNorm() const {
      return gravitational_constant_norm;
   };
   double const PotentialNorm() const { return potential_norm; };

   // model information
   // these all have the form (int idxelem, int idxpoly). They output the model
   // value within the idxelem element at the idxpoly node
   auto isSolid(int, int);
   auto Density(int, int);
   auto VP(int, int);
   auto VPV(int, int);
   auto VPH(int, int);
   auto VS(int, int);
   auto VSV(int, int);
   auto VSH(int, int);
   auto Eta(int, int);
   auto A(int, int);
   auto C(int, int);
   auto N(int, int);
   auto L(int, int);
   auto F(int, int);
   auto Kappa(int, int);
   auto Mu(int, int);

 private:
   int _num_layers, _poly_ord;
   double _max_radial_step;
   double length_norm, mass_norm, time_norm, density_norm, inertia_norm,
       velocity_norm, acceleration_norm, force_norm, stress_norm,
       gravitational_constant, gravitational_constant_norm, potential_norm;
   bool _isisotropic;
   std::vector<std::vector<bool>> _vec_issolid;
   std::vector<int> _vec_discontinuity_indices, _vec_laynum;
   std::vector<double> _vec_nodes, _vec_allnodes;
   std::vector<std::vector<double>> _vec_density, _vec_VP, _vec_VPV, _vec_VPH,
       _vec_VS, _vec_VSV, _vec_VSH, _vec_Eta, _vec_A, _vec_C, _vec_N, _vec_L,
       _vec_F, _vec_Kappa, _vec_Mu;
   GaussQuad::Quadrature1D<double> _q;
   Radial_Tools::node_info node_data;
   Radial_Tools::RadialMesh _radial_mesh;
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _mat_gaussderiv;
   SpectralElementTools::MatrixWeakForm _spectral_info;
};

template <class model>
   requires PlanetaryModel::SphericalElasticModel<model>
spherical_1D::spherical_1D(const model &inp_model,
                           const GaussQuad::Quadrature1D<double> &q,
                           double max_radial_step, double maxrad)
    : _q{q}, length_norm(inp_model.LengthNorm()),
      mass_norm(inp_model.MassNorm()), time_norm(inp_model.TimeNorm()),
      density_norm(inp_model.DensityNorm()),
      velocity_norm(inp_model.VelocityNorm()),
      acceleration_norm(inp_model.AccelerationNorm()),
      force_norm(inp_model.ForceNorm()), stress_norm(inp_model.StressNorm()),
      inertia_norm(inp_model.InertiaNorm()),
      gravitational_constant_norm(inp_model.GravitationalConstant()),
      gravitational_constant(6.6743015 * std::pow(10.0, -11.0) /
                             inp_model.GravitationalConstant()),
      potential_norm(
          std::pow(inp_model.LengthNorm() / inp_model.TimeNorm(), 2.0)) {

   // polynomial order
   _poly_ord = _q.N() - 1;

   // find nodes and number of layers in final stored model
   _vec_nodes = Radial_Tools::Radial_Node(inp_model, max_radial_step, maxrad);
   _vec_allnodes = Radial_Tools::All_Node(_vec_nodes, q);
   _num_layers = _vec_nodes.size() - 1;

   // constructing the class instance that stores all the information about
   // the nodes
   node_data = Radial_Tools::node_info(inp_model, q, max_radial_step, maxrad);
   _radial_mesh =
       Radial_Tools::RadialMesh(inp_model, q, max_radial_step, maxrad);

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

   // find all model information
   _vec_density.resize(_num_layers);
   _vec_VP.resize(_num_layers);
   _vec_VPV.resize(_num_layers);
   _vec_VPH.resize(_num_layers);
   _vec_VS.resize(_num_layers);
   _vec_VSV.resize(_num_layers);
   _vec_VSH.resize(_num_layers);
   _vec_Eta.resize(_num_layers);
   _vec_A.resize(_num_layers);
   _vec_C.resize(_num_layers);
   _vec_N.resize(_num_layers);
   _vec_L.resize(_num_layers);
   _vec_F.resize(_num_layers);
   _vec_Kappa.resize(_num_layers);
   _vec_Mu.resize(_num_layers);
   _vec_issolid.resize(_num_layers);
   //   _vec_isisotropic.resize(_num_layers);

   int laynum = 0;   // layer #
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      // temporary vectors
      std::vector<double> tmp_density(_poly_ord + 1, 0.0),
          tmp_VP(_poly_ord + 1, 0.0), tmp_VPV(_poly_ord + 1, 0.0),
          tmp_VPH(_poly_ord + 1, 0.0), tmp_VS(_poly_ord + 1, 0.0),
          tmp_VSV(_poly_ord + 1, 0.0), tmp_VSH(_poly_ord + 1, 0.0),
          tmp_Eta(_poly_ord + 1, 0.0), tmp_A(_poly_ord + 1, 0.0),
          tmp_C(_poly_ord + 1, 0.0), tmp_N(_poly_ord + 1, 0.0),
          tmp_L(_poly_ord + 1, 0.0), tmp_F(_poly_ord + 1, 0.0),
          tmp_Kappa(_poly_ord + 1, 0.0), tmp_Mu(_poly_ord + 1, 0.0);
      std::vector<bool> tmp_issolid(_poly_ord + 1, 0.0);   // tmp vector
      std::size_t idx_global = idxelem * _poly_ord;        // global idx

      // loop through polynomial indices
      for (int idxpoly = 0; idxpoly < _poly_ord + 1; ++idxpoly) {
         double rad_current = _vec_allnodes[idx_global];
         if (laynum < inp_model.NumberOfLayers()) {
            tmp_density[idxpoly] = inp_model.Density(laynum)(rad_current);
            tmp_VP[idxpoly] = inp_model.VP(laynum)(rad_current);
            tmp_VPV[idxpoly] = inp_model.VPV(laynum)(rad_current);
            tmp_VPH[idxpoly] = inp_model.VPH(laynum)(rad_current);
            tmp_VS[idxpoly] = inp_model.VS(laynum)(rad_current);
            tmp_VSV[idxpoly] = inp_model.VSV(laynum)(rad_current);
            tmp_VSH[idxpoly] = inp_model.VSH(laynum)(rad_current);
            tmp_Eta[idxpoly] = inp_model.Eta(laynum)(rad_current);
            tmp_A[idxpoly] = inp_model.A(laynum)(rad_current);
            tmp_C[idxpoly] = inp_model.C(laynum)(rad_current);
            tmp_N[idxpoly] = inp_model.N(laynum)(rad_current);
            tmp_L[idxpoly] = inp_model.L(laynum)(rad_current);
            tmp_F[idxpoly] = inp_model.F(laynum)(rad_current);
            tmp_Kappa[idxpoly] = inp_model.Kappa(laynum)(rad_current);
            tmp_Mu[idxpoly] = inp_model.Mu(laynum)(rad_current);
            tmp_issolid[idxpoly] = inp_model.IsSolid(laynum);
         }
         ++idx_global;
      }

      // increment layer
      if (idxelem + 1 == node_data.IdxUpperRadius(laynum)) {
         ++laynum;
      };

      // push back
      _vec_density[idxelem] = tmp_density;
      _vec_VP[idxelem] = tmp_VP;
      _vec_VPV[idxelem] = tmp_VPV;
      _vec_VPH[idxelem] = tmp_VPH;
      _vec_VS[idxelem] = tmp_VS;
      _vec_VSV[idxelem] = tmp_VSV;
      _vec_VSH[idxelem] = tmp_VSH;
      _vec_Eta[idxelem] = tmp_Eta;
      _vec_A[idxelem] = tmp_A;
      _vec_C[idxelem] = tmp_C;
      _vec_N[idxelem] = tmp_N;
      _vec_L[idxelem] = tmp_L;
      _vec_F[idxelem] = tmp_F;
      _vec_Kappa[idxelem] = tmp_Kappa;
      _vec_Mu[idxelem] = tmp_Mu;
      _vec_issolid[idxelem] = tmp_issolid;
   }

   // finding the spectral element grid
   // _spectral_info =
   //     SpectralElementTools::MatrixWeakForm(node_data.ElementNodes(), q);
   _spectral_info = SpectralElementTools::MatrixWeakForm(_radial_mesh, q);
};

// with path to file
spherical_1D::spherical_1D(const std::string &pathtofile,
                           const GaussQuad::Quadrature1D<double> &q,
                           double max_radial_step, double maxrad)
    : spherical_1D::spherical_1D(TestTools::EarthModel(pathtofile), q,
                                 max_radial_step, maxrad) {};

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
auto
spherical_1D::isSolid(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_issolid[idxelem][idxpoly];
};

auto
spherical_1D::Density(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_density[idxelem][idxpoly];
};
auto
spherical_1D::VP(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_VP[idxelem][idxpoly];
};
auto
spherical_1D::VPV(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_VPV[idxelem][idxpoly];
};
auto
spherical_1D::VPH(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_VPH[idxelem][idxpoly];
};
auto
spherical_1D::VS(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_VS[idxelem][idxpoly];
};
auto
spherical_1D::VSV(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_VSV[idxelem][idxpoly];
};
auto
spherical_1D::VSH(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_VSH[idxelem][idxpoly];
};
auto
spherical_1D::Eta(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_Eta[idxelem][idxpoly];
};
auto
spherical_1D::A(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_A[idxelem][idxpoly];
};
auto
spherical_1D::C(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_C[idxelem][idxpoly];
};
auto
spherical_1D::N(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_N[idxelem][idxpoly];
};
auto
spherical_1D::L(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_L[idxelem][idxpoly];
};
auto
spherical_1D::F(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_F[idxelem][idxpoly];
};
auto
spherical_1D::Kappa(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_Kappa[idxelem][idxpoly];
};
auto
spherical_1D::Mu(int idxelem, int idxpoly) {
   assert(idxelem < _num_layers && "input layer number outside model");
   assert(idxpoly < _poly_ord + 1 && "input polynomial number too large");
   return _vec_Mu[idxelem][idxpoly];
};
};   // namespace GeneralEarthModels
#endif