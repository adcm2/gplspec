#ifndef GENERAL_EARTH_MODEL_3D_H
#define GENERAL_EARTH_MODEL_3D_H

#include "../Radial_Tools.h"
#include "../Spectral_Element_Tools.h"
#include "../testtools.h"
#include <GSHTrans/All>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>

// This namespace includes several different types of models. There are multiple
// constructors in general and this file has the class for the most general 3D
// model
namespace GeneralEarthModels {
class General3D {
 public:
   //    using namespace GSHTrans;
   using MRange = GSHTrans::All;
   using NRange = GSHTrans::All;
   using Grid = GSHTrans::GaussLegendreGrid<double, MRange, NRange>;

   General3D() {};

   // with 1D input model that satisfies requirements for 1D model
   template <class model>
      requires PlanetaryModel::SphericalElasticModel<model>
   General3D(const model &, const GaussQuad::Quadrature1D<double> &,
             const Grid &, double, double);

   /////////////////////////////////////////////////////
   /////////////////////////////////////////////////////
   //////////////////    outputs    ////////////////////
   /////////////////////////////////////////////////////
   /////////////////////////////////////////////////////

   // info on nodes
   auto Num_Elements() const { return _num_layers; };
   auto Poly_Order() const { return _poly_ord; };
   auto Node_Information() { return node_data; };

   // interpolation polynomial derivative
   auto GaussDerivative() const { return _mat_gaussderiv; };
   auto q() const { return _q; };
   auto SpectralElementInformation() { return _spectral_info; };
   auto GSH_Grid() const { return _grid; };

   // density
   auto Density() const { return _vec_density; };
   auto Density_Radius(int, int) const;
   auto Density_Point(int, int, int) const;

   // VP
   auto VP() const { return _vec_VP; };
   auto VP_Radius(int, int) const;
   auto VP_Point(int, int, int) const;

   // VPV
   auto VPV() const { return _vec_VPV; };
   auto VPV_Radius(int, int) const;
   auto VPV_Point(int, int, int) const;

   // VPH
   auto VPH() const { return _vec_VPH; };
   auto VPH_Radius(int, int) const;
   auto VPH_Point(int, int, int) const;

   // VS
   auto VS() const { return _vec_VS; };
   auto VS_Radius(int, int) const;
   auto VS_Point(int, int, int) const;

   // VSV
   auto VSV() const { return _vec_VSV; };
   auto VSV_Radius(int, int) const;
   auto VSV_Point(int, int, int) const;

   // VSH
   auto VSH() const { return _vec_VSH; };
   auto VSH_Radius(int, int) const;
   auto VSH_Point(int, int, int) const;

   // Eta
   auto Eta() const { return _vec_Eta; };
   auto Eta_Radius(int, int) const;
   auto Eta_Point(int, int, int) const;

   // A
   auto A() const { return _vec_A; };
   auto A_Radius(int, int) const;
   auto A_Point(int, int, int) const;

   // C
   auto C() const { return _vec_C; };
   auto C_Radius(int, int) const;
   auto C_Point(int, int, int) const;

   // N
   auto N() const { return _vec_N; };
   auto N_Radius(int, int) const;
   auto N_Point(int, int, int) const;

   // L
   auto L() const { return _vec_L; };
   auto L_Radius(int, int) const;
   auto L_Point(int, int, int) const;

   // F
   auto F() const { return _vec_F; };
   auto F_Radius(int, int) const;
   auto F_Point(int, int, int) const;

   // Kappa
   auto Kappa() const { return _vec_Kappa; };
   auto Kappa_Radius(int, int) const;
   auto Kappa_Point(int, int, int) const;
   // Mu
   auto Mu() const { return _vec_Mu; };
   auto Mu_Radius(int, int) const;
   auto Mu_Point(int, int, int) const;

   // h
   auto Mapping() const { return _vec_h; };
   auto Mapping_Radius(int, int) const;
   auto Mapping_Point(int, int, int) const;

   auto LaplaceTensor() const { return _vec_a; };
   auto LaplaceTensor_Radius(int, int) const;
   auto LaplaceTensor_Point(int, int, int) const;

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

 private:
   using vecdb = std::vector<double>;
   using vvecdb = std::vector<vecdb>;
   using vvvecdb = std::vector<vvecdb>;
   using veceig = std::vector<Eigen::Matrix3cd>;
   using vveceig = std::vector<veceig>;
   using vvveceig = std::vector<std::vector<std::vector<Eigen::Matrix3cd>>>;
   int _num_layers, _poly_ord;
   double _max_radial_step;
   double length_norm, mass_norm, time_norm, density_norm, inertia_norm,
       velocity_norm, acceleration_norm, force_norm, stress_norm,
       gravitational_constant;
   bool _isisotropic;
   std::vector<std::vector<bool>> _vec_issolid;
   std::vector<int> _vec_discontinuity_indices, _vec_laynum;
   //    vecdb _vec_nodes, _vec_allnodes;

   // model parameters
   vvvecdb _vec_h, _vec_density, _vec_VP, _vec_VPV, _vec_VPH, _vec_VS, _vec_VSV,
       _vec_VSH, _vec_Eta, _vec_A, _vec_C, _vec_N, _vec_L, _vec_F, _vec_Kappa,
       _vec_Mu;
   ;   //, _vec_h, _vec_a;
   vvveceig _vec_a;

   GaussQuad::Quadrature1D<double> _q;
   Radial_Tools::node_info node_data;
   Radial_Tools::RadialMesh _radial_mesh;
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _mat_gaussderiv;
   SpectralElementTools::MatrixWeakForm _spectral_info;

   // grid
   Grid _grid;
};

template <class model>
   requires PlanetaryModel::SphericalElasticModel<model>
General3D::General3D(const model &inp_model,
                     const GaussQuad::Quadrature1D<double> &q, const Grid &grid,
                     double max_radial_step, double maxrad)
    : _q{q}, _grid{grid}, length_norm(inp_model.LengthNorm()),
      mass_norm(inp_model.MassNorm()), time_norm(inp_model.TimeNorm()),
      density_norm(inp_model.DensityNorm()),
      velocity_norm(inp_model.VelocityNorm()),
      acceleration_norm(inp_model.AccelerationNorm()),
      force_norm(inp_model.ForceNorm()), stress_norm(inp_model.StressNorm()),
      inertia_norm(inp_model.InertiaNorm()),
      gravitational_constant(inp_model.GravitationalConstant()) {

   // polynomial order
   _poly_ord = _q.N() - 1;

   // find nodes and number of layers in final stored model
   //    _vec_nodes = Radial_Tools::Radial_Node(inp_model, max_radial_step,
   //    maxrad); _vec_allnodes = Radial_Tools::All_Node(_vec_nodes, q);
   // constructing the class instance that stores all the information about
   // the nodes
   node_data = Radial_Tools::node_info(inp_model, q, max_radial_step, maxrad);
   _radial_mesh =
       Radial_Tools::RadialMesh(inp_model, q, max_radial_step, maxrad);
   _num_layers = node_data.NumberOfElements();

   // default h
   auto spatialsize = _grid.Longitudes().size() * _grid.CoLatitudes().size();
   _vec_h = vvvecdb(_num_layers, vvecdb(_q.N(), vecdb(spatialsize, 0.0)));
   {
      Eigen::Matrix3cd mat_a0 = Eigen::Matrix3cd::Zero();
      _vec_a =
          vvveceig(_num_layers, vveceig(_q.N(), veceig(spatialsize, mat_a0)));
   }

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

   // finding the spectral element grid
   // _spectral_info =
   //     SpectralElementTools::MatrixWeakForm(node_data.ElementNodes(), q);
   _spectral_info = SpectralElementTools::MatrixWeakForm(_radial_mesh, q);

   // density
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
   for (int idxelem = 0; idxelem < _num_layers; ++idxelem) {

      // temporary vectors
      vvecdb tmp_density, tmp_VP, tmp_VPV, tmp_VPH, tmp_VS, tmp_VSV, tmp_VSH,
          tmp_Eta, tmp_A, tmp_C, tmp_N, tmp_L, tmp_F, tmp_Kappa,
          tmp_Mu;   // to contain all values within the element
      int laynum = node_data.LayerNumber(idxelem);

      // loop through polynomial indices
      for (int idxpoly = 0; idxpoly < _poly_ord + 1; ++idxpoly) {
         double rad_current =
             node_data.NodeRadius(idxelem, idxpoly);   // current radius
         double density_1D = 0.0, VP_1D = 0.0, VPV_1D = 0.0, VPH_1D = 0.0,
                VS_1D = 0.0, VSV_1D = 0.0, VSH_1D = 0.0, Eta_1D = 0.0,
                A_1D = 0.0, C_1D = 0.0, N_1D = 0.0, L_1D = 0.0, F_1D = 0.0,
                Kappa_1D = 0.0,
                Mu_1D = 0.0;   // initialise density
                               //  std::cout << "seg fault check #1\n";

         // set density iff inside input model
         if (laynum < inp_model.NumberOfLayers()) {
            density_1D = inp_model.Density(laynum)(rad_current);
            VP_1D = inp_model.VP(laynum)(rad_current);
            VPV_1D = inp_model.VPV(laynum)(rad_current);
            VPH_1D = inp_model.VPH(laynum)(rad_current);
            VS_1D = inp_model.VS(laynum)(rad_current);
            VSV_1D = inp_model.VSV(laynum)(rad_current);
            VSH_1D = inp_model.VSH(laynum)(rad_current);
            Eta_1D = inp_model.Eta(laynum)(rad_current);
            A_1D = inp_model.A(laynum)(rad_current);
            C_1D = inp_model.C(laynum)(rad_current);
            N_1D = inp_model.N(laynum)(rad_current);
            L_1D = inp_model.L(laynum)(rad_current);
            F_1D = inp_model.F(laynum)(rad_current);
            Kappa_1D = inp_model.Kappa(laynum)(rad_current);
            Mu_1D = inp_model.Mu(laynum)(rad_current);
         }
         //  std::cout << "seg fault check #2\n";

         // fill out all values at particular radius
         vecdb tmp_vec_level_density(spatialsize, 0.0),
             tmp_vec_level_VP(spatialsize, 0.0),
             tmp_vec_level_VPV(spatialsize, 0.0),
             tmp_vec_level_VPH(spatialsize, 0.0),
             tmp_vec_level_VS(spatialsize, 0.0),
             tmp_vec_level_VSV(spatialsize, 0.0),
             tmp_vec_level_VSH(spatialsize, 0.0),
             tmp_vec_level_Eta(spatialsize, 0.0),
             tmp_vec_level_A(spatialsize, 0.0),
             tmp_vec_level_C(spatialsize, 0.0),
             tmp_vec_level_N(spatialsize, 0.0),
             tmp_vec_level_L(spatialsize, 0.0),
             tmp_vec_level_F(spatialsize, 0.0),
             tmp_vec_level_Kappa(spatialsize, 0.0),
             tmp_vec_level_Mu(spatialsize, 0.0);
         //  std::cout << "seg fault check #3\n";
         for (int idxspatial = 0; idxspatial < spatialsize; ++idxspatial) {
            tmp_vec_level_density[idxspatial] = density_1D;
            tmp_vec_level_VP[idxspatial] = VP_1D;
            tmp_vec_level_VPV[idxspatial] = VPV_1D;
            tmp_vec_level_VPH[idxspatial] = VPH_1D;
            tmp_vec_level_VS[idxspatial] = VS_1D;
            tmp_vec_level_VSV[idxspatial] = VSV_1D;
            tmp_vec_level_VSH[idxspatial] = VSH_1D;
            tmp_vec_level_Eta[idxspatial] = Eta_1D;
            tmp_vec_level_A[idxspatial] = A_1D;
            tmp_vec_level_C[idxspatial] = C_1D;
            tmp_vec_level_N[idxspatial] = N_1D;
            tmp_vec_level_L[idxspatial] = L_1D;
            tmp_vec_level_F[idxspatial] = F_1D;
            tmp_vec_level_Kappa[idxspatial] = Kappa_1D;
            tmp_vec_level_Mu[idxspatial] = Mu_1D;
         }
         //  std::cout << "seg fault check #4\n";
         tmp_density.push_back(tmp_vec_level_density);
         tmp_VP.push_back(tmp_vec_level_VP);
         tmp_VPV.push_back(tmp_vec_level_VPV);
         tmp_VPH.push_back(tmp_vec_level_VPH);
         tmp_VS.push_back(tmp_vec_level_VS);
         tmp_VSV.push_back(tmp_vec_level_VSV);
         tmp_VSH.push_back(tmp_vec_level_VSH);
         tmp_Eta.push_back(tmp_vec_level_Eta);
         tmp_A.push_back(tmp_vec_level_A);
         tmp_C.push_back(tmp_vec_level_C);
         tmp_N.push_back(tmp_vec_level_N);
         tmp_L.push_back(tmp_vec_level_L);
         tmp_F.push_back(tmp_vec_level_F);
         tmp_Kappa.push_back(tmp_vec_level_Kappa);
         tmp_Mu.push_back(tmp_vec_level_Mu);
      }

      // assign to the idxelem element
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
   }
};

// density
auto
General3D::Density_Radius(int idxelem, int idxpoly) const {
   return _vec_density[idxelem][idxpoly];
};
auto
General3D::Density_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_density[idxelem][idxpoly][idxspatial];
};

// mapping
auto
General3D::Mapping_Radius(int idxelem, int idxpoly) const {
   return _vec_h[idxelem][idxpoly];
};
auto
General3D::Mapping_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_h[idxelem][idxpoly][idxspatial];
};

// return Laplace tensor a = J C^{-1}
auto
General3D::LaplaceTensor_Radius(int idxelem, int idxpoly) const {
   return _vec_a[idxelem][idxpoly];
};
auto
General3D::LaplaceTensor_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_a[idxelem][idxpoly][idxspatial];
};

// VP
auto
General3D::VP_Radius(int idxelem, int idxpoly) const {
   return _vec_VP[idxelem][idxpoly];
};
auto
General3D::VP_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_VP[idxelem][idxpoly][idxspatial];
};

// VPV
auto
General3D::VPV_Radius(int idxelem, int idxpoly) const {
   return _vec_VPV[idxelem][idxpoly];
};
auto
General3D::VPV_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_VPV[idxelem][idxpoly][idxspatial];
};

// VPH
auto
General3D::VPH_Radius(int idxelem, int idxpoly) const {
   return _vec_VPH[idxelem][idxpoly];
};
auto
General3D::VPH_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_VPH[idxelem][idxpoly][idxspatial];
};

// VS
auto
General3D::VS_Radius(int idxelem, int idxpoly) const {
   return _vec_VS[idxelem][idxpoly];
};
auto
General3D::VS_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_VS[idxelem][idxpoly][idxspatial];
};

// VSV
auto
General3D::VSV_Radius(int idxelem, int idxpoly) const {
   return _vec_VSV[idxelem][idxpoly];
};
auto
General3D::VSV_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_VSV[idxelem][idxpoly][idxspatial];
};

// VSH
auto
General3D::VSH_Radius(int idxelem, int idxpoly) const {
   return _vec_VSH[idxelem][idxpoly];
};
auto
General3D::VSH_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_VSH[idxelem][idxpoly][idxspatial];
};

// Eta
auto
General3D::Eta_Radius(int idxelem, int idxpoly) const {
   return _vec_Eta[idxelem][idxpoly];
};
auto
General3D::Eta_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_Eta[idxelem][idxpoly][idxspatial];
};

// A
auto
General3D::A_Radius(int idxelem, int idxpoly) const {
   return _vec_A[idxelem][idxpoly];
};
auto
General3D::A_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_A[idxelem][idxpoly][idxspatial];
};

// C
auto
General3D::C_Radius(int idxelem, int idxpoly) const {
   return _vec_C[idxelem][idxpoly];
};
auto
General3D::C_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_C[idxelem][idxpoly][idxspatial];
};

// N
auto
General3D::N_Radius(int idxelem, int idxpoly) const {
   return _vec_N[idxelem][idxpoly];
};
auto
General3D::N_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_N[idxelem][idxpoly][idxspatial];
};

// L
auto
General3D::L_Radius(int idxelem, int idxpoly) const {
   return _vec_C[idxelem][idxpoly];
};
auto
General3D::L_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_C[idxelem][idxpoly][idxspatial];
};

// F
auto
General3D::F_Radius(int idxelem, int idxpoly) const {
   return _vec_F[idxelem][idxpoly];
};
auto
General3D::F_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_F[idxelem][idxpoly][idxspatial];
};

// Kappa
auto
General3D::Kappa_Radius(int idxelem, int idxpoly) const {
   return _vec_Kappa[idxelem][idxpoly];
};
auto
General3D::Kappa_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_Kappa[idxelem][idxpoly][idxspatial];
};
// Mu
auto
General3D::Mu_Radius(int idxelem, int idxpoly) const {
   return _vec_Mu[idxelem][idxpoly];
};
auto
General3D::Mu_Point(int idxelem, int idxpoly, int idxspatial) const {
   return _vec_Mu[idxelem][idxpoly][idxspatial];
};
};   // namespace GeneralEarthModels
#endif