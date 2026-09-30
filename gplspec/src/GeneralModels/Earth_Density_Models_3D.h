#ifndef DENSITY_EARTH_MODEL_3D_H
#define DENSITY_EARTH_MODEL_3D_H

// #include "../SimpleModels"
#include "../Radial_Tools.h"
// #include "Simple_Models.h"
#include "../Spectral_Element_Tools.h"
#include "../testtools.h"
#include <GSHTrans/All>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>
#include <TomographyModels/All>

// This namespace includes several different types of models. There are multiple
// constructors in general and this file has the class for the most general 3D
// model
namespace GeneralEarthModels {
class Density3D {
 public:
   //    using namespace GSHTrans;
   using MRange = GSHTrans::All;
   using NRange = GSHTrans::All;
   using Grid = GSHTrans::GaussLegendreGrid<double, MRange, NRange>;

   Density3D() {};

   // template <class model>
   //    requires PlanetaryModel::BasicSphericalDensityModel<model, int, double>
   // Density3D(const model &, const int, const int, double, double);

   // with path to file for 1D model
   // Density3D(const std::string &, const int, const int, double, double);

   // with path to file for 1D and tomography models
   // Density3D(const std::string &, const std::string &, const int, const int,
   //           double, double);

   inline static Density3D SphericalHomogeneousPlanet(double, double, double, double,
                                               double, double, double);

   inline static Density3D SphericalLayeredPlanet(std::vector<double>,
                                           std::vector<double>, double, double,
                                           double, double, double);

   inline static Density3D OneDimensionalPlanetFromFile(const std::string &, const int,
                                                 const int, double, double);
   inline static Density3D SphericalThreeDimensionalPlanetFromFile(const std::string &,
                                                            const std::string &,
                                                            const int,
                                                            const int, double,
                                                            double);
   template <class model>
      requires PlanetaryModel::BasicSphericalDensityModel<model, int, double>
   static Density3D OneDimensionalPlanetFromModel(const model &, const int,
                                                  const int, double, double);
   template <class model, class tomomodel>
      requires PlanetaryModel::BasicSphericalDensityModel<model, int, double> &&
               PlanetaryModel::TomographyModel<tomomodel>
   static Density3D SphericalThreeDimensionalPlanetFromModel(
       const model &, const tomomodel &, const int, const int, double, double);

   template <class model, class tomomodel, class mapclass>
      requires PlanetaryModel::BasicSphericalDensityModel<model, int, double> &&
               PlanetaryModel::TomographyModel<tomomodel> &&
               PlanetaryModel::RadialMappingClass<mapclass>
   static Density3D
   PhysicalPlanetModelWithProvidedMapping(const model &, const tomomodel &,
                                          const mapclass &, const int,
                                          const int, double, double);

   // template <class model, class tomomodel>
   //       requires PlanetaryModel::BasicSphericalDensityModel<model, int,
   //       double> &&
   //                PlanetaryModel::TomographyModel<tomomodel>
   //    static Density3D
   //    AsphericalTest(const model &, const tomomodel &, const int,
   //                                    const int, double, double);

   // full constructor
   template <class model, class tomomodel, class mapclass>
      requires PlanetaryModel::BasicSphericalDensityModel<model, int, double> &&
               PlanetaryModel::TomographyModel<tomomodel> &&
               PlanetaryModel::RadialMappingClass<mapclass>
   Density3D(const model &, const tomomodel &, const mapclass &, const int,
             const int, double, double);

   // constructor for homogeneous physical planet with mapping/layers
   template <class mapclass>
      requires PlanetaryModel::RadialMappingClass<mapclass>
   Density3D(
       double physicalradius, double physicaldensity, double referentialradius,
       const mapclass &inp_map, const int npoly, const int lMax,
       double ilength_norm = EarthModels::EarthConstants<double>().LengthNorm(),
       double itime_norm = EarthModels::EarthConstants<double>().TimeNorm(),
       double imass_norm = EarthModels::EarthConstants<double>().MassNorm(),
       double maxradialstep = 0.01, double ballrad = 1.2);

   // constructor for homogeneous physical planet with mapping/layers
   // this is used for Phobos
   inline Density3D(
       double physicaldensity, std::string pathtofile, const int lMaxMod,
       const int npoly, const int lMax,
       double ilength_norm = EarthModels::EarthConstants<double>().LengthNorm(),
       double itime_norm = EarthModels::EarthConstants<double>().TimeNorm(),
       double imass_norm = EarthModels::EarthConstants<double>().MassNorm(),
       double maxradialstep = 0.01, double ballrad = 1.2,
       bool randdensity = false);

   inline Density3D(
       double physicaldensity, std::string pathtofile, const int npoly,
       const int lMax,
       double ilength_norm = EarthModels::EarthConstants<double>().LengthNorm(),
       double itime_norm = EarthModels::EarthConstants<double>().TimeNorm(),
       double imass_norm = EarthModels::EarthConstants<double>().MassNorm(),
       double maxradialstep = 0.01, double ballrad = 1.2,
       bool randdensity = false);

   /////////////////////////////////////////////////////
   /////////////////////////////////////////////////////
   //////////////////    outputs    ////////////////////
   /////////////////////////////////////////////////////
   /////////////////////////////////////////////////////

   // info on nodes
   inline auto Num_Elements() const;
   inline auto Poly_Order() const;
   inline auto Node_Information() const;
   inline const Radial_Tools::RadialMesh &Node_InformationP() const;

   // interpolation polynomial derivative
   inline auto GaussDerivative() const;
   inline auto GaussDerivative(int, int) const;
   inline auto q() const;
   const GaussQuad::Quadrature1D<double> &qP() const { return _q; };

   inline auto SpectralElementInformation();
   inline auto GSH_Grid() const;
   inline const Grid &GSH_GridP() const;

   // density
   inline auto Density() const;
   inline auto DensityAtRadialNode(int, int) const;
   inline auto Density_Point(int, int, int) const;

   // h
   inline auto Mapping() const;
   inline auto MappingAtRadialNode(int, int) const;
   inline auto Mapping_Point(int, int, int) const;

   // j
   inline auto Jacobian() const;
   inline auto JacobianAtRadialNode(int, int) const;
   inline auto Jacobian_Point(int, int, int) const;

   // F^{-1}
   inline auto InverseF() const;
   inline const std::vector<std::vector<std::vector<Eigen::Matrix3cd>>> &
   InverseFRef() const;
   inline auto InverseFAtRadialNode(int, int) const;
   inline auto InverseF_Point(int, int, int) const;

   // a
   inline auto LaplaceTensor() const;
   inline const std::vector<std::vector<std::vector<Eigen::Matrix3cd>>> &
   LaplaceTensorRef() const;
   inline auto LaplaceTensorAtRadialNode(int, int) const;
   inline auto LaplaceTensor_Point(int, int, int) const;

   // physical radius
   inline auto PhysicalRadius_Line(int) const;
   inline auto PhysicalRadius_Point(int, int, int) const;

   // volume and mass etc
   inline auto Volume() const;
   inline auto Mass() const;
   // converting between form of solution as single Eigen::Vector to the
   // "standard" vector of vectors
   inline auto SingleEigenVectorToGeneralFormat(const Eigen::VectorXcd &) const;
   inline auto
   PowerSTD(const std::vector<std::vector<std::vector<std::complex<double>>>>
                &vec_pot) const;
   // lower triangle
   inline void SetLowerTriangle();
   //    void SetLowerTriangle() { _spectral_info.set_left_lower(); };

   // rotation to output a specific slice
   inline auto RotateSliceToEquator(
       std::vector<double> &, std::vector<double> &,
       const std::vector<std::vector<std::vector<std::complex<double>>>> &);
   inline void PhysicalOutputRotated(
       const std::string, std::vector<double> &, std::vector<double> &,
       const std::vector<std::vector<std::vector<std::complex<double>>>> &)
       const;
   inline void ReferentialOutputRotated(
       const std::string, std::vector<double> &, std::vector<double> &,
       const std::vector<std::vector<std::vector<std::complex<double>>>> &)
       const;
   // output
   // void OutputFullInformation(
   //     const std::string,
   //     const std::vector<std::vector<std::vector<std::complex<double>>>>
   //     &) const;
   inline void ReferentialOutputAtElement(
       const std::string,
       const std::vector<std::vector<std::vector<std::complex<double>>>> &,
       bool) const;
   inline void PhysicalOutputAtElement(
       const std::string,
       const std::vector<std::vector<std::vector<std::complex<double>>>> &,
       bool) const;
   inline void PhysicalOutputSlice(
       const std::string,
       const std::vector<std::vector<std::vector<std::complex<double>>>> &,
       bool) const;
   //    void SphericalHarmonicOutputAtPhysicalRadius(
   //        const std::string,
   //        const std::vector<std::vector<std::vector<std::complex<double>>>>
   //        &) const;
   inline void ReferentialOutputSlice(
       const std::string,
       const std::vector<std::vector<std::vector<std::complex<double>>>> &,
       bool) const;
   inline void OutputAtElement(
       const std::string,
       const std::vector<std::vector<std::complex<double>>> &) const;
   inline void PhysicalOutputAtElement(
       const std::string,
       const std::vector<std::vector<std::complex<double>>> &) const;
   inline void OutputAtElement(const std::string, const std::vector<double> &) const;
   inline void OutputAtElement(const std::string, const std::vector<double> &,
                        const std::vector<double> &) const;
   inline void OutputAtElement(const std::string,
                        const std::vector<std::vector<double>> &,
                        const std::vector<std::vector<double>> &) const;
   inline void CartesianOutputAtElement(
       const std::string,
       const std::vector<std::vector<std::complex<double>>> &) const;
   inline void CartesianOutputAtElement(
       const std::string,
       const std::vector<std::vector<std::vector<std::complex<double>>>> &)
       const;
   inline void ModelDensityOutputRotated(const std::string, std::vector<double> &,
                                  std::vector<double> &, bool) const;

   // norms
   inline double const LengthNorm() const;
   inline double const MassNorm() const;
   inline double const TimeNorm() const;
   inline double const DensityNorm() const;
   inline double const InertiaNorm() const;
   inline double const VelocityNorm() const;
   inline double const AccelerationNorm() const;
   inline double const ForceNorm() const;
   inline double const StressNorm() const;
   inline double const GravitationalConstantNorm() const;
   inline double const GravitationalConstant() const;
   inline double const PotentialNorm() const;

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
       gravitational_constant, gravitational_constant_norm, potential_norm;
   bool _isisotropic;
   std::vector<std::vector<bool>> _vec_issolid;
   std::vector<int> _vec_discontinuity_indices, _vec_laynum;
   //    vecdb _vec_nodes, _vec_allnodes;

   // model parameters
   vvvecdb _vec_h, _vec_density;
   std::vector<std::vector<std::vector<double>>> _vec_j;
   vvveceig _vec_invF, _vec_a;

   GaussQuad::Quadrature1D<double> _q;
   Radial_Tools::RadialMesh node_data;
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _mat_gaussderiv;
   SpectralElementTools::MatrixWeakForm _spectral_info;

   // grid
   Grid _grid;

   // full constructor
   template <class model, class tomomodel>
      requires PlanetaryModel::BasicSphericalDensityModel<model, int, double> &&
               PlanetaryModel::TomographyModel<tomomodel>
   Density3D(const model &, const tomomodel &, const int, const int, double,
             double);
};

};   // namespace GeneralEarthModels
#endif