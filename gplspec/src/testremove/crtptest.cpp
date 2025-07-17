#include <cmath>
#include <functional>
#include <iostream>

#include "Interpolation/All"

// The Curiously Recurring Template Pattern (CRTP)

// Declare a base_traits traits class template:
// Declare a base_traits traits class template:
template <typename Derived> struct base_traits;

// Define the base class that uses the traits:
template <typename Derived> struct NormBase {
   typedef typename base_traits<Derived>::value_type value_type;

   Derived &derived() { return *static_cast<Derived *>(this); };
   const Derived &derived() const {
      return *static_cast<const Derived *>(this);
   };

   value_type LengthNorm() { return derived().LengthNorm(); };
   value_type MassNorm() const { return derived().MassNorm(); };
   value_type TimeNorm() const { return derived().TimeNorm(); }

   value_type DensityNorm() const { return derived().DensityNorm(); };
   value_type InertiaNorm() const { return derived().InertiaNorm(); };
   value_type VelocityNorm() const { return derived().VelocityNorm(); };
   value_type AccelerationNorm() const { return derived().AccelerationNorm(); };
   value_type ForceNorm() const { return derived().ForceNorm(); };
   value_type StressNorm() const { return derived().StressNorm(); };
   value_type GravitationalConstant() const {
      return derived().GravitationalConstant();
   };
};

template <typename Derived> class SphericalGeometry : public NormBase<Derived> {
 protected:
   typedef typename base_traits<Derived>::value_type value_type;
   typedef NormBase<Derived> Base;

 public:
   using Base::derived;
   using retfunc = std::function<value_type(int)>;

   // model layers
   int NumberOfLayers() const { return derived().NumberOfLayers(); };
   value_type LowerRadius(int i) const { return derived().LowerRadius(i); }
   value_type UpperRadius(int i) const { return derived().UpperRadius(i); }
   value_type OuterRadius() const { return derived().OuterRadius(); }
};

template <typename Derived>
class SphericalDensityModel : public SphericalGeometry<Derived> {
 protected:
   typedef typename base_traits<Derived>::value_type value_type;
   typedef SphericalGeometry<Derived> Base;

 public:
   using Base::derived;
   using retfunc = std::function<value_type(value_type)>;

   // function returns
   retfunc Density(int i) const { return derived().Density(i); }
};

template <typename Derived>
class SphericalSeismicModel : public SphericalDensityModel<Derived> {
 protected:
   typedef typename base_traits<Derived>::value_type value_type;
   typedef SphericalDensityModel<Derived> Base;

 public:
   using Base::derived;
   using retfunc = std::function<value_type(value_type)>;

   // boolean return types
   bool IsIsotropic() const { return derived().IsIsotropic(); };
   bool IsSolid(int i) const { return derived().IsSolid(i); };
   bool IsFluid(int i) const { return derived().IsFluid(i); };

   // function returns
   retfunc VP(int i) const { return derived().VP(i); }
   retfunc VPV(int i) const { return derived().VPV(i); }
   retfunc VPH(int i) const { return derived().VPH(i); }
   retfunc VS(int i) const { return derived().VS(i); }
   retfunc VSV(int i) const { return derived().VSV(i); }
   retfunc VSH(int i) const { return derived().VSH(i); }

   // anisotropic functions
   retfunc Eta(int i) const { return derived().Eta(i); }
   retfunc Kappa(int i) const { return derived().Kappa(i); }
   retfunc Mu(int i) const { return derived().Mu(i); }
   retfunc A(int i) const { return derived().A(i); }
   retfunc C(int i) const { return derived().C(i); }
   retfunc N(int i) const { return derived().N(i); }
   retfunc L(int i) const { return derived().L(i); }
   retfunc F(int i) const { return derived().F(i); }
};

template <typename Derived>
class AsphericalSeismicModel : public SphericalSeismicModel<Derived> {
 protected:
   typedef typename base_traits<Derived>::value_type value_type;
   typedef SphericalSeismicModel<Derived> Base;

 public:
   using Base::derived;

   std::function<value_type(value_type, value_type)>
   RadialMap(value_type r) const {
      return derived().RadialMap(r);
   };
};

template <typename FLOAT>
class EarthNorms : public NormBase<EarthNorms<FLOAT>> {
 public:
   typedef typename base_traits<EarthNorms>::value_type value_type;
   value_type LengthNorm() const { return length_norm; };
   value_type MassNorm() const { return mass_norm; };
   value_type TimeNorm() const { return time_norm; }

   value_type DensityNorm() const { return density_norm; };
   value_type InertiaNorm() const { return inertia_norm; };
   value_type VelocityNorm() const { return velocity_norm; };
   value_type AccelerationNorm() const { return acceleration_norm; };
   value_type ForceNorm() const { return force_norm; };
   value_type StressNorm() const { return stress_norm; };
   value_type GravitationalConstant() const { return gravitational_constant; };

 private:
   const value_type length_norm = 6.371 * std::pow(10.0, 6.0);
   const value_type mass_norm = 5.972 * std::pow(10.0, 24.0);
   const value_type time_norm = 3600.0;

   const value_type density_norm = mass_norm / std::pow(length_norm, 3.0);
   const value_type inertia_norm = mass_norm * std::pow(length_norm, 2.0);
   const value_type velocity_norm = length_norm / time_norm;
   const value_type acceleration_norm = length_norm / std::pow(time_norm, 2.0);
   const value_type force_norm =
       mass_norm * length_norm / std::pow(time_norm, 2.0);
   const value_type stress_norm =
       mass_norm / (std::pow(time_norm, 2.0) * length_norm);
   const value_type gravitational_constant =
       std::pow(length_norm, 3.0) / (mass_norm * std::pow(time_norm, 2.0));
};

template <typename FLOAT> class PREM : public SphericalGeometry<PREM<FLOAT>> {
 public:
   typedef typename base_traits<PREM>::value_type value_type;
   using retfunc = std::function<value_type(int)>;

   // norms
   // using _normvalues::LengthNorm;
   value_type LengthNorm() const { return _normvalues.LengthNorm(); };
   value_type MassNorm() const { return _normvalues.MassNorm(); };
   value_type TimeNorm() const { return _normvalues.TimeNorm(); }
   value_type DensityNorm() const { return _normvalues.TimeNorm(); };
   value_type InertiaNorm() const { return _normvalues.InertiaNorm(); };
   value_type VelocityNorm() const { return _normvalues.VelocityNorm(); };
   value_type AccelerationNorm() const {
      return _normvalues.AccelerationNorm();
   };
   value_type ForceNorm() const { return _normvalues.ForceNorm(); };
   value_type StressNorm() const { return _normvalues.StressNorm(); };
   value_type GravitationalConstant() const {
      return _normvalues.GravitationalConstant();
   };
   Interpolation::Polynomial1D<FLOAT> Density(int i) const {
      return 1000.0 * vec_density[i];
   };
   // Geometry of PREM
   int NumberOfLayers() const { return 13; };
   value_type LowerRadius(int i) const { return vec_radii[i]; }
   value_type UpperRadius(int i) const { return vec_radii[i + 1]; }
   value_type OuterRadius() const { return vec_radii[13]; }

 private:
   std::vector<FLOAT> vec_radii{0.0,       1221500.0, 3480000.0, 3630000.0,
                                5600000.0, 5701000.0, 5771000.0, 5971000.0,
                                6151000.0, 6291000.0, 6346600.0, 6356000.0,
                                6368000.0, 6371000.0};
   EarthNorms<FLOAT> _normvalues;
   std::vector<Interpolation::Polynomial1D<FLOAT>> vec_density{
       {13.0855, 0, -8.8381},
       {12.5815, -1.2638, -3.6426, -5.5281},
       {7.9565, -6.4761, 5.5283, -3.0807},
       {7.9565, -6.4761, 5.5283, -3.0807},
       {7.9565, -6.4761, 5.5283, -3.0807},
       {5.3197, -1.4836},
       {11.2494, -8.0298},
       {7.1089, -3.8045},
       {2.6910, 0.6924},
       {2.6910, 0.6924},
       {2.900},
       {2.600},
       {1.020}};

   std::vector<Interpolation::Polynomial1D<FLOAT>> vec_p_velocity{
       {11.2622, 0, -6.3640},
       {11.0487, -4.0362, 4.8023, -13.5732},
       {15.3891, -5.3181, 5.5242, -2.5514},
       {24.9520, -40.4673, 51.4832, -26.6419},
       {29.2766, -23.6027, 5.5242, -2.5514},
       {19.0957, -9.8672},
       {39.7027, -32.6166},
       {20.3926, -12.2569},
       {4.1875, 3.9382},
       {4.1875, 3.9382},
       {6.800},
       {5.800},
       {1.450}};
   std::vector<Interpolation::Polynomial1D<FLOAT>> vec_pv_velocity{
       {11.2622, 0, -6.3640},
       {11.0487, -4.0362, 4.8023, -13.5732},
       {15.3891, -5.3181, 5.5242, -2.5514},
       {24.9520, -40.4673, 51.4832, -26.6419},
       {29.2766, -23.6027, 5.5242, -2.5514},
       {19.0957, -9.8672},
       {39.7027, -32.6166},
       {20.3926, -12.2569},
       {0.8317, 7.2180},
       {0.8317, 7.2180},
       {6.800},
       {5.800},
       {1.450}};
   std::vector<Interpolation::Polynomial1D<FLOAT>> vec_ph_velocity{
       {11.2622, 0, -6.3640},
       {11.0487, -4.0362, 4.8023, -13.5732},
       {15.3891, -5.3181, 5.5242, -2.5514},
       {24.9520, -40.4673, 51.4832, -26.6419},
       {29.2766, -23.6027, 5.5242, -2.5514},
       {19.0957, -9.8672},
       {39.7027, -32.6166},
       {20.3926, -12.2569},
       {3.5908, 4.6172},
       {3.5908, 4.6172},
       {6.800},
       {5.800},
       {1.450}};

   std::vector<Interpolation::Polynomial1D<FLOAT>> vec_s_velocity{
       {3.6678, 0, -4.4475},
       {0.0},
       {6.9254, 1.4672, -2.0834, 0.9783},
       {11.1671, -13.7818, 17.4575, -9.2777},
       {22.3459, -17.2473, -2.0834, 0.9783},
       {9.9839, -4.9324},
       {22.3512, -18.5856},
       {8.9496, -4.4597},
       {2.1519, 2.3481},
       {2.1519, 2.3481},
       {3.900},
       {3.200},
       {0}};
   std::vector<Interpolation::Polynomial1D<FLOAT>> vec_sv_velocity{
       {3.6678, 0, -4.4475},
       {0.0},
       {6.9254, 1.4672, -2.0834, 0.9783},
       {11.1671, -13.7818, 17.4575, -9.2777},
       {22.3459, -17.2473, -2.0834, 0.9783},
       {9.9839, -4.9324},
       {22.3512, -18.5856},
       {8.9496, -4.4597},
       {5.8582, -1.4678},
       {5.8582, -1.4678},
       {3.900},
       {3.200},
       {0}};
   std::vector<Interpolation::Polynomial1D<FLOAT>> vec_sh_velocity{
       {3.6678, 0, -4.4475},
       {0.0},
       {6.9254, 1.4672, -2.0834, 0.9783},
       {11.1671, -13.7818, 17.4575, -9.2777},
       {22.3459, -17.2473, -2.0834, 0.9783},
       {9.9839, -4.9324},
       {22.3512, -18.5856},
       {8.9496, -4.4597},
       {-1.0839, 5.7176},
       {-1.0839, 5.7176},
       {3.900},
       {3.200},
       {0}};
   std::vector<Interpolation::Polynomial1D<FLOAT>> vec_Qmu{
       {84.6}, {std::pow(10.0, 10.0)},
       {312},  {312},
       {312},  {143},
       {143},  {143},
       {80},   {600},
       {600},  {std::pow(10.0, 10.0)}};

   std::vector<Interpolation::Polynomial1D<FLOAT>> vec_QKappa{
       {1327.7}, {57823}, {57823}, {57823}, {57823}, {57823},
       {57823},  {57823}, {57823}, {57823}, {57823}, {57823}};

   std::vector<Interpolation::Polynomial1D<FLOAT>> vec_eta{{1},
                                                           {1},
                                                           {1},
                                                           {1},
                                                           {1},
                                                           {1},
                                                           {1},
                                                           {1},
                                                           {1},
                                                           {1},
                                                           {3.3687, -2.4778},
                                                           {3.3687, -2.4778},
                                                           {1},
                                                           {1},
                                                           {1}};

   std::vector<Interpolation::Polynomial1D<FLOAT>> vec_A;
};

// Declare and define a base_traits specialization for derived:
template <typename FLOAT> struct base_traits<EarthNorms<FLOAT>> {
   typedef FLOAT value_type;
};

// Declare and define a base_traits specialization for derived:
template <typename FLOAT> struct base_traits<PREM<FLOAT>> {
   typedef FLOAT value_type;
};

// template <typename FLOAT> struct base_traits<SphericalGeometry<FLOAT>> {
//    typedef FLOAT value_type;
// };

int
main() {

   EarthNorms<double> myEarth;
   PREM<double> SphEarth;

   std::cout << " " << myEarth.LengthNorm() << std::endl;
   std::cout << " " << SphEarth.LengthNorm() << " "
             << SphEarth.GravitationalConstant() << std::endl;
   std::cout << " " << SphEarth.NumberOfLayers() << " "
             << SphEarth.LowerRadius(2) << " " << SphEarth.UpperRadius(2) << " "
             << SphEarth.OuterRadius() << std::endl;

   return 0;
}
