#ifndef GPLSPEC_FIELD_BUNDLE_H
#define GPLSPEC_FIELD_BUNDLE_H

#include "Field_Evaluation.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <utility>

namespace GPLSpec::FieldBundle {

struct InputIdentity {
   std::string path;
   std::string sha256;
};

inline InputIdentity InputPath(const std::string &path) {
   if (path.empty()) throw std::invalid_argument("input path is required");
   return {path, {}};
}

struct RunMetadata {
   std::string workflow;
   std::string dataset;
   std::string solver_method;
   double solver_tolerance = 0.0;
   std::vector<InputIdentity> inputs;
   std::vector<std::pair<std::string, std::string>> parameters;
   std::vector<int> selector_l;
   std::vector<int> selector_m;
   std::vector<double> selector_weight;
   std::optional<FieldEvaluation::Rotation> rotation;
   struct Solve {
      int status = -1;
      long iterations = 0;
      double tolerance = 0.0;
      double recursive_error = 0.0;
      double residual_norm = 0.0;
      double relative_residual = 0.0;
   };
   std::vector<Solve> solves;
};

inline void CopyInputProvenance(const RunMetadata &source,
                                RunMetadata &derived) {
   derived.inputs = source.inputs;
   derived.parameters = source.parameters;
}

inline std::string JsonString(const std::string &value) {
   std::ostringstream out;
   out << '"';
   for (const unsigned char c : value) {
      switch (c) {
      case '"': out << "\\\""; break;
      case '\\': out << "\\\\"; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
         if (c < 0x20) throw std::invalid_argument("control byte in metadata string");
         out << static_cast<char>(c);
      }
   }
   out << '"';
   return out.str();
}

inline std::string ParameterNumber(double value) {
   if (!std::isfinite(value))
      throw std::invalid_argument("non-finite run parameter");
   std::ostringstream out;
   out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
   return out.str();
}

inline void Validate(const RunMetadata &metadata) {
   if (metadata.workflow.empty() || metadata.dataset.empty())
      throw std::invalid_argument("bundle workflow and dataset are required");
   if (!metadata.solver_method.empty() &&
       (!std::isfinite(metadata.solver_tolerance) || metadata.solver_tolerance <= 0.0))
      throw std::invalid_argument("solver tolerance must be finite and positive");
   if (metadata.selector_l.size() != metadata.selector_m.size() ||
       metadata.selector_l.size() != metadata.selector_weight.size())
      throw std::invalid_argument("sensitivity selector columns must have equal length");
   if (metadata.rotation)
      FieldEvaluation::detail::ValidateRotation(*metadata.rotation);
   for (const auto &input : metadata.inputs) {
      if (input.path.empty() || (!input.sha256.empty() &&
          (input.sha256.size() != 64 ||
          !std::all_of(input.sha256.begin(), input.sha256.end(), [](unsigned char c) {
             return std::isxdigit(c) != 0;
          }))))
         throw std::invalid_argument("input identity requires a path and optional SHA-256 digest");
   }
   for (double weight : metadata.selector_weight)
      if (!std::isfinite(weight)) throw std::invalid_argument("non-finite selector weight");
   for (const auto &solve : metadata.solves)
      if (solve.status < 0 || solve.iterations < 0 ||
          !std::isfinite(solve.tolerance) || solve.tolerance <= 0.0 ||
          !std::isfinite(solve.recursive_error) ||
          !std::isfinite(solve.residual_norm) ||
          !std::isfinite(solve.relative_residual))
         throw std::invalid_argument("invalid solve diagnostics in bundle metadata");
}

inline const char *QuantityName(FieldEvaluation::Quantity q) {
   using FieldEvaluation::Quantity;
   switch (q) {
   case Quantity::Potential: return "potential";
   case Quantity::SensitivityKernel: return "sensitivity_kernel";
   case Quantity::ReferentialDensity: return "density_referential";
   case Quantity::PhysicalDensity: return "density_physical";
   case Quantity::DegreePower: return "power_std";
   case Quantity::SurfaceRadius: return "surface_radius";
   }
   throw std::invalid_argument("unknown quantity");
}

inline const char *RepresentationName(FieldEvaluation::Representation r) {
   using FieldEvaluation::Representation;
   switch (r) {
   case Representation::AngularGrid: return "angular_grid";
   case Representation::HarmonicCoefficients: return "harmonic_coefficients";
   case Representation::RotatedSlice: return "rotated_slice";
   case Representation::DegreePower: return "degree_power";
   case Representation::RadialProfile: return "radial_profile";
   }
   throw std::invalid_argument("unknown representation");
}

inline const char *BoundaryName(FieldEvaluation::BoundaryKind b) {
   using FieldEvaluation::BoundaryKind;
   switch (b) {
   case BoundaryKind::None: return "none";
   case BoundaryKind::Numerical: return "numerical";
   case BoundaryKind::Material: return "material";
   case BoundaryKind::BodyVacuum: return "body_vacuum";
   case BoundaryKind::SurfaceWithoutExterior: return "surface_without_exterior";
   }
   throw std::invalid_argument("unknown boundary kind");
}

inline const char *SideName(FieldEvaluation::TraceSide side) {
   using FieldEvaluation::TraceSide;
   switch (side) {
   case TraceSide::None: return "none";
   case TraceSide::Inner: return "inner";
   case TraceSide::Outer: return "outer";
   }
   throw std::invalid_argument("unknown trace side");
}

inline std::string OptionalIndex(std::size_t value) {
   return value == FieldEvaluation::NoIndex ? "" : std::to_string(value);
}

inline void WriteSample(std::ostream &out, const FieldEvaluation::Sample &s) {
   const double values[] = {s.reference_radius_m, s.physical_radius_m,
                            s.theta_rad, s.phi_rad,
                            s.value.real(), s.value.imag()};
   for (double value : values)
      if (!std::isfinite(value) && !std::isnan(value))
         throw std::invalid_argument("infinite sample value cannot be serialized");
   if (!std::isfinite(s.value.real()) || !std::isfinite(s.value.imag()))
      throw std::invalid_argument("sample value must be finite");
   if (s.radial_node < -1)
      throw std::invalid_argument("radial node index must be nonnegative or absent");
   out << QuantityName(s.quantity) << ',' << RepresentationName(s.representation)
       << ',' << BoundaryName(s.boundary) << ',' << SideName(s.side) << ','
       << OptionalIndex(s.element) << ',' << OptionalIndex(s.numerical_interface)
       << ',' << OptionalIndex(s.material_interface) << ','
       << OptionalIndex(s.inner_element) << ',' << OptionalIndex(s.outer_element)
       << ',' << s.material << ',' << s.inner_material << ',' << s.outer_material
       << ',' << s.degree << ',' << s.order << ',';
   if (s.radial_node >= 0) out << s.radial_node;
   out << ','
       << std::setprecision(std::numeric_limits<double>::max_digits10);
   auto number = [&](double value) { if (std::isfinite(value)) out << value; };
   number(s.reference_radius_m); out << ','; number(s.physical_radius_m); out << ',';
   number(s.theta_rad); out << ','; number(s.phi_rad); out << ','
                            << s.value.real() << ',' << s.value.imag() << '\n';
}

class Writer {
public:
   Writer(const std::filesystem::path &directory,
          const GeneralEarthModels::Density3D &model, RunMetadata metadata)
       : directory_(directory), model_(model), metadata_(std::move(metadata)) {
      Validate(metadata_);
      std::filesystem::create_directories(directory_);
   }

   template<class Evaluator>
   void Write(Evaluator &&evaluate) const {
      const auto manifest_tmp = directory_ / "manifest.json.tmp";
      const auto samples_tmp = directory_ / "samples.csv.tmp";
      std::ofstream manifest(manifest_tmp, std::ios::binary | std::ios::trunc);
      std::ofstream samples(samples_tmp, std::ios::binary | std::ios::trunc);
      if (!manifest || !samples) throw std::runtime_error("cannot create field bundle");
      WriteManifest(manifest);
      samples << "quantity,representation,boundary,side,element,numerical_interface,"
                 "material_interface,inner_element,outer_element,material,inner_material,"
                 "outer_material,degree,order,radial_node,reference_radius_m,"
                 "physical_radius_m,theta_rad,phi_rad,value_real,value_imag\n";
      evaluate([&](const FieldEvaluation::Sample &sample) { WriteSample(samples, sample); });
      manifest.close(); samples.close();
      if (!manifest || !samples) throw std::runtime_error("failed writing field bundle");
      std::filesystem::rename(samples_tmp, directory_ / "samples.csv");
      std::filesystem::rename(manifest_tmp, directory_ / "manifest.json");
   }

private:
   void WriteManifest(std::ostream &out) const {
      out << std::setprecision(std::numeric_limits<double>::max_digits10)
          << "{\n  \"schema_version\": 1,\n  \"workflow\": "
          << JsonString(metadata_.workflow) << ",\n  \"dataset\": "
          << JsonString(metadata_.dataset)
          << ",\n  \"quantity_contracts\": {\n"
             "    \"potential\": {\"unit\": \"m^2 s^-2\", \"meaning\": \"Gravitational potential; harmonic coefficients are multiplied by PotentialNorm.\", \"representations\": {\"angular_grid\": \"Mapped physical spherical samples; reference and physical radii are both recorded.\", \"harmonic_coefficients\": \"GSHTrans orthonormal complex coefficients at the named radial boundary.\", \"rotated_slice\": \"Physical potential samples in the rotated equatorial plane.\", \"radial_profile\": \"Physical scalar potential on the named radial ray; both radii and angles are recorded.\"}},\n"
             "    \"sensitivity_kernel\": {\"unit\": \"m^2 s^-2 kg^-1\", \"meaning\": \"Derivative of the selected potential harmonic with respect to referential source density, integrated with dV_ref; this is a kernel, not a potential field.\", \"representations\": {\"angular_grid\": \"Mapped physical spherical locations with the referential integration frame retained.\", \"harmonic_coefficients\": \"GSHTrans orthonormal complex kernel coefficients at the named radial boundary.\", \"rotated_slice\": \"Kernel samples in the rotated equatorial plane; referential integration meaning is unchanged.\"}},\n"
             "    \"density_referential\": {\"unit\": \"kg m^-3\", \"meaning\": \"Mass density per reference volume.\", \"representations\": {\"angular_grid\": \"Element-local density at mapped physical spherical sample locations.\", \"harmonic_coefficients\": \"GSHTrans coefficients of referential density at the named radial boundary.\", \"rotated_slice\": \"Referential density evaluated at rotated physical sample locations.\"}},\n"
             "    \"density_physical\": {\"unit\": \"kg m^-3\", \"meaning\": \"Mass density per mapped physical volume; rho_physical=rho_referential/J.\", \"representations\": {\"angular_grid\": \"Element-local physical density at mapped physical spherical sample locations.\", \"harmonic_coefficients\": \"GSHTrans coefficients of physical density at the named radial boundary.\", \"rotated_slice\": \"Physical density evaluated at rotated physical sample locations.\"}},\n"
             "    \"power_std\": {\"unit\": \"1\", \"meaning\": \"PowerSTD_l=sum_m |Phi_lm|^2/|Phi_00|^2 at each radial sample.\", \"representations\": {\"degree_power\": \"Dimensionless normalized power by degree, not potential amplitude.\"}},\n"
             "    \"surface_radius\": {\"unit\": \"m\", \"meaning\": \"Physical surface radius at the requested angular sample.\", \"representations\": {\"angular_grid\": \"Explicit latitude/longitude or GSH-grid angles in the physical surface frame.\"}}\n"
             "  },\n"
             "  \"coordinates\": {\"reference\": \"spherical radius m\", \"physical\": \"spherical radius m\", \"angles\": \"theta colatitude and phi longitude in rad\", \"cartesian\": \"derive x=r sin(theta) cos(phi), y=r sin(theta) sin(phi), z=r cos(theta)\"},\n"
             "  \"harmonics\": {\"convention\": \"GSHTrans orthonormal complex Y_lm\", \"index\": \"l*l+l+m\", \"degree_max\": "
          << model_.GSH_Grid().MaxDegree() << ", \"length_norm_m\": "
          << model_.LengthNorm() << ", \"mass_norm_kg\": " << model_.MassNorm()
          << ", \"time_norm_s\": " << model_.TimeNorm()
          << ", \"potential_norm_m2_s-2\": " << model_.PotentialNorm()
          << ", \"density_norm_kg_m-3\": " << model_.DensityNorm()
          << "},\n  \"solver\": {\"method\": "
          << JsonString(metadata_.solver_method) << ", \"requested_tolerance\": ";
      if (metadata_.solver_method.empty()) out << "null";
      else out << metadata_.solver_tolerance;
      out << "},\n  \"inputs\": [";
      for (std::size_t i = 0; i < metadata_.inputs.size(); ++i) {
         if (i) out << ',';
         out << "{\"path\":" << JsonString(metadata_.inputs[i].path)
             << ",\"sha256\":";
         if (metadata_.inputs[i].sha256.empty()) out << "null";
         else out << JsonString(metadata_.inputs[i].sha256);
         out << '}';
      }
      out << "],\n  \"parameters\": [";
      for (std::size_t i = 0; i < metadata_.parameters.size(); ++i) {
         if (i) out << ',';
         out << "{\"name\":" << JsonString(metadata_.parameters[i].first)
             << ",\"value\":" << JsonString(metadata_.parameters[i].second) << '}';
      }
      out << "],\n  \"selector\": [";
      for (std::size_t i = 0; i < metadata_.selector_l.size(); ++i) {
         if (i) out << ',';
         out << "{\"l\":" << metadata_.selector_l[i] << ",\"m\":"
             << metadata_.selector_m[i] << ",\"weight\":"
             << metadata_.selector_weight[i] << '}';
      }
      out << "],\n  \"rotation\":{\"angles_rad\":";
      if (!metadata_.rotation) out << "null";
      else out << '[' << metadata_.rotation->alpha << ',' << metadata_.rotation->beta
               << ',' << metadata_.rotation->gamma << ']';
      out << ",\"convention\":\"F_rot(n)=F(Rz(alpha)Ry(beta)Rz(gamma)n)\","
             "\"frame\":\"right-handed body-fixed Cartesian frame\","
             "\"legacy_harmonic_coefficient_scale\":"
          << std::sqrt(3.1415926535 / std::numbers::pi_v<double>)
          << ",\"scale_note\":\"Preserved RotationOutputMatrices normalization; it differs from unity by the legacy pi approximation.\"},\n";
      out << "\n  \"solve_diagnostics\": [";
      for (std::size_t i = 0; i < metadata_.solves.size(); ++i) {
         if (i) out << ',';
         const auto &solve = metadata_.solves[i];
         out << "{\"status\":" << solve.status << ",\"iterations\":"
             << solve.iterations << ",\"tolerance\":" << solve.tolerance
             << ",\"recursive_error\":" << solve.recursive_error
             << ",\"residual_norm\":" << solve.residual_norm
             << ",\"relative_residual\":" << solve.relative_residual << '}';
      }
      out << "],\n  \"trace_contract\": \"side is by increasing reference radius; potential uses shared global DOFs; density is element-local\",\n"
             "  \"sample_columns\": \"See samples.csv header; empty optional geometry is inapplicable, never zero\"\n}\n";
   }

   std::filesystem::path directory_;
   const GeneralEarthModels::Density3D &model_;
   RunMetadata metadata_;
};

inline RunMetadata WithSolve(const Gravity_Tools::SolveResult &result,
                             std::string workflow, std::string dataset,
                             std::string method, double tolerance) {
   RunMetadata metadata;
   metadata.workflow = std::move(workflow);
   metadata.dataset = std::move(dataset);
   metadata.solver_method = std::move(method);
   metadata.solver_tolerance = tolerance;
   for (const auto &solve : result.solves) {
      metadata.solves.push_back({static_cast<int>(solve.status),
          static_cast<long>(solve.iterations), solve.requested_tolerance,
          solve.recursive_error, solve.residual_norm, solve.relative_residual});
   }
   return metadata;
}

inline void WriteAngular(const std::filesystem::path &directory,
                         const GeneralEarthModels::Density3D &model,
                         std::span<const std::complex<double>> solution,
                         FieldEvaluation::Quantity quantity,
                         RunMetadata metadata,
                         const FieldEvaluation::Rotation *rotation = nullptr,
                         FieldEvaluation::TraceSelection traces =
                             FieldEvaluation::TraceSelection::Both) {
   metadata.rotation = rotation ?
       std::optional<FieldEvaluation::Rotation>(*rotation) : std::nullopt;
   Writer(directory, model, std::move(metadata)).Write([&](auto &&emit) {
      FieldEvaluation::ForEachAngularGridSample(model, solution, quantity,
                                                std::forward<decltype(emit)>(emit),
                                                rotation, traces);
   });
}

inline void WriteAngular(const std::filesystem::path &directory,
                         const GeneralEarthModels::Density3D &model,
                         const Eigen::VectorXcd &solution,
                         FieldEvaluation::Quantity quantity,
                         RunMetadata metadata,
                         const FieldEvaluation::Rotation *rotation = nullptr,
                         FieldEvaluation::TraceSelection traces =
                             FieldEvaluation::TraceSelection::Both) {
   WriteAngular(directory, model,
       std::span<const std::complex<double>>(solution.data(),
           static_cast<std::size_t>(solution.size())), quantity,
       std::move(metadata), rotation, traces);
}

inline void WriteRadialHarmonics(const std::filesystem::path &directory,
                                const GeneralEarthModels::Density3D &model,
                                std::span<const std::complex<double>> solution,
                                FieldEvaluation::Quantity quantity,
                                RunMetadata metadata) {
   Writer(directory, model, std::move(metadata)).Write([&](auto &&emit) {
      FieldEvaluation::ForEachRadialHarmonicSample(
          model, solution, quantity, std::forward<decltype(emit)>(emit));
   });
}

inline void WriteRadialHarmonics(const std::filesystem::path &directory,
                                const GeneralEarthModels::Density3D &model,
                                const Eigen::VectorXcd &solution,
                                FieldEvaluation::Quantity quantity,
                                RunMetadata metadata) {
   WriteRadialHarmonics(directory, model,
       std::span<const std::complex<double>>(solution.data(),
           static_cast<std::size_t>(solution.size())), quantity,
       std::move(metadata));
}

inline void WriteSlice(const std::filesystem::path &directory,
                       const GeneralEarthModels::Density3D &model,
                       std::span<const std::complex<double>> solution,
                       FieldEvaluation::Quantity quantity,
                       const FieldEvaluation::Rotation &rotation,
                       const std::vector<double> &azimuths,
                       RunMetadata metadata) {
   metadata.rotation = rotation;
   Writer(directory, model, std::move(metadata)).Write([&](auto &&emit) {
      FieldEvaluation::ForEachRotatedSliceSample(
          model, solution, quantity, rotation, azimuths,
          std::forward<decltype(emit)>(emit));
   });
}

inline void WriteEquatorialSlice(const std::filesystem::path &directory,
                                 const GeneralEarthModels::Density3D &model,
                                 std::span<const std::complex<double>> solution,
                                 FieldEvaluation::Quantity quantity,
                                 const std::vector<double> &azimuths,
                                 RunMetadata metadata) {
   metadata.rotation.reset();
   Writer(directory, model, std::move(metadata)).Write([&](auto &&emit) {
      FieldEvaluation::ForEachEquatorialSliceSample(
          model, solution, quantity, azimuths,
          std::forward<decltype(emit)>(emit));
   });
}

inline void WriteEquatorialSlice(const std::filesystem::path &directory,
                                 const GeneralEarthModels::Density3D &model,
                                 const Eigen::VectorXcd &solution,
                                 FieldEvaluation::Quantity quantity,
                                 const std::vector<double> &azimuths,
                                 RunMetadata metadata) {
   WriteEquatorialSlice(directory, model,
       std::span<const std::complex<double>>(solution.data(),
           static_cast<std::size_t>(solution.size())), quantity, azimuths,
       std::move(metadata));
}

inline void WriteSlice(const std::filesystem::path &directory,
                       const GeneralEarthModels::Density3D &model,
                       const Eigen::VectorXcd &solution,
                       FieldEvaluation::Quantity quantity,
                       const FieldEvaluation::Rotation &rotation,
                       const std::vector<double> &azimuths,
                       RunMetadata metadata) {
   WriteSlice(directory, model,
       std::span<const std::complex<double>>(solution.data(),
           static_cast<std::size_t>(solution.size())), quantity, rotation,
       azimuths, std::move(metadata));
}

} // namespace GPLSpec::FieldBundle

#endif
