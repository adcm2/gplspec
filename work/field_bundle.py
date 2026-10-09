"""Strict streaming reader for GPLSpec field bundles (schema version 1)."""

import csv
import json
import math
from pathlib import Path

HEADER = [
    "quantity", "representation", "boundary", "side", "element",
    "numerical_interface", "material_interface", "inner_element",
    "outer_element", "material", "inner_material", "outer_material",
    "degree", "order", "radial_node", "reference_radius_m",
    "physical_radius_m", "theta_rad", "phi_rad", "value_real", "value_imag",
]
CONTRACTS = {
    "potential": ("m^2 s^-2", {"angular_grid", "harmonic_coefficients", "rotated_slice", "radial_profile"}),
    "sensitivity_kernel": ("m^2 s^-2 kg^-1", {"angular_grid", "harmonic_coefficients", "rotated_slice"}),
    "density_referential": ("kg m^-3", {"angular_grid", "harmonic_coefficients", "rotated_slice"}),
    "density_physical": ("kg m^-3", {"angular_grid", "harmonic_coefficients", "rotated_slice"}),
    "power_std": ("1", {"degree_power"}),
    "surface_radius": ("m", {"angular_grid"}),
}
QUANTITIES = set(CONTRACTS)
REPRESENTATIONS = {rep for _, reps in CONTRACTS.values() for rep in reps}
COMPATIBLE = {quantity: reps for quantity, (_, reps) in CONTRACTS.items()}
BOUNDARIES = {"none", "numerical", "material", "body_vacuum", "surface_without_exterior"}
SIDES = {"none", "inner", "outer"}
OPTIONAL_INDEX_COLUMNS = {
    "element", "numerical_interface", "material_interface", "inner_element",
    "outer_element", "radial_node",
}
INTEGER_COLUMNS = {"material", "inner_material", "outer_material", "degree", "order"}
FLOAT_COLUMNS = set(HEADER[15:])
EXPECTED_ROTATION = "F_rot(n)=F(Rz(alpha)Ry(beta)Rz(gamma)n)"
EXPECTED_ROTATION_FRAME = "right-handed body-fixed Cartesian frame"


def _finite_number(value, name, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"invalid numeric metadata field {name}")
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0.0):
        raise ValueError(f"invalid numeric metadata field {name}")
    return value


def _validate_manifest(manifest):
    if manifest.get("schema_version") != 1:
        raise ValueError("unsupported field-bundle schema version")
    for name in ("workflow", "dataset", "quantity_contracts", "coordinates", "harmonics",
                 "solver", "inputs", "parameters", "selector", "rotation",
                 "solve_diagnostics", "trace_contract", "sample_columns"):
        if name not in manifest:
            raise ValueError(f"manifest is missing {name}")
    if not isinstance(manifest["workflow"], str) or not manifest["workflow"]:
        raise ValueError("manifest workflow must be a non-empty string")
    if not isinstance(manifest["dataset"], str) or not manifest["dataset"]:
        raise ValueError("manifest dataset must be a non-empty string")

    contracts = manifest["quantity_contracts"]
    if not isinstance(contracts, dict) or set(contracts) != QUANTITIES:
        raise ValueError("manifest quantity contracts are incomplete")
    for quantity, (unit, representations) in CONTRACTS.items():
        contract = contracts[quantity]
        if not isinstance(contract, dict) or contract.get("unit") != unit:
            raise ValueError(f"invalid unit contract for {quantity}")
        if not isinstance(contract.get("meaning"), str) or not contract["meaning"].strip():
            raise ValueError(f"missing quantity meaning for {quantity}")
        declared = contract.get("representations")
        if (not isinstance(declared, dict) or set(declared) != representations or
                any(not isinstance(text, str) or not text.strip() for text in declared.values())):
            raise ValueError(f"invalid representation contract for {quantity}")
    if "rho_physical=rho_referential/J" not in contracts["density_physical"]["meaning"]:
        raise ValueError("physical-density contract must state rho_physical=rho_referential/J")
    if "PowerSTD_l=sum_m |Phi_lm|^2/|Phi_00|^2" not in contracts["power_std"]["meaning"]:
        raise ValueError("power contract must state the PowerSTD formula")
    if "dV_ref" not in contracts["sensitivity_kernel"]["meaning"]:
        raise ValueError("sensitivity contract must state its referential-volume measure")

    coordinates = manifest["coordinates"]
    required_coordinates = {"reference", "physical", "angles", "cartesian"}
    if not isinstance(coordinates, dict) or set(coordinates) != required_coordinates:
        raise ValueError("invalid coordinate contract")
    if any(not isinstance(value, str) or not value.strip() for value in coordinates.values()):
        raise ValueError("coordinate descriptions must be non-empty strings")

    harmonics = manifest["harmonics"]
    if not isinstance(harmonics, dict):
        raise ValueError("harmonics metadata must be an object")
    if harmonics.get("convention") != "GSHTrans orthonormal complex Y_lm" or harmonics.get("index") != "l*l+l+m":
        raise ValueError("unsupported harmonic convention or index")
    degree_max = harmonics.get("degree_max")
    if isinstance(degree_max, bool) or not isinstance(degree_max, int) or degree_max < 0:
        raise ValueError("invalid maximum harmonic degree")
    for name in ("length_norm_m", "mass_norm_kg", "time_norm_s",
                 "potential_norm_m2_s-2", "density_norm_kg_m-3"):
        _finite_number(harmonics.get(name), name, positive=True)

    rotation = manifest["rotation"]
    required_rotation = {"angles_rad", "convention", "frame", "legacy_harmonic_coefficient_scale", "scale_note"}
    if not isinstance(rotation, dict) or set(rotation) != required_rotation:
        raise ValueError("invalid rotation metadata")
    angles = rotation["angles_rad"]
    if angles is not None:
        if not isinstance(angles, list) or len(angles) != 3:
            raise ValueError("rotation angles must be null or three radians")
        for angle in angles:
            _finite_number(angle, "rotation angle")
    if rotation["convention"] != EXPECTED_ROTATION or rotation["frame"] != EXPECTED_ROTATION_FRAME:
        raise ValueError("unsupported rotation convention or frame")
    expected_scale = math.sqrt(3.1415926535 / math.pi)
    if not math.isclose(_finite_number(rotation["legacy_harmonic_coefficient_scale"], "rotation scale"),
                        expected_scale, rel_tol=0.0, abs_tol=2e-16):
        raise ValueError("unexpected legacy rotation coefficient scale")
    if not isinstance(rotation["scale_note"], str) or not rotation["scale_note"].strip():
        raise ValueError("missing rotation scale note")

    solver = manifest["solver"]
    if not isinstance(solver, dict) or set(solver) != {"method", "requested_tolerance"}:
        raise ValueError("invalid solver metadata")
    if not isinstance(solver["method"], str):
        raise ValueError("solver method must be a string")
    if solver["method"]:
        _finite_number(solver["requested_tolerance"], "requested solver tolerance", positive=True)
    elif solver["requested_tolerance"] is not None:
        raise ValueError("solver tolerance must be null when no solver method is declared")

    if not isinstance(manifest["inputs"], list) or not isinstance(manifest["selector"], list):
        raise ValueError("manifest inputs and selector must be lists")
    for source in manifest["inputs"]:
        if not isinstance(source, dict) or set(source) != {"path", "sha256"} or not isinstance(source["path"], str) or not source["path"]:
            raise ValueError("invalid input identity")
        digest = source["sha256"]
        if digest is not None and (not isinstance(digest, str) or len(digest) != 64 or
                                   any(ch not in "0123456789abcdefABCDEF" for ch in digest)):
            raise ValueError("invalid input SHA-256")
    if not isinstance(manifest["parameters"], list):
        raise ValueError("manifest parameters must be a list")
    for parameter in manifest["parameters"]:
        if (not isinstance(parameter, dict) or set(parameter) != {"name", "value"} or
                not isinstance(parameter["name"], str) or not parameter["name"] or
                not isinstance(parameter["value"], str)):
            raise ValueError("invalid run parameter")
    for selector in manifest["selector"]:
        if not isinstance(selector, dict) or set(selector) != {"l", "m", "weight"}:
            raise ValueError("invalid sensitivity selector")
        degree, order = selector["l"], selector["m"]
        if (isinstance(degree, bool) or not isinstance(degree, int) or degree < 0 or degree > degree_max or
                isinstance(order, bool) or not isinstance(order, int) or abs(order) > degree):
            raise ValueError("sensitivity selector is outside harmonic range")
        _finite_number(selector["weight"], "selector weight")
    if not isinstance(manifest["solve_diagnostics"], list):
        raise ValueError("solve diagnostics must be a list")
    for diagnostic in manifest["solve_diagnostics"]:
        fields = {"status", "iterations", "tolerance", "recursive_error", "residual_norm", "relative_residual"}
        if not isinstance(diagnostic, dict) or set(diagnostic) != fields:
            raise ValueError("invalid solve diagnostics")
        if (isinstance(diagnostic["status"], bool) or not isinstance(diagnostic["status"], int) or diagnostic["status"] < 0 or
                isinstance(diagnostic["iterations"], bool) or not isinstance(diagnostic["iterations"], int) or diagnostic["iterations"] < 0):
            raise ValueError("invalid solve status or iteration count")
        _finite_number(diagnostic["tolerance"], "solve tolerance", positive=True)
        for name in ("recursive_error", "residual_norm", "relative_residual"):
            _finite_number(diagnostic[name], name)
    if manifest["trace_contract"] != "side is by increasing reference radius; potential uses shared global DOFs; density is element-local":
        raise ValueError("unsupported trace contract")


def _row_key(row):
    def optional_integer(value):
        return None if value == "" else int(value)

    def number(value):
        if value == "":
            return None
        result = float(value)
        return 0.0 if result == 0.0 else result

    return (
        row["quantity"], row["representation"], row["boundary"], row["side"],
        optional_integer(row["element"]), optional_integer(row["numerical_interface"]),
        optional_integer(row["material_interface"]), optional_integer(row["inner_element"]),
        optional_integer(row["outer_element"]), int(row["degree"]), int(row["order"]),
        optional_integer(row["radial_node"]), number(row["reference_radius_m"]),
        number(row["physical_radius_m"]), number(row["theta_rad"]), number(row["phi_rad"]),
    )


def open_bundle(directory):
    """Validate the manifest and return (manifest, streaming row iterator)."""
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    _validate_manifest(manifest)
    degree_max = manifest["harmonics"]["degree_max"]
    stream = (directory / "samples.csv").open(newline="", encoding="utf-8")
    reader = csv.DictReader(stream)
    if reader.fieldnames != HEADER:
        stream.close()
        raise ValueError("sample CSV header does not match schema version 1")

    def rows():
        seen = set()
        try:
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError("malformed sample row")
                quantity, representation = row["quantity"], row["representation"]
                if quantity not in QUANTITIES or representation not in REPRESENTATIONS:
                    raise ValueError("unknown sample quantity or representation")
                if representation not in COMPATIBLE[quantity]:
                    raise ValueError("quantity and representation are incompatible")
                if row["boundary"] not in BOUNDARIES or row["side"] not in SIDES:
                    raise ValueError("unknown boundary or trace side")
                for column in OPTIONAL_INDEX_COLUMNS:
                    value = row[column]
                    if value and (not value.isdecimal()):
                        raise ValueError(f"invalid optional index {column}")
                for column in INTEGER_COLUMNS:
                    try:
                        numeric = int(row[column])
                    except ValueError as exc:
                        raise ValueError(f"invalid integer {column}") from exc
                    if str(numeric) != row[column]:
                        raise ValueError(f"non-canonical integer {column}")
                degree, order = int(row["degree"]), int(row["order"])
                if representation == "harmonic_coefficients":
                    if degree < 0 or degree > degree_max or abs(order) > degree:
                        raise ValueError("harmonic coefficient index is outside declared range")
                elif representation == "degree_power":
                    if degree < 0 or degree > degree_max or order != 0:
                        raise ValueError("degree-power index is outside declared range")
                elif degree != -1 or order != 0:
                    raise ValueError("spatial samples must use degree=-1 and order=0")
                numeric_values = {}
                for column in FLOAT_COLUMNS:
                    if row[column] == "":
                        numeric_values[column] = None
                        continue
                    try:
                        number = float(row[column])
                    except ValueError as exc:
                        raise ValueError(f"invalid number {column}") from exc
                    if not math.isfinite(number):
                        raise ValueError(f"non-finite number {column}")
                    numeric_values[column] = number
                if numeric_values["value_real"] is None or numeric_values["value_imag"] is None:
                    raise ValueError("sample value is required")
                if representation in {"angular_grid", "rotated_slice"}:
                    for coordinate in ("reference_radius_m", "physical_radius_m", "theta_rad", "phi_rad"):
                        if numeric_values[coordinate] is None:
                            raise ValueError(f"spatial sample is missing {coordinate}")
                elif numeric_values["reference_radius_m"] is None:
                    raise ValueError("radial/harmonic sample is missing reference radius")
                # Serialize the canonical typed identity compactly. Keeping a tuple
                # of row strings for every Phobos sample costs several hundred MB.
                key = json.dumps(_row_key(row), separators=(",", ":"),
                                 ensure_ascii=True).encode("ascii")
                if key in seen:
                    raise ValueError("duplicate sample key")
                seen.add(key)
                yield row
        finally:
            stream.close()

    return manifest, rows()
