import hashlib
import importlib.util
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
import csv
import json

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("field_bundle", ROOT / "work" / "field_bundle.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def rejected(directory):
    try:
        _, rows = MODULE.open_bundle(directory)
        list(rows)
    except (ValueError, KeyError):
        return True
    return False


def copy_bundle(source, destination):
    shutil.copytree(source, destination)


def edit_manifest(source, destination, edit):
    copy_bundle(source, destination)
    path = destination / "manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    edit(manifest)
    path.write_text(json.dumps(manifest), encoding="utf-8")


def edit_first_row(source, destination, edit):
    copy_bundle(source, destination)
    path = destination / "samples.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        header = reader.fieldnames
    edit(rows[0])
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


def main():
    executable = Path(sys.argv[1])
    with tempfile.TemporaryDirectory(prefix="gplspec-field-bundle-") as temp:
        root = Path(temp)
        bundle = root / "bundle"
        rotated_bundle = root / "rotated_bundle"
        input_file = root / "input.txt"
        input_file.write_bytes(b"verified input bytes\n")
        (root / "auxiliary.txt").write_bytes(b"another runtime input\n")
        runner = ROOT / "work" / "run_field_bundle.py"
        subprocess.run([sys.executable, str(runner), "--workdir", str(root),
                        "--input", "input.txt", "--input", "auxiliary.txt",
                        "--bundle", "bundle", "--bundle", "rotated_bundle", "--",
                        str(executable), str(bundle), str(rotated_bundle)], check=True)
        manifest, rows = MODULE.open_bundle(bundle)
        data = list(rows)
        assert manifest["schema_version"] == 1
        assert manifest["harmonics"]["index"] == "l*l+l+m"
        assert manifest["quantity_contracts"]["density_physical"]["meaning"].endswith(
            "rho_physical=rho_referential/J.")
        assert manifest["rotation"]["convention"] == MODULE.EXPECTED_ROTATION
        assert manifest["rotation"]["frame"] == MODULE.EXPECTED_ROTATION_FRAME
        assert manifest["inputs"][0]["path"] == "input.txt"
        assert manifest["inputs"][0]["sha256"] == hashlib.sha256(
            b"verified input bytes\n").hexdigest()
        assert len(data) > 0
        assert {row["representation"] for row in data} == {
            "harmonic_coefficients", "radial_profile"}
        assert {row["quantity"] for row in data} == {"potential", "sensitivity_kernel"}
        assert {row["side"] for row in data if row["boundary"] == "body_vacuum"} == {"inner", "outer"}
        assert all(float(row["value_imag"]) == 0.0 for row in data)
        rotated_manifest, rotated_rows = MODULE.open_bundle(rotated_bundle)
        rotated_data = list(rotated_rows)
        assert rotated_manifest["rotation"]["angles_rad"] == [0.2, 0.7, -0.1]
        assert rotated_manifest["rotation"]["legacy_harmonic_coefficient_scale"] < 1.0
        assert rotated_data and {row["representation"] for row in rotated_data} == {"angular_grid"}
        assert {row["quantity"] for row in rotated_data} == {
            "potential", "density_referential", "density_physical"}
        grid_groups = {}
        for row in rotated_data:
            grid_groups.setdefault((row["quantity"], row["element"], row["radial_node"]), []).append(row)
        assert len(grid_groups) >= 12 and len({len(rows) for rows in grid_groups.values()}) == 1
        assert all(len({row["theta_rad"] for row in rows}) > 1 and
                   len({row["phi_rad"] for row in rows}) > 1
                   for rows in grid_groups.values())

        malformed = root / "malformed"
        copy_bundle(bundle, malformed)
        path = malformed / "samples.csv"
        path.write_text("wrong,header\n", encoding="utf-8")
        assert rejected(malformed)

        nonfinite = root / "nonfinite"
        copy_bundle(bundle, nonfinite)
        path = nonfinite / "samples.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        columns = lines[0].split(",")
        lines[1] = lines[1].split(",")[:-2] + ["nan", "0"]
        lines[1] = ",".join(lines[1])
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        assert rejected(nonfinite)

        duplicate = root / "duplicate"
        copy_bundle(bundle, duplicate)
        path = duplicate / "samples.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        lines.append(lines[1])
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        assert rejected(duplicate)

        canonical_duplicate = root / "canonical-duplicate"
        copy_bundle(bundle, canonical_duplicate)
        path = canonical_duplicate / "samples.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        fields = lines[1].split(",")
        radius_index = MODULE.HEADER.index("reference_radius_m")
        original = float(fields[radius_index])
        alternative = format(original, ".17e")
        if alternative == fields[radius_index]:
            alternative = format(original, ".17g")
        fields[radius_index] = alternative
        lines.append(",".join(fields))
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        assert rejected(canonical_duplicate)

        incompatible = root / "incompatible"
        edit_first_row(bundle, incompatible,
                       lambda row: row.update(quantity="power_std"))
        assert rejected(incompatible)

        bad_harmonic = root / "bad-harmonic"
        edit_first_row(bundle, bad_harmonic,
                       lambda row: row.update(order="99"))
        assert rejected(bad_harmonic)

        bad_scale = root / "bad-scale"
        edit_manifest(bundle, bad_scale,
                      lambda value: value["harmonics"].update(length_norm_m=float("inf")))
        assert rejected(bad_scale)

        bad_rotation = root / "bad-rotation"
        edit_manifest(bundle, bad_rotation,
                      lambda value: value["rotation"].update(angles_rad=[0.0, float("nan"), 0.0]))
        assert rejected(bad_rotation)

        bad_diagnostic = root / "bad-diagnostic"
        edit_manifest(bundle, bad_diagnostic,
                      lambda value: value.update(solve_diagnostics=[{
                          "status": 0, "iterations": 1, "tolerance": 1e-12,
                          "recursive_error": 0.0, "residual_norm": float("nan"),
                          "relative_residual": 0.0}]))
        assert rejected(bad_diagnostic)

    print(f"validated {len(data)} field-bundle rows and metadata, compatibility, range, canonical-key, malformed, nonfinite, and duplicate rejection")


if __name__ == "__main__":
    main()
