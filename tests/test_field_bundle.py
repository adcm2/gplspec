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
    model_data = Path(sys.argv[2])
    with tempfile.TemporaryDirectory(prefix="gplspec-field-bundle-") as temp:
        root = Path(temp)
        bundle = root / "bundle"
        rotated_bundle = root / "rotated_bundle"
        surface_bundle = root / "surface_bundle"
        rotated_authoritative_bundle = root / "rotated_authoritative_bundle"
        unrotated_angular_bundle = root / "unrotated_angular_bundle"
        surface_slice_bundle = root / "surface_slice_bundle"
        full_geometry = root / "full_geometry"
        truncated_geometry = root / "truncated_geometry"
        full_slice_geometry = root / "full_slice_geometry"
        truncated_slice_geometry = root / "truncated_slice_geometry"
        surface_equatorial_slice_bundle = root / "surface_equatorial_slice_bundle"
        input_file = root / "input.txt"
        input_file.write_bytes(b"verified input bytes\n")
        (root / "auxiliary.txt").write_bytes(b"another runtime input\n")
        runner = ROOT / "work" / "run_field_bundle.py"
        subprocess.run([sys.executable, str(runner), "--workdir", str(root),
                        "--input", "input.txt", "--input", "auxiliary.txt",
                        "--bundle", "bundle", "--bundle", "rotated_bundle",
                        "--bundle", "surface_bundle", "--bundle", "rotated_authoritative_bundle",
                        "--bundle", "unrotated_angular_bundle",
                        "--bundle", "surface_slice_bundle",
                        "--bundle", "full_geometry", "--bundle", "truncated_geometry",
                        "--bundle", "full_slice_geometry", "--bundle", "truncated_slice_geometry",
                        "--bundle", "surface_equatorial_slice_bundle",
                        "--", str(executable),
                        str(bundle), str(rotated_bundle), str(surface_bundle),
                        str(rotated_authoritative_bundle), str(unrotated_angular_bundle),
                        str(surface_slice_bundle), str(full_geometry), str(truncated_geometry),
                        str(full_slice_geometry), str(truncated_slice_geometry),
                        str(model_data), str(surface_equatorial_slice_bundle)], check=True)
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
        authoritative_manifest, authoritative_rows = MODULE.open_bundle(
            rotated_authoritative_bundle)
        assert authoritative_manifest["rotation"]["angles_rad"] == [0.2, 0.7, -0.1]
        assert list(authoritative_rows)
        unrotated_manifest, unrotated_rows = MODULE.open_bundle(unrotated_angular_bundle)
        assert unrotated_manifest["rotation"]["angles_rad"] is None
        assert list(unrotated_rows)
        surface_manifest, surface_rows = MODULE.open_bundle(surface_bundle)
        surface_data = list(surface_rows)
        assert surface_manifest["schema_version"] == 1
        assert len(surface_data) == 1
        assert surface_data[0]["quantity"] == "surface_radius"
        assert surface_data[0]["radial_node"] == ""
        with (surface_bundle / "samples.csv").open(newline="", encoding="utf-8") as stream:
            serialized_surface = next(csv.DictReader(stream))
        assert serialized_surface["radial_node"] == ""

        surface_slice_manifest, surface_slice_rows = MODULE.open_bundle(surface_slice_bundle)
        surface_slice_data = list(surface_slice_rows)
        assert surface_slice_manifest["rotation"]["angles_rad"] == [0.2, 0.7, -0.1]
        assert surface_slice_data
        assert {row["representation"] for row in surface_slice_data} == {"rotated_slice"}
        assert {row["quantity"] for row in surface_slice_data} == {"surface_radius"}
        assert len(surface_slice_data) == 12
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        sys.path.insert(0, str(ROOT / "work"))
        import plot_fields
        from types import SimpleNamespace
        output = root / "surface-slice.png"
        plot_fields.slice_plot(SimpleNamespace(
            bundle=[str(surface_slice_bundle)], quantity="surface_radius",
            element=None, output=str(output)))
        assert output.is_file() and output.stat().st_size > 0
        rendered_lines = plt.gcf().axes[0].lines
        assert len(rendered_lines) == 2
        for line in rendered_lines:
            rendered_radii = np.hypot(line.get_xdata(), line.get_ydata())
            expected = float(surface_slice_data[0]["value_real"])
            np.testing.assert_allclose(rendered_radii, expected, rtol=1e-12, atol=1e-8)
        plt.close("all")
        equatorial_manifest, equatorial_rows = MODULE.open_bundle(surface_equatorial_slice_bundle)
        assert equatorial_manifest["rotation"]["angles_rad"] is None
        assert {row["representation"] for row in equatorial_rows} == {"rotated_slice"}
        plot_fields.slice_plot(SimpleNamespace(
            bundle=[str(surface_equatorial_slice_bundle)], quantity="surface_radius",
            element=None, output=None))
        assert "body-fixed frame" in plt.gcf().axes[0].get_xlabel()
        plt.close("all")
        try:
            plot_fields.slice_plot(SimpleNamespace(
                bundle=[str(surface_slice_bundle)], quantity="surface_radius",
                element=999, output=None))
        except ValueError as error:
            assert "element is absent" in str(error)
        else:
            raise AssertionError("missing surface slice element was accepted")

        def physical_radii(directory):
            _, samples = MODULE.open_bundle(directory)
            return [float(row["physical_radius_m"]) for row in samples]

        for full_path, truncated_path in ((full_geometry, truncated_geometry),
                                          (full_slice_geometry, truncated_slice_geometry)):
            full_radii = physical_radii(full_path)
            truncated_radii = physical_radii(truncated_path)
            assert len(full_radii) == len(truncated_radii) and full_radii
            assert max(abs(a - b) for a, b in zip(full_radii, truncated_radii)) > 1.0e-8

        example_source = (ROOT / "examples" / "phobos_gravity_perturbation.cpp").read_text()
        for product, writer in (("full_rotated_grid_bundle", "WriteAngular"),
                                ("full_rotated_slice_bundle", "WriteSlice")):
            import re
            call = re.search(rf"{writer}\(\s*std::filesystem::path\(pathtofolder3\) / \"{product}\"(.*?)\);",
                             example_source, re.S)
            assert call and re.search(r"\bphobosfull\s*,\s*full_result\.solution", call.group(1))

        legacy_surface = root / "legacy_surface_manifest"
        edit_manifest(surface_bundle, legacy_surface, lambda value:
                      value["quantity_contracts"]["surface_radius"]["representations"].pop("rotated_slice"))
        assert list(MODULE.open_bundle(legacy_surface)[1])
        undeclared_surface_slice = root / "undeclared_surface_slice"
        edit_manifest(surface_slice_bundle, undeclared_surface_slice, lambda value:
                      value["quantity_contracts"]["surface_radius"]["representations"].pop("rotated_slice"))
        assert rejected(undeclared_surface_slice)

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
