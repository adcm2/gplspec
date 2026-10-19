"""Plot GPLSpec field bundles.

Commands:
  radial       compare radial profiles on matching rays, split at element traces
  map          plot one angular-grid quantity on its stored angular grid
  slice        plot one explicitly evaluated rotated slice
  power        compare normalized degree power at selected radial samples
  topography   plot physical surface radius from TABLEA1 samples
  sweep        plot the CleanBench 7 centre-potential diagnostic table
"""
import argparse
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from field_bundle import open_bundle


PLOT_COLUMNS = (
    "representation", "element", "side", "material", "degree", "order",
    "radial_node", "reference_radius_m", "physical_radius_m", "theta_rad",
    "phi_rad", "value_real", "value_imag",
)


def load(directory, quantity):
    """Validate each row and retain only typed columns required by plotting."""
    manifest, rows = open_bundle(directory)
    columns = {name: [] for name in PLOT_COLUMNS}
    for row in rows:
        if row["quantity"] != quantity:
            continue
        for name in PLOT_COLUMNS:
            raw = row[name]
            if name in {"element", "material", "degree", "order", "radial_node"}:
                columns[name].append(-1 if raw == "" else int(raw))
            elif name in {"reference_radius_m", "physical_radius_m", "theta_rad", "phi_rad", "value_real", "value_imag"}:
                columns[name].append(np.nan if raw == "" else float(raw))
            else:
                columns[name].append(raw)
    if not columns["representation"]:
        raise ValueError(f"{directory}: no {quantity} samples")
    for name in columns:
        dtype = str if name in {"representation", "side"} else None
        columns[name] = np.asarray(columns[name], dtype=dtype)
    return manifest, columns


def selected(columns, representation):
    return columns["representation"] == representation


def radial(args):
    compare = args.difference and len(args.bundle) > 1
    fig, ax = plt.subplots(2 if compare else 1, 1, squeeze=False,
                           figsize=(9, 7 if compare else 5), sharex=True)
    axes = ax[:, 0]
    datasets = []
    for directory in args.bundle:
        manifest, c = load(directory, args.quantity)
        datasets.append((manifest, c))
    if compare:
        frame = (datasets[0][0]["coordinates"], datasets[0][0]["rotation"])
        if any((manifest["coordinates"], manifest["rotation"])
               != frame for manifest, _ in datasets[1:]):
            raise ValueError(
                "difference requires matching coordinate and rotation frames")
    available = [set(c["representation"]) for _, c in datasets]
    if args.representation != "auto":
        reps = [args.representation] * len(datasets)
        if any(rep not in choices for rep, choices in zip(reps, available)):
            raise ValueError("requested radial representation is absent from a bundle")
    elif all(("angular_grid" in choices or "radial_profile" in choices)
             for choices in available):
        reps = ["angular_grid" if "angular_grid" in choices else "radial_profile"
                for choices in available]
    elif all("harmonic_coefficients" in choices for choices in available):
        reps = ["harmonic_coefficients"] * len(datasets)
    else:
        raise ValueError("bundles have no compatible radial representation; compare separately")
    ray = None
    if any(rep in {"angular_grid", "radial_profile"} for rep in reps):
        ray_sets = []
        for i, rep in enumerate(reps):
            if rep not in {"angular_grid", "radial_profile"}:
                continue
            _, c = datasets[i]
            mask = selected(c, rep)
            angle_data = np.isfinite(c["theta_rad"][mask]) & np.isfinite(c["phi_rad"][mask])
            if np.any(angle_data):
                ray_sets.append(set(zip(c["theta_rad"][mask][angle_data],
                                        c["phi_rad"][mask][angle_data])))
        if ray_sets:
            common_rays = set.intersection(*ray_sets)
            if not common_rays:
                raise ValueError("bundles have no exactly matching sampled ray")
            requested = args.theta, args.phi
            ray = min(common_rays, key=lambda item: (item[0]-requested[0])**2 + (item[1]-requested[1])**2)
    reference = None
    component_column = "value_real" if args.component == "real" else "value_imag"
    for dataset_index, (manifest, c) in enumerate(datasets):
        rep = reps[dataset_index]
        mask = selected(c, rep)
        if rep in {"angular_grid", "radial_profile"} and ray is not None:
            has_angles = np.isfinite(c["theta_rad"]) & np.isfinite(c["phi_rad"])
            mask &= ~has_angles | ((c["theta_rad"] == ray[0]) &
                                  (c["phi_rad"] == ray[1]))
        elif rep == "harmonic_coefficients":
            mask &= (c["degree"] == args.degree) & (c["order"] == args.order)
        else:
            mask &= ((c["degree"] == -1) |
                     ((c["degree"] == args.degree) & (c["order"] == args.order)))
        groups = defaultdict(list)
        for i in np.flatnonzero(mask):
            # Both limits remain in one element segment; adjacent elements stay separate.
            groups[(int(c["element"][i]), int(c["material"][i]))].append(i)
        current = {}
        color = axes[0]._get_lines.get_next_color()
        first_group = True
        for key, indices in groups.items():
            indices = np.asarray(indices)
            order = np.lexsort((c["side"][indices], c["reference_radius_m"][indices]))
            indices = indices[order]
            # Harmonic coefficients and some exact radial profiles are defined
            # on the reference sphere; use physical radius only when supplied.
            physical = c["physical_radius_m"][indices]
            x = np.where(np.isfinite(physical), physical,
                         c["reference_radius_m"][indices])
            y = c[component_column][indices]
            current[key] = (c["reference_radius_m"][indices], x, y,
                            c["side"][indices])
            axes[0].plot(x, y, marker=".", color=color,
                         label=manifest["dataset"] if first_group else "_nolegend_")
            first_group = False
        if reference is not None and compare:
            if set(reference) != set(current):
                raise ValueError("difference requested but element/material trace sets differ")
            for key in sorted(set(reference), key=str):
                rr, rx, ry, rs = reference[key]
                cr, x, y, cs = current[key]
                same = (np.array_equal(rr, cr) and np.array_equal(rs, cs) and
                        np.allclose(rx, x, rtol=1e-12, atol=1e-14))
                if not same:
                    raise ValueError("difference requested for incompatible radial coordinates/traces")
                axes[1].plot(x, np.abs(y-ry), marker=".",
                             label=f"{manifest['dataset']} element {key[0]} material {key[1]}")
        if reference is None:
            reference = current
    units = datasets[0][0]["quantity_contracts"][args.quantity]["unit"]
    if all(rep == "harmonic_coefficients" for rep in reps):
        axes[0].set_ylabel(
            f"{args.quantity} {args.component} coefficient "
            f"(l={args.degree}, m={args.order}) ({units})")
        radius_label = "Reference radius (m)"
    else:
        axes[0].set_ylabel(f"{args.quantity} {args.component} component ({units})")
        if ray is not None:
            axes[0].set_title(f"Ray θ={ray[0]:.6g} rad, φ={ray[1]:.6g} rad")
        radius_label = "Radius (physical where available, otherwise reference; m)"
    axes[0].legend(fontsize=7)
    if compare:
        axes[1].set_ylabel("Absolute difference")
        axes[1].set_xlabel(radius_label)
    else:
        axes[0].set_xlabel(radius_label)
    fig.tight_layout()
    finish(args.output)


def _spatial(columns, representation, element, side):
    mask = selected(columns, representation)
    mask &= (columns["element"] == element) & (columns["side"] == side)
    indices = np.flatnonzero(mask)
    if not len(indices):
        raise ValueError(f"requested {representation} element trace is absent")
    return indices


def _frame_label(manifest):
    return "rotated output frame" if manifest["rotation"]["angles_rad"] is not None else "body-fixed frame"


def map_plot(args):
    fig, axes = plt.subplots(1, len(args.bundle), squeeze=False,
                             figsize=(7*len(args.bundle), 5))
    for directory, ax in zip(args.bundle, axes[0]):
        manifest, c = load(directory, args.quantity)
        indices = _spatial(c, "angular_grid", args.element, args.side)
        theta, phi = c["theta_rad"][indices], c["phi_rad"][indices]
        vals = c["value_real"][indices]
        image = ax.tricontourf(np.degrees(phi), 90.0-np.degrees(theta), vals, levels=100)
        fig.colorbar(image, ax=ax, label=f"{args.quantity} ({manifest['quantity_contracts'][args.quantity]['unit']})")
        frame = _frame_label(manifest)
        ax.set(xlabel=f"Longitude in {frame} (degrees)", ylabel=f"Latitude in {frame} (degrees)",
               title=f"{manifest['dataset']} element {args.element} {args.side} trace")
    finish(args.output)


def _slice_mesh(columns, elements):
    points_x, points_y, point_values, cells = [], [], [], []
    offset = 0
    for element in elements:
        indices = np.flatnonzero(selected(columns, "rotated_slice") &
                                  (columns["element"] == element) &
                                  np.isin(columns["side"], ["none", "outer", "inner"]))
        sides = columns["side"][indices]
        lower_side = "none" if element == 0 else "outer"
        lower = indices[sides == lower_side]
        upper = indices[sides == "inner"]
        if not len(upper):
            upper = indices[sides == "none"]
        if not len(lower) or not len(upper):
            raise ValueError(f"element {element} lacks both slice boundary traces")
        lower = lower[np.argsort(columns["phi_rad"][lower])]
        upper = upper[np.argsort(columns["phi_rad"][upper])]
        if len(lower) != len(upper):
            raise ValueError(f"element {element} boundary rings have different angular sampling")
        n = len(lower)
        def xy(idx):
            theta = columns["theta_rad"][idx]
            phi = columns["phi_rad"][idx]
            radius = columns["physical_radius_m"][idx]
            return radius*np.sin(theta)*np.cos(phi), radius*np.sin(theta)*np.sin(phi)
        xl, yl = xy(lower)
        xu, yu = xy(upper)
        if element == 0:
            x = np.r_[float(np.mean(xl)), xu]
            y = np.r_[float(np.mean(yl)), yu]
            values = np.r_[float(np.mean(columns["value_real"][lower])),
                           columns["value_real"][upper]]
            triangles = [[0, 1+i, 1+(i+1) % n] for i in range(n)]
        else:
            x, y = np.r_[xl, xu], np.r_[yl, yu]
            values = np.r_[columns["value_real"][lower],
                           columns["value_real"][upper]]
            triangles = []
            for i in range(n):
                j = (i+1) % n
                triangles.extend(([i, j, n+i], [j, n+j, n+i]))
        points_x.extend(x)
        points_y.extend(y)
        point_values.extend(values)
        cells.extend([[offset + vertex for vertex in triangle]
                      for triangle in triangles])
        offset += len(x)
    import matplotlib.tri as tri
    return tri.Triangulation(np.asarray(points_x), np.asarray(points_y),
                             triangles=np.asarray(cells)), np.asarray(point_values)


def slice_plot(args):
    fig, axes = plt.subplots(1, len(args.bundle), squeeze=False,
                             figsize=(7*len(args.bundle), 7))
    for directory, ax in zip(args.bundle, axes[0]):
        manifest, c = load(directory, args.quantity)
        if args.quantity == "surface_radius":
            indices = np.flatnonzero(selected(c, "rotated_slice"))
            if args.element is not None:
                indices = indices[c["element"][indices] == args.element]
            if not len(indices):
                raise ValueError("requested surface-radius slice element is absent")
            for element in np.unique(c["element"][indices]):
                ring = indices[c["element"][indices] == element]
                ring = ring[np.argsort(c["phi_rad"][ring])]
                radius = c["value_real"][ring]
                theta, phi = c["theta_rad"][ring], c["phi_rad"][ring]
                x = radius * np.sin(theta) * np.cos(phi)
                y = radius * np.sin(theta) * np.sin(phi)
                ax.plot(np.r_[x, x[0]], np.r_[y, y[0]], marker=".",
                        label=f"element {element} surface")
            ax.set_aspect("equal")
            frame = _frame_label(manifest)
            ax.set(xlabel=f"x in {frame} (m)", ylabel=f"y in {frame} (m)",
                   title=f"{manifest['dataset']} surface-radius slice")
            ax.legend(fontsize=7)
            continue
        available = np.unique(c["element"][selected(c, "rotated_slice")])
        elements = available if args.element is None else np.array([args.element])
        mesh, values = _slice_mesh(c, elements)
        image = ax.tripcolor(mesh, values, shading="gouraud")
        fig.colorbar(image, ax=ax, label=f"{args.quantity} ({manifest['quantity_contracts'][args.quantity]['unit']})")
        ax.set_aspect("equal")
        label = "all elements" if args.element is None else f"element {args.element}"
        ax.set(xlabel="Output-frame x (m)", ylabel="Output-frame y (m)",
               title=f"{manifest['dataset']} {label}")
    finish(args.output)

def power(args):
    fig, ax = plt.subplots()
    for directory in args.bundle:
        manifest, c = load(directory, "power_std")
        mask = selected(c, "degree_power") & (c["element"] == args.element)
        mask &= c["radial_node"] == args.radial_node
        indices = np.flatnonzero(mask)
        if not len(indices):
            raise ValueError(f"{directory}: selected power trace is absent")
        indices = indices[np.argsort(c["degree"][indices])]
        ax.plot(c["degree"][indices], c["value_real"][indices], marker=".", label=manifest["dataset"])
    ax.set_yscale("log")
    ax.set(xlabel="Harmonic degree l", ylabel="Normalized power (1)",
           title=f"Element {args.element}, local radial node {args.radial_node}")
    ax.legend()
    finish(args.output)


def topography(args):
    manifest, c = load(args.bundle, "surface_radius")
    indices = np.flatnonzero(selected(c, "angular_grid"))
    if not len(indices):
        raise ValueError("topography needs angular-grid samples")
    theta, phi = c["theta_rad"][indices], c["phi_rad"][indices]
    radius = c["value_real"][indices]
    fig, ax = plt.subplots(figsize=(8, 5))
    image = ax.tricontourf(np.degrees(phi), 90-np.degrees(theta), radius, levels=100)
    fig.colorbar(image, ax=ax, label="Physical surface radius (m)")
    frame = _frame_label(manifest)
    ax.set(xlabel=f"Longitude in {frame} (degrees)", ylabel=f"Latitude in {frame} (degrees)",
           title=f"{manifest['dataset']} surface")
    finish(args.output)

def load_sweep_table(path):
    """Load the retained three-column CleanBench 7 diagnostic table."""
    try:
        values = np.loadtxt(path, delimiter=";", ndmin=2, dtype=float)
    except (OSError, ValueError) as exc:
        raise ValueError(f"cannot read numeric semicolon table {path}: {exc}") from exc
    if values.ndim != 2 or values.shape[1] != 3 or values.shape[0] < 4:
        raise ValueError("sweep table must have at least four rows and exactly three columns")
    if not np.all(np.isfinite(values)):
        raise ValueError("sweep table contains a non-finite value")
    h, approximate, exact = values.T
    if np.any(h < 0) or np.any(np.diff(h) <= 0):
        raise ValueError("sweep amplitudes must be nonnegative and strictly increasing")
    normalization = float(np.max(np.abs(approximate)))
    if not np.isfinite(normalization) or normalization <= 0:
        raise ValueError("mapping-perturbation centre potentials must have nonzero scale")
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        relative_error = np.abs(exact - approximate) / normalization
    if not np.all(np.isfinite(relative_error)):
        raise ValueError("sweep relative errors are not finite")
    anchor = 3
    if h[anchor] <= 0 or relative_error[anchor] <= 0:
        raise ValueError("fourth sample must have positive h and nonzero relative error for the quadratic guide")
    with np.errstate(over="ignore", invalid="ignore"):
        quadratic_prediction = relative_error[anchor] * (h / h[anchor]) ** 2
    if not np.all(np.isfinite(quadratic_prediction)):
        raise ValueError("quadratic guide values are not finite")
    return h, approximate, exact, relative_error, quadratic_prediction


def sweep_plot(args):
    h, approximate, exact, relative_error, quadratic_prediction = load_sweep_table(args.table)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    axes[0].plot(h, exact, "b-", linewidth=2, label="Analytic exact")
    axes[0].plot(h, approximate, "r--", linewidth=2,
                 label="Mapping-perturbation approximation")
    axes[0].set(xlabel="Mapping amplitude h",
                ylabel="Centre potential (m² s⁻²)",
                title="CleanBench 7 centre potential")
    axes[0].legend(fontsize=8)

    # h=0 belongs in the potential panel; zero and zero-error values cannot
    # appear on logarithmic axes, so mask them explicitly in the error panel.
    valid_error = (h > 0) & (relative_error > 0)
    valid_prediction = (h > 0) & (quadratic_prediction > 0)
    axes[1].loglog(h[valid_error], relative_error[valid_error] * 100.0,
                   "b.", label="Relative error")
    axes[1].loglog(h[valid_prediction], quadratic_prediction[valid_prediction] * 100.0,
                   "r--", label="Quadratic guide (fourth sample)")
    axes[1].set(xlabel="Mapping amplitude h", ylabel="Relative error (%)",
                title="Relative error and quadratic guide")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    finish(args.output)


def finish(output):
    if output:
        plt.savefig(output, dpi=180, bbox_inches="tight")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("radial", help="compare radial profile bundles")
    p.add_argument("bundle", nargs="+")
    p.add_argument("--quantity", default="potential")
    p.add_argument("--representation", choices=("auto", "angular_grid", "harmonic_coefficients", "radial_profile"), default="auto")
    p.add_argument("--degree", type=int, default=0)
    p.add_argument("--order", type=int, default=0)
    p.add_argument("--component", choices=("real", "imag"), default="real",
                   help="plot the real or imaginary part of the complex field value")
    p.add_argument("--theta", type=float, default=0.0)
    p.add_argument("--phi", type=float, default=0.0)
    p.add_argument("--difference", action="store_true",
                   help="add absolute differences after strict trace/coordinate matching")
    p.add_argument("--output")
    p.set_defaults(run=radial)
    for name, function in (("map", map_plot), ("slice", slice_plot)):
        p = commands.add_parser(name)
        p.add_argument("bundle", nargs="+")
        p.add_argument("--quantity", default="potential")
        p.add_argument("--element", type=int, required=(name == "map"), default=None)
        if name == "map":
            p.add_argument("--side", choices=("none", "inner", "outer"), required=True)
        p.add_argument("--output")
        p.set_defaults(run=function)
    p = commands.add_parser("power")
    p.add_argument("bundle", nargs="+")
    p.add_argument("--element", type=int, required=True)
    p.add_argument("--radial-node", type=int, required=True)
    p.add_argument("--output")
    p.set_defaults(run=power)
    p = commands.add_parser("topography")
    p.add_argument("bundle")
    p.add_argument("--output")
    p.set_defaults(run=topography)
    p = commands.add_parser("sweep", help="plot the CleanBench 7 centre-potential sweep table")
    p.add_argument("table", help="three-column semicolon table: h;approximate;exact")
    p.add_argument("--output")
    p.set_defaults(run=sweep_plot)
    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
