#!/usr/bin/env python3
"""Day 1 extension - wetting dependence of Kapitza conductance.

Reads the case files written by conductance_wetting.in and compares the
interfacial conductance G = J / dT as the wall-fluid LJ attraction eps_wf is
varied. Larger eps_wf is the simple Day-1 proxy for a more wetting wall.

Run after:

    ../submit.sh conductance_wetting.in
    python analyse_conductance_wetting.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lammps_io import read_params, read_profile, read_timeseries

SURF = 1.2


def linfit(x, y):
    """Least-squares slope, intercept, R^2."""
    if len(x) < 2 or np.allclose(x, x[0]):
        return np.nan, np.nan, np.nan
    s, b = np.polyfit(x, y, 1)
    yhat = s * x + b
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return s, b, r2


def discover_cases():
    """Return [(name, eps_wf), ...] from day1_wetting_cases.txt."""
    fn = "day1_wetting_cases.txt"
    if not os.path.exists(fn):
        sys.exit(
            "day1_wetting_cases.txt not found. Run `../submit.sh conductance_wetting.in` first."
        )
    cases = []
    with open(fn, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            name, eps, *_ = line.split()
            cases.append((name, float(eps)))
    if not cases:
        sys.exit("no wetting cases found in day1_wetting_cases.txt")
    return cases


def analyse_case(name, eps_wf):
    """Compute per-wall conductance for one wetting case."""
    profile = f"day1_wetting_{name}_Tz.profile"
    heat = f"day1_wetting_{name}_heat.dat"
    params = f"day1_wetting_{name}_params.txt"
    for fn in (profile, heat, params):
        if not os.path.exists(fn):
            sys.exit(f"{fn} not found; wait for the wetting job to finish.")

    P = read_params(params)
    Lz, wall_th, area = P["Lz"], P["wall_th"], P["area"]
    dt, Tbot, Ttop = P["dt"], P["Tbot"], P["Ttop"]
    z_wall_bot, z_wall_top = wall_th, Lz - wall_th

    z, N, T = read_profile(profile)
    m = (N > 0.5) & (T > 1e-6)
    z, N, T = z[m], N[m], T[m]

    bot_win = (z < z_wall_bot) & (z >= z_wall_bot - SURF)
    top_win = (z > z_wall_top) & (z <= z_wall_top + SURF)
    if not bot_win.any() or not top_win.any():
        sys.exit(f"{name}: no wall-surface bins found")
    Tw_bot = np.average(T[bot_win], weights=N[bot_win])
    Tw_top = np.average(T[top_win], weights=N[top_win])

    layer = 2.0
    fluid = (z > z_wall_bot + layer) & (z < z_wall_top - layer)
    if fluid.sum() < 3:
        sys.exit(f"{name}: fluid-interior window < 3 bins")
    sT, bT, r2T = linfit(z[fluid], T[fluid])
    Tf_bot = sT * z_wall_bot + bT
    Tf_top = sT * z_wall_top + bT
    dT_bot = abs(Tf_bot - Tw_bot)
    dT_top = abs(Tf_top - Tw_top)

    ts, q_bot, q_top = read_timeseries(heat)
    t = ts * dt
    p_bot, _, r2lo = linfit(t, q_bot)
    p_top, _, r2hi = linfit(t, q_top)
    J_bot, J_top = abs(p_bot) / area, abs(p_top) / area
    G_bot = J_bot / dT_bot if dT_bot > 1e-9 else float("nan")
    G_top = J_top / dT_top if dT_top > 1e-9 else float("nan")
    imbal = abs(abs(p_bot) - abs(p_top)) / max(abs(p_bot), abs(p_top), 1e-12)

    hot = "bottom" if Tbot > Ttop else "top"
    return {
        "case": name,
        "eps_wf": eps_wf,
        "hot": hot,
        "sT": sT,
        "r2T": r2T,
        "dT_bot": dT_bot,
        "dT_top": dT_top,
        "J_bot": J_bot,
        "J_top": J_top,
        "G_bot": G_bot,
        "G_top": G_top,
        "G_mean": 0.5 * (G_bot + G_top),
        "dT_mean": 0.5 * (dT_bot + dT_top),
        "imbalance": imbal,
        "r2_heat_bot": r2lo,
        "r2_heat_top": r2hi,
        "z": z,
        "N": N,
        "T": T,
        "z_wall_bot": z_wall_bot,
        "z_wall_top": z_wall_top,
        "Tw_bot": Tw_bot,
        "Tw_top": Tw_top,
        "Tf_bot": Tf_bot,
        "Tf_top": Tf_top,
        "fit_slope": sT,
        "fit_intercept": bT,
    }


def write_table(results, out="day1_wetting_summary.csv"):
    headers = [
        "case",
        "eps_wf",
        "G_bot",
        "G_top",
        "G_mean",
        "dT_bot",
        "dT_top",
        "dT_mean",
        "J_bot",
        "J_top",
        "imbalance",
        "r2T",
    ]
    with open(out, "w", encoding="utf-8") as handle:
        handle.write(",".join(headers) + "\n")
        for row in results:
            handle.write(
                ",".join(
                    [
                        row["case"],
                        f"{row['eps_wf']:.6g}",
                        f"{row['G_bot']:.6g}",
                        f"{row['G_top']:.6g}",
                        f"{row['G_mean']:.6g}",
                        f"{row['dT_bot']:.6g}",
                        f"{row['dT_top']:.6g}",
                        f"{row['dT_mean']:.6g}",
                        f"{row['J_bot']:.6g}",
                        f"{row['J_top']:.6g}",
                        f"{row['imbalance']:.6g}",
                        f"{row['r2T']:.6g}",
                    ]
                )
                + "\n"
            )
    print(f"    summary -> {out}")


def import_pyplot():
    try:
        import matplotlib

        if not os.environ.get("DISPLAY"):
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    return plt


def maybe_show(plt):
    if os.environ.get("DISPLAY"):
        try:
            plt.show()
        except KeyboardInterrupt:
            pass


def plot_conductance(results, out="day1_wetting_conductance.png"):
    plt = import_pyplot()
    if plt is None:
        print("    (matplotlib not found - skipping the conductance plot)")
        return
    eps = np.array([r["eps_wf"] for r in results])
    g_bot = np.array([r["G_bot"] for r in results])
    g_top = np.array([r["G_top"] for r in results])
    g_mean = np.array([r["G_mean"] for r in results])

    order = np.argsort(eps)
    eps, g_bot, g_top, g_mean = eps[order], g_bot[order], g_top[order], g_mean[order]

    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    ax.plot(eps, g_mean, "o-", color="#1F3A5F", label="mean")
    ax.plot(eps, g_bot, "s--", color="#B23A2E", alpha=0.8, label="bottom wall")
    ax.plot(eps, g_top, "^--", color="#3C7A3F", alpha=0.8, label="top wall")
    ax.set_xlabel(r"wall-fluid attraction $\epsilon_{wf}$")
    ax.set_ylabel(r"interfacial conductance $G=J/\Delta T$")
    ax.set_title("wetting dependence of Kapitza conductance")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    print(f"    plot -> {out}")
    maybe_show(plt)


def plot_temperature_profiles(results, out="day1_wetting_Tz_profiles.png"):
    plt = import_pyplot()
    if plt is None:
        print("    (matplotlib not found - skipping the T(z) plot)")
        return

    colors = ["#4B5563", "#1F3A5F", "#B23A2E", "#3C7A3F", "#7C3AED"]
    fig, ax = plt.subplots(figsize=(4.9, 4.2))
    for i, row in enumerate(sorted(results, key=lambda r: r["eps_wf"])):
        color = colors[i % len(colors)]
        z, T = row["z"], row["T"]
        label = rf"$\epsilon_{{wf}}={row['eps_wf']:.3g}$"
        ax.scatter(T, z, s=10, alpha=0.34, color=color)
        zf = np.linspace(row["z_wall_bot"], row["z_wall_top"], 80)
        ax.plot(
            row["fit_slope"] * zf + row["fit_intercept"],
            zf,
            color=color,
            lw=1.4,
            label=label,
        )
        for Tw, Tf, zw in (
            (row["Tw_bot"], row["Tf_bot"], row["z_wall_bot"]),
            (row["Tw_top"], row["Tf_top"], row["z_wall_top"]),
        ):
            ax.plot([Tw], [zw], "o", color=color, ms=4.8)
            ax.plot([Tf], [zw], "s", color=color, ms=4.4)
            ax.plot([Tw, Tf], [zw, zw], color=color, lw=0.9, alpha=0.8)

    ax.axhline(results[0]["z_wall_bot"], color="#999", ls=":", lw=0.8)
    ax.axhline(results[0]["z_wall_top"], color="#999", ls=":", lw=0.8)
    ax.set_xlabel(r"$T(z)$")
    ax.set_ylabel(r"$z\ (\sigma)$")
    ax.set_title("temperature profiles for wetting sweep")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    print(f"    T(z) plot -> {out}")
    maybe_show(plt)


def plot_temperature_jumps(results, out="day1_wetting_temperature_jump.png"):
    plt = import_pyplot()
    if plt is None:
        print("    (matplotlib not found - skipping the temperature-jump plot)")
        return

    rows = sorted(results, key=lambda r: r["eps_wf"])
    eps = np.array([r["eps_wf"] for r in rows])
    d_bot = np.array([r["dT_bot"] for r in rows])
    d_top = np.array([r["dT_top"] for r in rows])
    d_mean = np.array([r["dT_mean"] for r in rows])

    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    ax.plot(eps, d_mean, "o-", color="#1F3A5F", label="mean")
    ax.plot(eps, d_bot, "s--", color="#B23A2E", alpha=0.8, label="bottom wall")
    ax.plot(eps, d_top, "^--", color="#3C7A3F", alpha=0.8, label="top wall")
    ax.set_xlabel(r"wall-fluid attraction $\epsilon_{wf}$")
    ax.set_ylabel(r"temperature jump $\Delta T$")
    ax.set_title("Kapitza jump vs wetting strength")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    print(f"    jump plot -> {out}")
    maybe_show(plt)


def main():
    results = [analyse_case(name, eps) for name, eps in discover_cases()]
    results.sort(key=lambda row: row["eps_wf"])

    print("\nDay 1 extension: wetting dependence of interfacial conductance")
    print("    eps_wf is the wall-fluid LJ attraction; larger eps_wf = more wetting.")
    for row in results:
        print(
            f"    {row['case']} eps_wf={row['eps_wf']:.3g}: "
            f"G_bot={row['G_bot']:.3f}  G_top={row['G_top']:.3f}  "
            f"G_mean={row['G_mean']:.3f}  imbalance={100*row['imbalance']:.1f}%  "
            f"T-fit R^2={row['r2T']:.3f}"
        )
        if row["imbalance"] > 0.05:
            print("      WARNING: not steady enough; lengthen nequil/nprod.")
        if row["r2T"] < 0.95:
            print("      WARNING: poor linear T(z) fit; G is unreliable.")

    write_table(results)
    plot_conductance(results)
    plot_temperature_profiles(results)
    plot_temperature_jumps(results)


if __name__ == "__main__":
    main()
