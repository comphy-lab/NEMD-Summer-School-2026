#!/usr/bin/env python3
"""Analyse vanilla conductance equilibration-length matrix.

This reads the files produced by submit_nequil_matrix.sh and asks whether the
heat-bath imbalance decreases monotonically with nequil, or only on average
across stochastic runs.
"""
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lammps_io import read_params, read_profile, read_timeseries

SURF = 1.2


def linfit(x, y):
    if len(x) < 2 or np.allclose(x, x[0]):
        return np.nan, np.nan, np.nan
    s, b = np.polyfit(x, y, 1)
    yhat = s * x + b
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return s, b, r2


def read_cases(path="day1_nequil_cases.txt"):
    if not os.path.exists(path):
        sys.exit("day1_nequil_cases.txt not found. Run ./submit_nequil_matrix.sh first.")
    cases = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            case, nequil, seed, nprod = line.split()
            cases.append(
                {
                    "case": case,
                    "nequil": int(nequil),
                    "seed": int(seed),
                    "nprod": int(nprod),
                }
            )
    return cases


def analyse_case(meta):
    case = meta["case"]
    profile = f"day1_nequil_{case}_Tz.profile"
    heat = f"day1_nequil_{case}_heat.dat"
    params = f"day1_nequil_{case}_params.txt"
    for fn in (profile, heat, params):
        if not os.path.exists(fn):
            sys.exit(f"{fn} not found; wait for the matrix job to finish.")

    P = read_params(params)
    Lz, wall_th, area = P["Lz"], P["wall_th"], P["area"]
    dt = P["dt"]
    z_wall_bot, z_wall_top = wall_th, Lz - wall_th

    z, N, T = read_profile(profile)
    m = (N > 0.5) & (T > 1e-6)
    z, N, T = z[m], N[m], T[m]

    bot_win = (z < z_wall_bot) & (z >= z_wall_bot - SURF)
    top_win = (z > z_wall_top) & (z <= z_wall_top + SURF)
    Tw_bot = np.average(T[bot_win], weights=N[bot_win])
    Tw_top = np.average(T[top_win], weights=N[top_win])

    fluid = (z > z_wall_bot + 2.0) & (z < z_wall_top - 2.0)
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
    imbal = abs(abs(p_bot) - abs(p_top)) / max(abs(p_bot), abs(p_top), 1e-12)

    return {
        **meta,
        "imbalance": imbal,
        "r2T": r2T,
        "r2_heat_bot": r2lo,
        "r2_heat_top": r2hi,
        "G_bot": J_bot / dT_bot,
        "G_top": J_top / dT_top,
        "G_mean": 0.5 * (J_bot / dT_bot + J_top / dT_top),
        "dT_bot": dT_bot,
        "dT_top": dT_top,
        "J_bot": J_bot,
        "J_top": J_top,
    }


def write_summary(rows, out="day1_nequil_summary.csv"):
    headers = [
        "case",
        "nequil",
        "seed",
        "nprod",
        "imbalance",
        "G_bot",
        "G_top",
        "G_mean",
        "dT_bot",
        "dT_top",
        "J_bot",
        "J_top",
        "r2T",
    ]
    with open(out, "w", encoding="utf-8") as handle:
        handle.write(",".join(headers) + "\n")
        for row in rows:
            handle.write(
                ",".join(
                    [
                        row["case"],
                        str(row["nequil"]),
                        str(row["seed"]),
                        str(row["nprod"]),
                        f"{row['imbalance']:.8g}",
                        f"{row['G_bot']:.8g}",
                        f"{row['G_top']:.8g}",
                        f"{row['G_mean']:.8g}",
                        f"{row['dT_bot']:.8g}",
                        f"{row['dT_top']:.8g}",
                        f"{row['J_bot']:.8g}",
                        f"{row['J_top']:.8g}",
                        f"{row['r2T']:.8g}",
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


def plot_imbalance(rows, out="day1_nequil_imbalance.png"):
    plt = import_pyplot()
    if plt is None:
        print("    (matplotlib not found - skipping imbalance plot)")
        return

    by_seed = defaultdict(list)
    by_nequil = defaultdict(list)
    for row in rows:
        by_seed[row["seed"]].append(row)
        by_nequil[row["nequil"]].append(row["imbalance"])

    fig, ax = plt.subplots(figsize=(5.2, 3.7))
    for seed, seed_rows in sorted(by_seed.items()):
        seed_rows = sorted(seed_rows, key=lambda r: r["nequil"])
        ax.plot(
            [r["nequil"] for r in seed_rows],
            [100 * r["imbalance"] for r in seed_rows],
            "o--",
            alpha=0.55,
            label=f"seed {seed}",
        )

    nq = np.array(sorted(by_nequil))
    means = np.array([100 * np.mean(by_nequil[n]) for n in nq])
    stds = np.array([100 * np.std(by_nequil[n], ddof=1) if len(by_nequil[n]) > 1 else 0 for n in nq])
    ax.errorbar(nq, means, yerr=stds, fmt="o-", color="#111827", lw=2.0, label="mean")
    ax.axhline(5.0, color="#B23A2E", ls=":", lw=1.0, label="5% warning threshold")
    ax.set_xscale("log")
    ax.set_xlabel("equilibration steps, nequil")
    ax.set_ylabel("heat-bath imbalance (%)")
    ax.set_title("vanilla conductance: imbalance vs equilibration")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=170)
    print(f"    imbalance plot -> {out}")


def plot_conductance(rows, out="day1_nequil_conductance.png"):
    plt = import_pyplot()
    if plt is None:
        print("    (matplotlib not found - skipping conductance plot)")
        return

    by_seed = defaultdict(list)
    for row in rows:
        by_seed[row["seed"]].append(row)

    fig, ax = plt.subplots(figsize=(5.2, 3.7))
    for seed, seed_rows in sorted(by_seed.items()):
        seed_rows = sorted(seed_rows, key=lambda r: r["nequil"])
        ax.plot(
            [r["nequil"] for r in seed_rows],
            [r["G_mean"] for r in seed_rows],
            "o--",
            alpha=0.7,
            label=f"seed {seed}",
        )
    ax.set_xscale("log")
    ax.set_xlabel("equilibration steps, nequil")
    ax.set_ylabel(r"mean conductance $G$")
    ax.set_title("vanilla conductance estimate vs equilibration")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=170)
    print(f"    conductance plot -> {out}")


def main():
    rows = [analyse_case(meta) for meta in read_cases()]
    rows.sort(key=lambda r: (r["nequil"], r["seed"]))

    print("\nVanilla conductance equilibration matrix")
    for row in rows:
        print(
            f"    nequil={row['nequil']:6d} seed={row['seed']}: "
            f"imbalance={100*row['imbalance']:5.1f}%  "
            f"G_mean={row['G_mean']:.3f}  T-fit R^2={row['r2T']:.3f}"
        )

    write_summary(rows)
    plot_imbalance(rows)
    plot_conductance(rows)


if __name__ == "__main__":
    main()
