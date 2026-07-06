"""Shared output parsers for the Day-2 slit-pore analyser.

`read_vx_rho` reads the `poiseuille_w<width>.profile` files written by
`poiseuille.in`. Each profile row contains two binned values (vx and mass
density) from a single `fix ave/chunk`, so `read_vx_rho` returns both.
`read_params` reads `poiseuille_w<width>.params.txt`. The physics (the
momentum balance, the strain rate, eta(z)) is in `analyse_poiseuille.py`.

The analyser runs from inside `poiseuille/` (`cd poiseuille &&
python analyse_poiseuille.py`) and imports this file from the parent
directory:

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from lammps_io import read_vx_rho, read_params
"""
import getpass
import os
import sys
import tempfile

import numpy as np

# Point matplotlib's cache at a writable dir BEFORE any analyser imports matplotlib,
# so matplotlib does not warn and rebuild the font cache (~10 s) each run when $HOME
# is not writable on HPC.
# Per-user path: on a shared login node a fixed /tmp name is owned by whoever ran first,
# and everyone else gets a not-writable warning on every run.
os.environ.setdefault("MPLCONFIGDIR",
                      os.path.join(tempfile.gettempdir(), f"day2_mplcache_{getpass.getuser()}"))


def read_vx_rho(path):
    """LAMMPS `fix ave/chunk ... vx density/mass` file -> (z, Ncount, vx, rho_m)
    of the LAST complete block.

    Each block is a 2+-token header (timestep n_chunks [total_count]; only
    n_chunks is used) followed by n_chunks rows: chunk Coord1 Ncount vx density.
    Truncated/short trailing rows (a run still writing or killed mid-write) are
    tolerated by falling back to the last COMPLETE block, exactly as Day-1's
    read_profile does."""
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#") and l.strip()]
    blocks, i = [], 0
    while i < len(lines):
        try:
            n = int(lines[i].split()[1])                                # header: timestep n_chunks [total]
        except (ValueError, IndexError):
            break                                                       # garbled/partial header -> stop
        if n <= 0 or i + 1 + n > len(lines):
            break                                                       # trailing block not fully written
        rows, ok = [], True
        for k in range(n):
            c = lines[i + 1 + k].split()
            if len(c) < 5:                                              # row truncated before the density column
                ok = False
                break
            try:
                rows.append([float(x) for x in c[:5]])
            except ValueError:                                          # a value truncated mid-number
                ok = False
                break
        if not ok:
            break                                                       # incomplete trailing block -> use last complete
        blocks.append(np.array(rows))
        i += 1 + n
    if not blocks:
        sys.exit(f"{path}: no complete data block - the run may still be writing, was killed "
                 f"mid-write, or is too short. Wait for the job to finish (squeue --me) or "
                 f"delete {path} and rerun.")
    a = blocks[-1]
    return a[:, 1], a[:, 2], a[:, 3], a[:, 4]                           # Coord1, Ncount, vx, density/mass


def read_params(path):
    """Read a `key value` text file (poiseuille_w<width>.params.txt) -> {key: float}.

    Lines whose second token is non-numeric are skipped, so header/comment lines
    are ignored. Identical to Day-1's read_params."""
    p = {}
    with open(path) as f:
        for l in f:
            t = l.split()
            if len(t) >= 2:
                try:
                    p[t[0]] = float(t[1])
                except ValueError:
                    pass
    return p
