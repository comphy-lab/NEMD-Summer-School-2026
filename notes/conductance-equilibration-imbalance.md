# Conductance equilibration and heat-bath imbalance

This note documents a local investigation in the CoMPhy fork of the Day 1
`conductance` exercise. It is not a replacement for the official course handout;
it records an extra diagnostic question that came up while testing the example
on Cirrus.

## Question

The Day 1 handout asks students to break the conductance measurement by reducing
the equilibration length:

```bash
../submit.sh conductance.in -var nequil 10000
python analyse_conductance.py
```

The expected failure mode is correct: with too little equilibration, the two
wall baths do not exchange heat at equal and opposite rates, so the bath-power
imbalance is large and the inferred conductances are unreliable.

The counterintuitive observation was that increasing `nequil` beyond the default
did not always reduce the reported imbalance. We tested whether that is a real
monotonic-equilibration issue or a finite-sampling issue in the heat-current
slope estimate.

## Test

The helper scripts added in this fork are:

- `day1/conductance/submit_nequil_matrix.sh`
- `day1/conductance/analyse_conductance_nequil.py`

The matrix runner executes independent vanilla `conductance.in` jobs for a small
set of equilibration lengths and seeds:

```text
nequil = 10000, 40000, 80000, 160000
seed   = 90187, 90287
nprod  = 40000
```

Each run uses the same base conductance case: `Tbot = 1.3`, `Ttop = 0.9`,
`eps_wf = 0.6`. The analyser computes the same bath-power imbalance used by
`analyse_conductance.py`:

```text
imbalance = | |P_bot| - |P_top| | / max(|P_bot|, |P_top|)
```

where `P_bot` and `P_top` are the fitted slopes of the cumulative heat exchanged
by the two Langevin wall baths.

## Result from Cirrus run 283712

With two seeds and `nprod = 40000`, the imbalance was:

| nequil | seed 90187 | seed 90287 |
|---:|---:|---:|
| 10000 | 51.3% | 53.4% |
| 40000 | 1.4% | 4.0% |
| 80000 | 6.0% | 7.4% |
| 160000 | 8.2% | 4.6% |

The result confirms the handout's main point: `nequil = 10000` is far too short.
It also shows that the imbalance is not a monotonic function of `nequil` in a
single finite stochastic trajectory. The default `nequil = 40000` was best in
this small test, while longer equilibration shifted the fixed-length production
window to a different stochastic segment whose heat-current slope estimate was
noisier.

## Interpretation

The imbalance is not itself a thermodynamic state variable. It is a finite-window
estimate derived from two noisy cumulative-heat slopes over the production
interval. Increasing `nequil` changes where the production window starts; it
does not make the slope estimate over `nprod = 40000` intrinsically less noisy.

Therefore:

- too little equilibration gives a real transient and a large imbalance;
- after the run is approximately steady, additional equilibration can still make
  the reported imbalance go up or down in one seed;
- the correct diagnostic is statistical: compare multiple seeds and, ideally,
  longer production windows.

The next useful test is larger:

```text
nequil = 10000, 20000, 40000, 80000, 160000
nprod  = 40000 and 160000
seeds  = 16 or 32
```

If longer `nprod` collapses the scatter, the main culprit is finite-window
slope noise. If the median imbalance genuinely rises with `nequil` even across
many seeds and long `nprod`, then there is a deeper slow-relaxation or sampling
issue worth investigating.

## Local artefacts

The local result bundle for the small matrix is outside the Git repository:

```text
/Users/vatsal/cowork-os/0-Projects/NEMD-2026-local/results/day1-conductance-nequil-matrix/
```

It contains:

- `day1_nequil_summary.csv`
- `day1_nequil_imbalance.png`
- `day1_nequil_conductance.png`
- per-case `day1_nequil_*_heat.dat`, `*_Tz.profile`, `*_params.txt`, and logs.
