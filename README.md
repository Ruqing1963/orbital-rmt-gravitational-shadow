# The Gravitational Shadow: Universal Level Repulsion in Planetary Orbital Architecture

**A Random Matrix Theory test from the Solar System to 320 Kepler systems**

Ruqing Chen, GUT Geoservice Inc., Montreal, Canada · June 2026

---

This is the **eighth** study in a unified Random Matrix Theory (RMT) program testing
whether single-source, charge-and-release systems show Wigner–Dyson **level repulsion**,
while superposed independent sources collapse to Poisson. It is the first to extend the
program from **time to space** and from **Earth to the cosmos**: the "events" are planets,
and the spacing is the logarithmic orbital separation Δln(a).

## The gravitational shadow hypothesis

When a planet forms it clears its feeding zone — its Hill sphere and resonance-swept
neighborhood — casting a **gravitational shadow** from which other bodies are excluded.
This cleared annulus is a hard minimum spacing, the spatial analog of a recharge time.
If the charge-and-release law generalizes, orbital architecture should show level
repulsion in Δln(a).

## Key results

| Target | n | ⟨r⟩ | 95% CI | CV | Class |
|---|---|---|---|---|---|
| Solar System (8 planets, unfolded) | 7 | **0.632** | [0.55, 0.87] | 0.34 | GOE–GUE |
| Galilean moons (4, unfolded) | 3 | (0.911) | — | **0.09** | ⚠ resonance-locked, excluded |
| **Kepler 320 systems (aligned pool)** | **810** | **0.674** | **[0.68, 0.72]** | 0.34 | **GUE** |
| Kepler global mix (frames scrambled) | 827 | 0.437 | [0.41, 0.45] | 0.86 | Poisson |
| Poisson null (independent clocks) | 810 | 0.396 | [0.36, 0.41] | 1.05 | Poisson |

Reference values: Poisson 0.386, GOE 0.531, GUE 0.603.

### What the numbers say

1. **Single gravitational wells repel.** The Solar System, after removing the
   Titius–Bode geometric trend, shows genuine GOE–GUE repulsion (⟨r⟩ = 0.632).

2. **The Galilean moons are a trap, and we don't fall for it.** Their raw ⟨r⟩ = 0.911
   looks "super-GUE," but after unfolding CV = 0.09 — they are **resonance-locked**
   (Laplace 1:2:4), pure periodicity, not RMT repulsion. The pipeline flags and excludes
   them automatically.

3. **Pooling 320 systems does NOT collapse to Poisson — it stays GUE (0.674).**
   This contradicts the naive superposition prediction. The reason: each planetary system
   is *internally* repulsive (dynamical packing / orbital clearing), and we unfold each
   system in its own frame, preserving that repulsion. Pooling aligned, internally
   repulsive systems reveals the **universality** of orbital repulsion across 320
   independent stars.

4. **The collapse is a local-spacetime property.** Destroy the per-system frames — throw
   all 1130 planets into one catalog ("global mix") — and ⟨r⟩ drops to 0.437, toward the
   Poisson null (0.396). Repulsion lives **within a single gravitational well**; mixing
   independent wells (different stellar masses, formation histories, dynamical clocks)
   erases it.

## Two statistical traps handled

- **Titius–Bode trap**: planetary ln(a) is near-geometric, so naive ⟨r⟩ is inflated by
  periodicity. Removed by **local unfolding** of each system.
- **Resonance trap**: mean-motion-resonant chains lock into exactly geometric spacing.
  Flagged by post-unfolding **CV < 0.15** and excluded from RMT interpretation.

## Reproduce

```bash
pip install -r requirements.txt
# Target A only (textbook constants, hardcoded):
python code/orbital_rmt_pipeline.py
# Targets A + B (supply the Kepler CSV):
python code/orbital_rmt_pipeline.py data/kepler_multiplanet.csv
```

## Data

`data/kepler_multiplanet.csv` — 1130 planets in 341 systems (≥3 planets each) from the
NASA Exoplanet Archive `pscomppars` table. TAP query:

```
https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+hostname,pl_name,pl_orbsmax,sy_pnum+from+pscomppars+where+pl_orbsmax>0+and+sy_pnum>2&format=csv
```

All 320 analyzed systems are Galactic (no resolved extragalactic orbital catalogs exist).

## Repository Structure

```
orbital-rmt-gravitational-shadow/
├── README.md
├── LICENSE                       # MIT
├── requirements.txt
├── .zenodo.json
├── CITATION.cff
├── paper/
│   ├── paper.tex
│   ├── paper.pdf
│   └── figs/
│       ├── fig1_orbital_repulsion.pdf
│       └── fig2_collapse_test.pdf
├── code/
│   └── orbital_rmt_pipeline.py   # unfolding, resonance flag, Targets A+B, collapse test
├── data/
│   └── kepler_multiplanet.csv    # 1130 planets, 341 systems (NASA Exoplanet Archive)
├── figures/                      # showcase PNGs + vector PDFs
└── results/
    └── orbital_rmt_results.json  # all locked statistics
```

## The unified program (8 systems)

| # | Domain | ⟨r⟩ | DOI |
|---|---|---|---|
| 1 | Stratigraphy | GOE | [20774581](https://zenodo.org/records/20774581) |
| 2 | Seismotectonics | scale-dependent | [20768130](https://zenodo.org/records/20768130) |
| 3 | Mantle plumes | 0.630 | [20768420](https://zenodo.org/records/20768420) |
| 4 | Metallogeny | 0.57–0.71 | [20768849](https://zenodo.org/records/20768849) |
| 5 | Evolution | 0.62 | [20783763](https://zenodo.org/records/20783763) |
| 6 | Hydrogeology | 0.55–0.74 | [20780389](https://zenodo.org/records/20780389) |
| 7 | Solar flares | Poisson (depletion limit) | [20784967](https://zenodo.org/records/20784967) |
| **8** | **Orbital architecture** | **0.674 (local-frame GUE)** | this work |

## License

MIT (code) · Data courtesy NASA Exoplanet Archive.
