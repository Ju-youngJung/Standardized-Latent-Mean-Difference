# Estimator Dependence of Standardized Latent Mean Differences

Code and output for

> Jung, J. (2026). *Evaluating estimator dependence of standardized latent mean
> differences under approximate measurement invariance with binary indicators.*
> Manuscript submitted to *Research Synthesis Methods*.

The repository contains one R script that reproduces every simulation result,
the empirical analysis of the Rosenberg Self-Esteem Scale, and all tables and
figures reported in the manuscript.

---

## Contents

```
.
├── latent_effect_size_gmm.R   # complete analysis (simulation, empirical, tables, figures)
├── README.md
├── LICENSE
└── results/                   # outputs reported in the manuscript (optional, see below)
    ├── output_simulation/summary.csv
    ├── output_empirical/
    ├── output_tables/
    └── output_figures/
```

Running the script creates four folders in the working directory:

| Folder | Contents |
|---|---|
| `output_simulation/` | One checkpoint per design cell (`cell_###.rds`) and the cell summary (`summary.csv`, `summary.rds`) |
| `output_empirical/` | Downloaded and cleaned RSE data, full-sample estimates (`effect_sizes.csv`), item tests (`item_tests.csv`), moment residuals (`moment_residuals.csv`), descriptives, and subsample estimates (`subsample_estimates.csv`) |
| `output_tables/` | Tables 1 to 9 as CSV and Markdown |
| `output_figures/` | Figures as PNG (300 dpi) |

---

## Requirements

- R (version 4.3 or later recommended)
- Required package: `pbivnorm`
- Comparator and figure packages: `lavaan` (WLSMV), `mirt` (MML), `ggplot2` (figures)
- `parallel` (included with R) for multicore execution on Linux and macOS

Missing packages are installed automatically on first run. The versions used for
the manuscript are listed under [Software versions](#software-versions).

---

## Quick start

### Full reproduction

```r
setwd("<folder containing latent_effect_size_gmm.R>")
source("latent_effect_size_gmm.R")
```

or from a terminal:

```bash
Rscript latent_effect_size_gmm.R > run_log.txt 2>&1
```

The full run (64 cells × 1,000 replications plus the empirical study) took about
17 hours on 9 cores. With `cfg$clean_start = TRUE`, the four output folders are
deleted first. An interrupted run resumes from the last completed cell if
`cfg$clean_start` is set to `FALSE` before the script is re-sourced.

### Quick check (about 1 hour)

To verify that the code runs and that results are in line with the manuscript
without the full run time, reduce the number of replications in Section 0 of the
script before sourcing it:

```r
cfg$n_rep <- 50
```

With 50 replications, Monte Carlo error is roughly 4.5 times larger than in the
manuscript (about .03 for coverage near .95), so values will differ from the
reported ones in the second decimal, but the qualitative pattern should match.

### Tables and figures only

If `output_simulation/summary.csv` and the files in `output_empirical/` are
already present (for example, copied from `results/`), the tables and figures can
be regenerated without rerunning the analyses:

```r
cfg$clean_start    <- FALSE
cfg$run_simulation <- FALSE
cfg$run_empirical  <- FALSE
cfg$make_outputs   <- TRUE
```

---

## Configuration

All settings are in the `cfg` list in Section 0 of the script.

| Setting | Default | Meaning |
|---|---|---|
| `n_rep` | 1000 | Replications per design cell |
| `seed` | 20261030 | Base seed; replication seeds are `seed + 100000 × cell + replication` |
| `n_grid` | 500, 1000, 2000 | Sample size per group |
| `J_grid` | 5, 10, 20 | Number of binary items |
| `d_grid` | 0, 0.4 | True standardized latent mean difference |
| `dif_grid` | none, uniform, nonuniform, skewlatent | Measurement conditions |
| `dif_prop` | 0.20 | Proportion of noninvariant items (rounded up) |
| `dif_shift` | 0.5 | Uniform DIF: latent response shift in the reference group |
| `dif_load` | 0.7 | Nonuniform DIF: loading multiplier in the reference group |
| `latent_skew` | 0.5 | Skewness of the latent trait in the skewed condition |
| `drop_small_n_long_test` | TRUE | Exclude the eight cells with 500 per group and 20 items |
| `emp_neg_items` | 3, 5, 8, 9, 10 | Negatively worded RSE items (reverse scored) |
| `emp_subsamples`, `emp_sub_n` | 200, 500 | Real-data Monte Carlo: number and size per group of subsamples |
| `cores` | detected cores − 1 | Parallel workers |

---

## Script structure

| Section | Content |
|---|---|
| 0 | Configuration |
| 1 | Population and data generation |
| 2 | Model-implied probabilities and analytic gradients |
| 3 | Model specification with equality constraints (full and partial invariance) |
| 4 | Estimators: iterated GMM, delta-method effect size, misspecification-robust variance, overidentification test, item-level subset moment (C) and restriction (D) tests, sequential exclusion and partial invariance procedures, sample splitting, and the WLSMV and MML comparators |
| 5 | Simulation driver with checkpointing and stage-wise failure recording |
| 6 | Empirical study: download, cleaning, estimation, moment residuals, real-data Monte Carlo |
| 7 | Tables |
| 8 | Figures |
| 9 | Main: consistency check, simulation, empirical study, outputs |

The consistency check at the start of a run fits WLSMV, MML, and GMM to one
correctly specified data set and stops if they disagree beyond sampling error,
which guards against specification errors in the comparators.

---

## Data

The empirical example uses the Rosenberg Self-Esteem Scale data from the Open
Psychometrics archive (<https://openpsychometrics.org/_rawdata/>). The script
downloads `RSE.zip` automatically from
<https://openpsychometrics.org/_rawdata/RSE.zip>. The data are anonymous,
publicly available responses to an online questionnaire; no new data were
collected. If the download fails, place `RSE.zip` in `output_empirical/`
manually and rerun.

Cleaning rules: respondents who reported their gender as female or male,
answered all ten items, and reported an age between 13 and 90 are retained
(N = 44,918). Negatively worded items are reverse scored, and each item is
dichotomized at the scale midpoint.

---

## Mapping to the manuscript

| Manuscript | Source |
|---|---|
| Tables 1 to 3 (performance under invariance) | `output_simulation/summary.csv`, `dif == "none"`, `d == 0.4` |
| Table 4 (uniform DIF) | `summary.csv`, `dif == "uniform"` |
| Table 5 (detection and localization) | `summary.csv`, `est_id == "excl"` (columns `J_reject`, `exact_*`, `bias`, `cover`) |
| Tables 6 and 7 (empirical) | `output_empirical/effect_sizes.csv`, `item_tests.csv` |
| Table B1 (pseudo-true values) | `summary.csv`, column `target` |
| Tables C1 and C2 (no true difference) | `summary.csv`, `dif == "none"`, `d == 0` |
| Figures | `output_figures/` |

In `summary.csv`, `est_id` identifies the estimator: `obs` (observed d),
`wlsmv`, `mml`, `conv` (GMM, conventional SE), `rob` (GMM, robust SE), `excl`
and `part` (exclusion and partial invariance), and `sx` and `sp` (their sample
splitting versions).

---

## Software versions

Fill in from `sessionInfo()` of the run used for the manuscript:

| Software | Version |
|---|---|
| R | x.y.z |
| pbivnorm | x.y |
| lavaan | x.y.z |
| mirt | x.y.z |
| ggplot2 | x.y.z |

---

## License

The code is released under the MIT License (see `LICENSE`).

## Citation

If you use this code, please cite the article above. Questions can be sent to
Juyoung Jung, The University of Iowa (juyoung-jung@uiowa.edu).
