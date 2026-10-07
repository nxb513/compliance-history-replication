# Replication package

**Whose Compliance History? Plants, Parent Firms, and the Information Used in Environmental Enforcement in the United States**

This package contains the analysis files, code and pre-specified analysis plan behind every table and figure in the article and its Supplementary material. All inputs are public data from United States federal agencies.

## Contents

| Folder | Content |
|---|---|
| `data/interim/` | `one_bad_plant_panel.parquet` (primary plant-year panel, 7,215 plants of 924 firms, 2010 to 2025); `formal_enf.parquet` (water program formal enforcement actions from ECHO) |
| `extension/` | Analysis files derived from public EPA data: plant groups, air and hazardous waste plant-year panels, Facility Registry Service links, plant attributes, ECHO Exporter fields for sample plants, penalty cases (`e21_units.parquet`) |
| `extension/scripts/` | Analysis scripts (see the map below); `ext_common.py` holds shared helpers |
| `extension/tables/` | Output tables (CSV) produced by the scripts |
| `figures/` | `make_figures.py` and Figs. 1, 2 and S1 (EPS, PDF and PNG) |
| `construction/` | Scripts that download the raw public data and build the primary panel |
| `docs/analysis_plan.md` | Analysis plan with dated addenda, each written before the corresponding analysis was run |

## Data sources (all public)

- United States Environmental Protection Agency, Enforcement and Compliance History Online (ECHO) data downloads, https://echo.epa.gov/tools/data-downloads: water program (ICIS-NPDES) downloads and quarterly noncompliance history; `ICIS-AIR_downloads.zip`; `rcra_downloads.zip`; `echo_exporter.zip`.
- Facility Registry Service national files (program links and NAICS codes).
- Toxics Release Inventory basic data files (parent company names and identification numbers).
- National Oceanic and Atmospheric Administration, Storm Events Database detail files 2009 to 2025 (Supplementary material S6 only).

Raw files are not redistributed because of their size. The download scripts in `construction/` record the files used; extracted raw tables are expected under `data/raw/` and `data/interim/` with the file names used in the scripts.

## How to run

Run every script from the root of this package (all paths are relative), for example `python extension/scripts/r4_linkage_strict.py`. Python 3.12 with pandas, numpy, statsmodels, scikit-learn, duckdb, pyarrow and matplotlib.

Scripts that run on the analysis files included here, without raw downloads: `r4_linkage_strict.py`, `cross_analysis.py`, `b1_targeting.py`, `b5_multiple_testing.py`, `e14_event.py`, `e15_distance_border.py`, `e16_targeting_allocation.py`, `e20_regimes.py`, `robustness_R1_R3.py` and `figures/make_figures.py`. Scripts that read raw ECHO or NOAA extracts (air violation history, stack tests and formal actions; hazardous waste enforcement; Facility Registry Service NAICS codes; Storm Events files) require the raw downloads.

## Map from article to code

| Article | Script | Output table(s) |
|---|---|---|
| Table S2, sample construction | `construction/build_primary_panel.py`, `extension/scripts/cross_prep_air.py`, `cross_prep_rcra.py`, `prep_groups_links.py` | panels in `extension/` |
| Table 2, variance components | `r4_linkage_strict.py`; `e18_e19.py` (stack tests); `reml_robustness_all.py` (optimizer check) | `R4_linkage_strict_reml.csv`, `E18_stacktest_partition.csv`, `REML_optimizer_robustness.csv` |
| Table 3, cross-program gap | `cross_analysis.py` | `E9_air_lpm.csv`, `E9_rcra_lpm.csv` |
| Section 4.2, covariance and worst plant | `e12_e13.py`, `cross_analysis.py` | `E13_crossmedia_decomposition.csv`, `E13_lead_lag.csv`, `Xb_worst_plant_coincidence.csv` |
| Fig. 1, Section 4.3 | `e15_distance_border.py`, `robustness_R1_R3.py` | `E15_pair_correlations_by_distance.csv`, `E15_contrasts_bootstrap.csv`, `R1_E8_bootstrap_ci.csv` |
| Fig. S1 | `e14_event.py` | `E14_event_coefficients.csv`, `E14_event_tests.csv` |
| Li and Lyon replication | `e17_spillover_liylon.py` | `E17_liylon_spillover.csv` |
| Fig. 2, Table 4 | `b1_targeting.py`, `e16_targeting_allocation.py` | `B1_auc_single.csv`, `B1_auc_combinations.csv`, `B1_auc_differences.csv`, `E16a_targeting_simulation.csv` |
| Table 5, Panel A | `e21_penalties.py` | `E21_penalty_history.csv`, `E21_units_descriptives.csv` |
| Table 5, Panel B | `e16_targeting_allocation.py` | `E16b_allocation_vs_risk.csv` |
| Supplementary material S3 | `robustness_R1_R3.py`, `b2_b4.py`, `b5_multiple_testing.py`, `e20_regimes.py` | `R2_*`, `R3_*`, `B2_*` to `B5_*`, `E20_*` |
| Supplementary material S6 | `e22_florida.py`, `e23_floods.py`, `e18_e19.py`, `e19_event.py` | `E22_*`, `E23_*`, `E19_*` |

## Notes

- Variance components use restricted maximum likelihood with the best of three optimizers (`reml_best` in `ext_common.py`).
- Firm bootstraps use fixed seeds; results reproduce exactly.
- In `docs/analysis_plan.md`, references to an earlier manuscript by the author were replaced by neutral wording when the plan was prepared for review; nothing else was changed.
