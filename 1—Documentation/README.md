
# Vegetation Index Mapping — Technical Documentation

Oct 23, 2026 · @Syed Jafri

## Overview

The project computes NDVI, NDRE and GNDVI maps from repeat multispectral drone flights and tests, stage by stage through the season, which index best separates stressed from healthy crop.

**Research question** — at each growth stage, which index gives the largest, most reliable difference between stressed and unstressed plots, and how early does each index detect stress?

**Working hypotheses** (to test, not assume)

- NDVI performs well early, when canopy is sparse, but saturates once canopy closes.
- NDRE stays sensitive in dense canopy because red-edge light penetrates deeper and tracks chlorophyll.
- GNDVI sits between the two, with less saturation than NDVI.

**Scope**

- In: multispectral orthomosaics (Green, Red, Red Edge, NIR), plot or zone boundaries, stress ground truth (treatments or scouting), growth stage records.
- Out: index rasters per flight, plot-level index tables, per-stage comparison statistics, ranked index recommendation per stage, and maps/figures.
- Not in scope: diagnosing the cause of stress, or yield prediction (possible extension).

## Pipeline

Imagery and ground truth meet in one plot × date × index table, so every index is judged on identical plots, dates and masks.

&#91;embedded content: pipeline · flights to stage-wise index ranking\]

Stack: Python 3.11, rasterio/GDAL, numpy, geopandas with exactextract for plot statistics, statsmodels for mixed models, scikit-learn for ROC metrics, and a Typer CLI `vimap`.

## Imagery acquisition and preprocessing

Indices are only comparable across dates if every flight is converted to calibrated surface reflectance with the same sensor and similar sun conditions.

| Parameter | Requirement |
| --- | --- |
| Sensor | One multispectral camera all season with Green, Red, Red Edge (\~717 or 730 nm) and NIR bands; record band centres and widths |
| Ground sample distance | 3–8 cm; consistent across flights (±20%) |
| Timing | Within ±2 h of solar noon; avoid broken cloud |
| Calibration | Reflectance panel imaged before and after each flight; downwelling light sensor on |
| Overlap | 75–80% front, 70% side |
| Revisit | At least one flight per key growth stage; weekly preferred |
| Positioning | RTK/PPK or ≥ 5 GCPs for plot-level extraction |

**Preprocessing**

1. **Orthomosaic + reflectance** — process in Pix4D, Metashape or OpenDroneMap with panel and DLS calibration to reflectance (0–1).
2. **Calibration QA** — check reflectance of stable targets (panel, bare road) across flights; flag flights with > 10% drift.
3. **Co-registration** — align all flights to a reference grid so plot extraction uses identical pixels.
4. **Soil and shadow masking** — exclude pixels with NDVI < 0.2 (tunable) and low NIR brightness shadow pixels, so indices describe the canopy rather than the canopy-soil mix. Results are also reported unmasked for comparison, because early-season masking changes NDVI's behaviour.

## Index computation

All three indices are normalized differences of NIR against one visible or red-edge band, computed per pixel on reflectance with identical masking, so differences between them come from the band, not the processing.

| Index | Formula | Contrast band | Expected strength |
| --- | --- | --- | --- |
| NDVI | (NIR − Red) / (NIR + Red) | Red, strongly absorbed by chlorophyll | Canopy cover and biomass; saturates in dense canopy |
| NDRE | (NIR − RedEdge) / (NIR + RedEdge) | Red edge, partly absorbed | Chlorophyll / nitrogen in dense canopy |
| GNDVI | (NIR − Green) / (NIR + Green) | Green, weakly absorbed | Chlorophyll; intermediate saturation |

**Implementation**

```python
def nd(a, b):
    num, den = a - b, a + b
    return np.where(den > 0.01, num / den, np.nan)   # guard near-zero denominators

ndvi  = nd(nir, red)
ndre  = nd(nir, red_edge)
gndvi = nd(nir, green)
```

- Output: float32 Cloud-Optimized GeoTIFFs, nodata = NaN, one file per index per flight.
- Values outside \[−1, 1\] are impossible after the guard; values outside the plausible range for canopy (e.g. NDVI < 0) are flagged in QA.
- Because index ranges differ (NDRE values are typically lower than NDVI), comparisons use standardized effect sizes, never raw differences.

**Plot extraction** — for each plot and flight: mean, median, SD, 10th and 90th percentiles of each index over valid pixels inside the plot, shrunk by an edge buffer (default 0.3 m).

## Stress ground truth and growth stages

The comparison needs plots whose stress status is known independently of the imagery, and a growth stage attached to every flight.

**Stress ground truth (choose one or combine)**

| Source | Design | Strength |
| --- | --- | --- |
| Controlled treatments | Replicated plots with e.g. nitrogen rates (0, 50, 100, 150% of recommended) or irrigation levels | Known cause and timing; cleanest comparison |
| Field scouting | Georeferenced points scored for stress severity on flight dates | Real-world conditions |
| Physiological measures | SPAD/chlorophyll meter, leaf N, stomatal conductance or leaf water potential on sampled plants | Continuous reference, not just a class |

Treatments are laid out in randomized complete blocks with ≥ 4 replicates, so stress is not confounded with field position.

**Growth stages** — record the stage on each flight date with a standard scale (e.g. V-stages / R-stages for maize and soybean, or BBCH for cereals). For analysis, flights are grouped into stage windows:

| Window | Example (maize) | Canopy condition |
| --- | --- | --- |
| Early vegetative | VE–V5 | Sparse, much soil visible |
| Mid vegetative | V6–V12 | Rapid growth, canopy closing |
| Late vegetative / reproductive onset | V13–R1 | Closed canopy |
| Reproductive | R2–R5 | Dense canopy, then senescence begins |

Adapt windows to the crop; the key is that each window has a distinct canopy density, since that drives the expected index behaviour.

## Comparison method

For every growth-stage window, each index is scored on the same plots with the same metrics, and differences between indices are tested on paired data, because all three come from the same pixels.

**Metrics per index per stage**

| Metric | Measures | Computed as |
| --- | --- | --- |
| Cohen's d | Separation of stressed vs healthy | (mean healthy − mean stressed) / pooled SD |
| ROC-AUC | Classification ability, threshold-free | Plot index value as the score, stress class as truth |
| Correlation with reference | Tracking of continuous stress | Spearman ρ vs SPAD, leaf N or treatment rate |
| Coefficient of variation in healthy plots | Noise floor | SD / mean of healthy plots |
| Saturation check | Remaining dynamic range | Share of healthy pixels within 0.02 of the index's seasonal 99th percentile |
| First detection date | Earliness | First flight where stressed vs healthy differ at p < 0.05 (mixed model below) |

**Statistical tests**

- **Treatment effect per index and stage** — linear mixed model: index \~ treatment × stage + (1 | block) + (1 | plot), on standardized index values.
- **Comparing indices** — paired bootstrap (resampling plots, 2 000 iterations) of the difference in AUC and in Cohen's d between each pair of indices; DeLong's test for AUC as a check. An index "wins" a stage only if its 95% interval excludes the others.
- **Multiple comparisons** — Holm correction across stages × index pairs.

**Outputs** — a stage × index results table, a line plot of Cohen's d (with intervals) per index across stages, ROC curves per stage, and side-by-side index maps for one representative flight per stage.

## Configuration, CLI and outputs

One YAML file per season defines bands, plots, stage windows and analysis settings; adding a new index is one line in the config plus its band formula.

**`config/season.yaml`**

```yaml
season: 2026_maize_trial
crop: maize
flights: data/ortho/{date}.tif
bands: {green: 2, red: 3, red_edge: 4, nir: 5}
plots: data/plots.gpkg              # plot_id, block, treatment
truth: data/truth.csv               # plot_id, date, stress_class, spad
stages: data/stages.csv             # date, stage
stage_windows:
  early_veg: [VE, V5]
  mid_veg:   [V6, V12]
  late_veg:  [V13, R1]
  reprod:    [R2, R5]
indices:
  NDVI:  [nir, red]
  NDRE:  [nir, red_edge]
  GNDVI: [nir, green]
mask: {ndvi_min: 0.2, report_unmasked: true}
plot_edge_buffer_m: 0.3
stats: {bootstrap: 2000, alpha: 0.05, correction: holm}
```

**CLI**

```bash
vimap compute  --config config/season.yaml          # index rasters for all flights
vimap extract  --config config/season.yaml          # plot table (parquet)
vimap compare  --config config/season.yaml          # stage-wise metrics and tests
vimap report   --config config/season.yaml          # HTML/PDF with maps and figures
```

**Output layout**

```
derived/indices/{date}_{index}.tif
derived/plots.parquet               # plot × date × index statistics
results/stage_index_metrics.csv     # one row per stage × index
results/pairwise_tests.csv
results/report.html
```

## Validation, limitations and roadmap

The stage-by-stage ranking is only trusted if it holds across replicates, and ideally across a second field or season.

**Validation**

- Calibration check: index values of stable targets constant across flights within ±0.02.
- Repeatability: a same-day repeat flight on at least one date; plot index values should correlate at r > 0.95.
- Robustness: rerun the comparison with and without soil masking, and with plot means vs medians; report whether the winning index changes.
- External check: repeat on a second field or season; agreement of the per-stage ranking is the main evidence for a general recommendation.

**Known limitations**

- Results depend on the sensor's exact red-edge band position and width; NDRE from different cameras is not directly comparable.
- Index differences can reflect canopy structure (leaf angle, row cover) as well as stress.
- Stress type matters: nitrogen, water and disease stress may favour different indices; conclusions are specific to the stress studied.
- Few flights per stage window limits statistical power; stage boundaries are partly arbitrary.

**Roadmap**

- [ ] Add indices: SAVI/OSAVI (soil-adjusted), CIre, MCARI
- [ ] Thermal band for water-stress comparison
- [ ] Combine indices in a multivariate model and compare with the best single index
- [ ] Recommend an index-by-stage schedule as a reusable field protocol
