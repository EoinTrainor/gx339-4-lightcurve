# Lightcurve of the Donor Star to GX 339-4

## HAWK-I Near-Infrared Image Reduction | ZOGY Difference Imaging | Ellipsoidal Lightcurve Modelling with ICARUS

---

### Project Overview

This repository contains the full data reduction and analysis pipeline for a research masters project aimed at constraining the black hole mass in the low-mass X-ray binary (LMXB) **GX 339-4**.

The approach uses ground-based near-infrared (NIR) imaging of the donor star obtained with **ESO HAWK-I** (High Acuity Wide field K-band Imager) at the VLT. By extracting a precision lightcurve of the donor star in the Ks band and modelling its ellipsoidal modulation, we derive the binary mass function and place constraints on the black hole mass.

---

### Science Background

GX 339-4 is one of the most active and well-studied black hole LMXBs, yet its black hole mass remains poorly constrained due to the system's persistent activity and the difficulty of isolating the donor star's flux. Two neighbouring stars in close angular proximity to GX 339-4 complicate aperture photometry, motivating the use of **ZOGY difference imaging** to isolate the variability of the donor star from contaminating sources.

---

### Data

- **Instrument:** ESO HAWK-I (VLT UT4, Yepun), Ks band
- **Detector of interest:** Bottom-left chip (Detector 1) — GX 339-4 imaged in all OBs
- **Dataset:** 12 Observing Blocks (OBs), 317 raw science frames
- **Frame format:** DIT=10s, NDIT=9 (9 sub-exposures averaged per FITS file; 90s effective sky time per frame)
- **Calibration frames:** 240 dark frames (EXPTIME=10s, NDIT=9), 234 Ks sky flats (DIT=1.68s and 3.5s)
- **Flat selection:** DIT=1.68s flats used for master flat construction — higher signal-to-noise than the DIT=3.5s group and sufficient for accurate flat-fielding

### Observing Log

| OB | Date (UT) | Frames | Total Exp (s) | Notes |
|----|-----------|--------|---------------|-------|
| GX339_Ks_Imaging_1  | 2025-05-17 | 26 | 2340 | |
| GX339_Ks_Imaging_2  | 2025-06-02 | 26 | 2340 | |
| GX339_Ks_Imaging_3  | 2025-06-04 | 26 | 2340 | |
| GX339_Ks_Imaging_4  | 2025-06-19 | 26 | 2340 | |
| GX339_Ks_Imaging_5  | 2025-07-03 | 26 | 2340 | |
| GX339_Ks_Imaging_6  | 2025-07-14 | 26 | 2340 | |
| GX339_Ks_Imaging_7  | 2025-07-18 | 26 | 2340 | |
| GX339_Ks_Imaging_8  | 2025-07-20 | 26 | 2340 | |
| GX339_Ks_Imaging_9  | 2025-07-23 | 31 | 810 + 1980 | **OB interrupted mid-execution** — ESO ADP split into two products |
| GX339_Ks_Imaging_10 | 2025-08-15 | 26 | 2340 | |
| GX339_Ks_Imaging_11 | 2025-09-20 | 26 | 2340 | Outburst |
| GX339_Ks_Imaging_12 | 2025-08-30 | 26 | 2340 | Outburst |
| **Total** | | **317** | **~28,530** | |

### ESO Archive Data Products (ADP)

Alongside the raw frames, ESO delivered **39 pipeline-reduced Archive Data Products** (ADPs) — one set of 3 files per OB:

| File type | Description |
|-----------|-------------|
| `TILED_IMAGE` | Dark-subtracted, flat-fielded, sky-subtracted, co-added mosaic of all 4 HAWK-I detectors (~4800×4800 px) |
| `TILED_VAR_MAP` | Per-pixel variance map (noise estimate) |
| `TILED_CONFIDENCE_MAP` | Per-pixel weight map (number of contributing frames) |

These products are stored in `GX 339-4 Output/eso_adp_fits/` and are useful as an independent quality reference and astrometric check. They **cannot be used for the lightcurve** directly — each ADP stacks an entire OB into a single image, discarding all time resolution. Per-frame photometry requires the custom reduction pipeline in this repository.

### Data Selection

Of the 12 OBs, the final 2 were obtained during a GX 339-4 outburst and are **excluded from lightcurve modelling**. Outburst state classification was performed using **Swift/BAT** hard X-ray monitoring data, with OB epochs cross-matched against the long-term BAT light curve.

| OBs | Status | Used in lightcurve modelling? | Used as ZOGY reference? |
|-----|--------|-------------------------------|-------------------------|
| 1 – 9 | Quiescent | ✅ Yes | ✅ Yes |
| 10 | Suspected non-quiescent | ✅ Yes (science frames) | ❌ No — excluded from reference coadd |
| 11 – 12 | Outburst (Swift/BAT confirmed) | ❌ No — excluded because of accretion disc activity | ❌ No |

During outburst, accretion disc and jet emission contribute significantly to the NIR flux, diluting the donor star's ellipsoidal modulation signal and rendering the mass function derivation unreliable. The X-ray classification and OB epoch cross-match are documented in `06_xray_analysis/outburst/`.

> **Note on campaign state:** Phase-folded lightcurve analysis (Stage 10) reveals a monotonic brightening of GX 339-4 over the 66-day observing window (factor ~4 rise in ZOGY differential flux), indicating the system was rising toward the confirmed outburst in OBs 11–12. OB7 (2025-07-18) shows anomalously elevated flux and may represent early outburst activity. The ellipsoidal modulation signal (~5–15% amplitude) cannot be cleanly separated from this long-term trend with 9 orbital phase measurements. Future observations should target confirmed quiescent epochs over 3–5 consecutive nights to achieve full phase coverage before variability builds up.

> **Note:** Raw data, calibration frames, and pipeline outputs are stored locally and are **not tracked by this repository**. See `config.py` for path configuration.

---

### Repository Structure

```
├── README.md               # This file
├── .gitignore              # Excludes data, outputs, and Python clutter
├── requirements.txt        # Python dependencies
├── config.py               # All data paths and parameters — edit for your local setup
│
├── 01_data_prep/           # One-time data sorting and subset selection scripts
│   ├── 01_organise_science.py
│   ├── 02_flatten_fits.py
│   ├── 03_build_header_inventory.py
│   ├── 04_sort_by_type_and_filter.py
│   ├── 05_calibration_summary.py
│   ├── 06_build_ks_working_subset.py
│   └── 07_inspect_ks_working_subset.py
│
├── 02_image_calibration/   # Image calibration pipeline
│   ├── 00_diagnostics.py           # Inspect raw FITS structure and headers
│   ├── 01_inspect_dark_frames.py   # Per-frame dark statistics and outlier check
│   ├── 02_build_master_dark.py     # Median-combine darks, build master dark
│   ├── 03_build_bad_pixel_mask.py  # Build bad pixel mask from dark + RMS map
│   ├── 04_split_flats_by_dit.py    # Split flat frames by DIT into groups
│   ├── 05_build_master_flat.py     # Build normalised master flat field
│   ├── 06_calibrate.py             # Dark subtract, flat field, sky subtract all 317 frames
│   └── 07_calibration_summary.py   # Global sky statistics plots across all OBs
│
├── 03_alignment/           # Astrometric alignment to common WCS
│
├── 04_zogy_difference_imaging/  # ZOGY difference imaging (PSF modelling + subtraction)
│   ├── 09_zogy.py                  # Build coadded reference, model PSFs, run ZOGY Eq. 13–17; outputs lightcurve_raw.csv
│   ├── 09b_visualise_zogy.py       # Report-quality figures: reference image, before/after, D histogram, quality timeline
│   └── 10_lightcurve.py            # Phase-fold S_target on Heida+2017 ephemeris; per-OB detrending; lightcurve plots
│
├── 05_photometry/          # Aperture/PSF flux extraction and lightcurve
│   ├── aperture_v1/            # First-pass circular aperture photometry (batch, imaging)
│   │   └── draft1/             # Early iterative draft scripts (5.1–5.8)
│   └── aperture_v2/            # WCS-tracked aperture pipeline with comparison stars,
│                               #   2D Gaussian centroid fitting, SNR optimisation, lightcurve
│
├── 06_xray_analysis/       # Swift/BAT X-ray monitoring and orbital phase analysis
│   ├── orbital_phase/          # MJD extraction, orbital phase (Heida+2017), coverage charts
│   └── outburst/               # Swift/BAT lightcurve overlay, OB state qualification
│
├── 07_modelling/           # ICARUS ellipsoidal lightcurve modelling
│
└── notebooks/              # Exploratory scripts for FITS inspection and visualisation
    ├── mef_science_images/     # HAWK-I MEF viewer, chronologiser, text reporter (4-chip)
    └── science_images/         # Single-detector image viewer, locator/timebar scripts
```

---

### Pipeline Stages

| Stage | Folder | Script | Status | Description |
|-------|--------|--------|--------|-------------|
| 1 | `02_image_calibration/` | `02_build_master_dark.py` | ✅ Done | Median-combine 240 dark frames, save master dark + RMS map |
| 2 | `02_image_calibration/` | `03_build_bad_pixel_mask.py` | ✅ Done | Build bad pixel mask from dark level + RMS (8σ threshold) |
| 3 | `02_image_calibration/` | `04_split_flats_by_dit.py` | ✅ Done | Split 234 flat frames into DIT=1.68s and DIT=3.5s groups |
| 4 | `02_image_calibration/` | `05_build_master_flat.py` | ✅ Done | Median-combine DIT=1.68s flats, normalise to median=1 |
| 5 | `02_image_calibration/` | `06_calibrate.py` | ✅ Done | Dark subtract, flat field, sky subtract all 317 science frames |
| 6 | `02_image_calibration/` | `07_calibration_summary.py` | ✅ Done | Global sky statistics and quality plots across all 12 OBs |
| 7 | `02_image_calibration/` | `08_reduction_summary.py` | ✅ Done | Global reduction CSV — all 317 frames, raw→aligned provenance |
| 8 | `03_alignment/` | `08_align.py` | ✅ Done | Align 314/317 frames to common reference (astroalign, star-triangle) |
| 9 | `04_zogy_difference_imaging/` | `09_zogy.py` | ✅ Done | ZOGY proper image subtraction (Zackay+2016) Eq. 13–17; S statistic + 8 reference stars extracted per frame; 236 frames, D_std(3σ-clip)=0.73, 5σ depth Ks≈17.1 |
| 10 | `04_zogy_difference_imaging/` | `10_lightcurve.py` | ✅ Done | Phase-fold S_target on Heida+2017 ephemeris (P=1.7587 d); per-OB detrending; 9 phase points recovered. Long-term brightening trend over 66 days dominates ellipsoidal signal — GX 339-4 was entering an active state during the campaign. New observations in confirmed quiescence required. |
| 11 | `06_xray_analysis/` | — | ✅ Done | Cross-match OB epochs with Swift/BAT; confirm outburst OBs |
| 12 | `07_modelling/` | `11_icarus_model.py` | 🔲 Pending | Ellipsoidal lightcurve modelling with ICARUS |

---

### Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/EoinTrainor/gx339-4-hawki-reduction.git
cd gx339-4-hawki-reduction
pip install -r requirements.txt
```

Then edit `config.py` to point to your local data directories before running any pipeline script.

---

### Dependencies

See `requirements.txt` for the full list. Core dependencies:

- `astropy` — FITS I/O, WCS, units
- `photutils` — aperture and PSF photometry
- `numpy` — array operations
- `matplotlib` — visualisation
- `scipy` — interpolation and convolution
- `astroalign` — image alignment
- `sep` — fast source extraction
- `astroquery` — 2MASS catalog queries for photometric calibration

---

### Usage

Edit `config.py` with your local paths, then run scripts in pipeline order:

```bash
python 02_image_calibration/00_diagnostics.py       # inspect raw data first
python 02_image_calibration/02_build_master_dark.py
python 02_image_calibration/03_build_bad_pixel_mask.py
python 02_image_calibration/04_split_flats_by_dit.py
python 02_image_calibration/05_build_master_flat.py
python 02_image_calibration/06_calibrate.py
python 02_image_calibration/07_calibration_summary.py
python 02_image_calibration/08_reduction_summary.py
python 03_alignment/08_align.py
python 04_zogy_difference_imaging/09_zogy.py
# photometry — in development
```

Each script saves outputs to the directories defined in `config.py` under `OUTPUT_ROOT`.

---

### Authors

- Eoin Trainor — Research Masters Student, University College Cork
- Supervisor: Dr. Mark Kennedy

---

### References

- Zackay, Ofek & Gal-Yam (2016) — ZOGY difference imaging algorithm
- Breton et al. (2012) — ICARUS lightcurve modelling code
- ESO HAWK-I instrument documentation

---

### License

This code is developed for academic research purposes. Please contact the author before reuse or reproduction.
