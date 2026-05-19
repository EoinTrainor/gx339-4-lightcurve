"""
config.py
---------
Central configuration for the GX 339-4 HAWK-I reduction pipeline.

Edit the paths in this file to match your local data directory structure.
This file IS tracked by git — do not store passwords or credentials here.
"""

from pathlib import Path

# ─── Root data directory ──────────────────────────────────────────────────────
# Change this to wherever your three data folders live on your machine
DATA_ROOT = Path(r"C:\Astronomy\GX 339-4 Raw Data\ESO_RAW_GX339_4")

# ─── Input data directories ───────────────────────────────────────────────────
DARKS_DIR   = DATA_ROOT / "Darks_match_10.0_NDIT9"   # Raw dark frames (Ks band, EXPTIME=10s, NDIT=9)
FLATS_DIR   = DATA_ROOT / "Flats_Ks"                  # Raw flat frames (Ks band)
SCIENCE_DIR = DATA_ROOT / "Science_Ks"                # Raw science frames (Ks band)

# ─── Flat subgroups (split by DIT for correct dark scaling) ──────────────────
FLAT_GROUPS_DIR = Path(r"C:\Astronomy\GX 339-4 Raw Data\ESO_RAW_GX339_4\Ks_Flat_Groups")
FLATS_DIT_1p68  = FLAT_GROUPS_DIR / "DIT_1p676206"   # Used for master flat
FLATS_DIT_3p5   = FLAT_GROUPS_DIR / "DIT_3p5"

# ─── Output directories ───────────────────────────────────────────────────────
OUTPUT_ROOT      = Path(r"C:\Astronomy\GX 339-4\GX 339-4 Output")
MASTERS_DIR      = OUTPUT_ROOT / "masters"        # Master dark, flat, pixel mask
CALIBRATED_DIR   = OUTPUT_ROOT / "calibrated"     # Calibrated science frames (per-OB subfolders)
ALIGNED_DIR      = OUTPUT_ROOT / "aligned"        # Aligned calibrated frames
DIFF_DIR         = OUTPUT_ROOT / "difference"     # ZOGY difference images
LIGHTCURVE_DIR   = OUTPUT_ROOT / "lightcurves"    # Final lightcurve files
LOGS_DIR         = OUTPUT_ROOT / "logs"           # Pipeline logs (root)

# ─── Log subfolders ───────────────────────────────────────────────────────────
LOGS_MASTERS_DIR     = LOGS_DIR / "masters"             # Dark, flat, BPM reports
LOGS_DATA_PREP_DIR   = LOGS_DIR / "data_prep"           # Data preparation logs
LOGS_CALIBRATION_DIR = LOGS_DIR / "calibration"         # Per-OB calibration reports
LOGS_CAL_SUMMARY_DIR = LOGS_CALIBRATION_DIR / "summary_plots"  # Global summary plots
LOGS_ALIGNMENT_DIR   = LOGS_DIR / "alignment"           # Per-OB alignment reports

# ─── Calibration product filenames (in MASTERS_DIR) ─────────────────────────
MASTER_DARK_FILE   = MASTERS_DIR / "master_dark_Ks_science_match.fits"
MASTER_DARK_RMS    = MASTERS_DIR / "master_dark_Ks_science_match_rms.fits"
BAD_PIXEL_MASK     = MASTERS_DIR / "bad_pixel_mask_from_dark.fits"
MASTER_FLAT_FILE   = MASTERS_DIR / "master_flat_Ks_DIT1p676206.fits"
MASTER_FLAT_RMS    = MASTERS_DIR / "master_flat_Ks_DIT1p676206_rms.fits"

# ─── HAWK-I instrument configuration ─────────────────────────────────────────
DETECTOR_OF_INTEREST = 1          # Bottom-left chip (1-indexed); GX 339-4 location
N_DETECTORS          = 4          # HAWK-I has 4 Hawaii-2RG detectors
FILTER               = "Ks"       # Observing filter

# ─── Target configuration ─────────────────────────────────────────────────────
TARGET_NAME = "GX339-4"
TARGET_RA   = 255.705780          # degrees (J2000) — NASA/SIMBAD
TARGET_DEC  = -48.789750          # degrees (J2000) — NASA/SIMBAD
# Approximate pixel seed for centroid fitting in the aligned reference frame.
# WCS carries ~31px residual in Y for this crowded field; this visual position
# is used as the centroid starting point with a tight (8px) search radius.
TARGET_PIX_SEED = (750, 1012)     # (x, y), 0-indexed — visually confirmed

# ─── Calibration parameters ───────────────────────────────────────────────────
DARK_SIGMA_CLIP      = 3.0        # Sigma threshold for hot pixel masking in dark
FLAT_SIGMA_CLIP      = 3.0        # Sigma threshold for flat normalisation
SKY_SIGMA_CLIP       = 3.0        # Sigma for sky frame sigma-clipping
MIN_FLAT_VALUE       = 0.5        # Fraction of median below which pixels are masked
MAX_FLAT_VALUE       = 1.5        # Fraction of median above which pixels are masked

# ─── ZOGY parameters (placeholders — tuned during development) ────────────────
ZOGY_SIG_PSF         = 2.0        # Initial PSF sigma estimate (pixels)
ZOGY_SIG_CLIPPING    = 3.0        # Sigma clipping for ZOGY background estimation

# ─── Photometry parameters ───────────────────────────────────────────────────
APERTURE_RADIUS      = 5.0        # Circular aperture radius (pixels) — TBD from FWHM
SKY_ANNULUS_INNER    = 8.0        # Inner sky annulus radius (pixels)
SKY_ANNULUS_OUTER    = 12.0       # Outer sky annulus radius (pixels)


# ─── Utility: create output directories if they don't exist ───────────────────
def make_output_dirs():
    for d in [MASTERS_DIR, CALIBRATED_DIR, ALIGNED_DIR,
              DIFF_DIR, LIGHTCURVE_DIR,
              LOGS_MASTERS_DIR, LOGS_DATA_PREP_DIR,
              LOGS_CALIBRATION_DIR, LOGS_CAL_SUMMARY_DIR,
              LOGS_ALIGNMENT_DIR]:
        d.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    make_output_dirs()
    print("Output directories created successfully.")
    print(f"  Data root : {DATA_ROOT}")
    print(f"  Output root: {OUTPUT_ROOT}")
