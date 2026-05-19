"""
09_zogy.py
----------
Stage 9: ZOGY proper image subtraction for GX 339-4 HAWK-I Ks-band data.

Phase 1  — Proper coaddition reference (Zackay & Ofek 2015b / ZOGY Eq. 22–24)
           All quiescent aligned frames (OBs 1–10).  Per-frame: PSF from 2D
           Moffat fit to isolated stars, σ from sigma-clipped background RMS,
           F from 2MASS Ks catalogue matching.  Accumulate in Fourier space.

Phase 2  — ZOGY per-frame subtraction (Zackay, Ofek & Gal-Yam 2016, Eq. 13–17)
           Each quiescent frame is differenced against the proper-coaddition
           reference.  Outputs: D (Eq. 13), P_D (Eq. 14), F_D (Eq. 15),
           and S (Eq. 17) — the matched-filter detection statistic.
           S is evaluated at GX 339-4's pixel position each frame;
           differential flux = S_target / F_D.

Phase 3  — Reference image quality metrics
           FWHM distribution, flux zero-point timeline, reference PSF profile,
           D pixel statistics, approximate depth estimate.

References
----------
Zackay, Ofek & Gal-Yam 2016, ApJ 830:27
Zackay & Ofek 2015b (proper coaddition)
"""

# ══════════════════════════════════════════════════════════════════════════════
# SETTINGS  — edit here before running
# ══════════════════════════════════════════════════════════════════════════════
REF_OBS        = [f"GX339_Ks_Imaging_{i}" for i in range(1, 10)]  # OBs 1–9 (exclude OB10–12: suspected non-quiescent)

PSF_STAMP_SIZE = 51       # pixels, must be odd
MAX_PSF_STARS  = 30       # max stars per frame used for PSF fit
DETECT_SIGMA   = 5.0      # sep source detection threshold (σ above background)
ISOLATION_FWHM = 6        # exclude stars with any neighbour within N × FWHM

TWOMASS_RADIUS = 7.0      # arcmin, 2MASS search radius around target
TWOMASS_KS_MIN = 9.5      # brightest allowed Ks (avoid HAWK-I saturation)
TWOMASS_KS_MAX = 14.5     # faintest allowed Ks (signal adequate for ZP fit)
APER_RADIUS_PX = 8        # circular aperture radius for instrumental flux (px)

PIXSCALE       = 0.106    # arcsec/pixel (HAWK-I)
EPSILON        = 1e-10    # FFT denominator regularisation

TEST_MODE      = False   # False: full run on all frames
TEST_N_FRAMES  = 5
# ══════════════════════════════════════════════════════════════════════════════

import sys, csv, json, warnings, time
from pathlib import Path
from astropy.time import Time

# Force UTF-8 output on Windows (avoids cp1252 errors with Greek/box chars)
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
from numpy.fft import fft2, ifft2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.modeling import models, fitting
import sep

from astroquery.vizier import Vizier

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config

# ── Output paths ──────────────────────────────────────────────────────────────
ALIGNED_DIR = config.ALIGNED_DIR
DIFF_DIR    = config.DIFF_DIR
LOGS_ZOGY   = config.LOGS_DIR / "zogy"
QUALITY_DIR = LOGS_ZOGY / "quality"
for d in [DIFF_DIR, LOGS_ZOGY, QUALITY_DIR]:
    d.mkdir(parents=True, exist_ok=True)

REF_FITS = DIFF_DIR / "reference_R.fits"
REF_PSF  = DIFF_DIR / "reference_Pr.fits"


# ══════════════════════════════════════════════════════════════════════════════
# HELPER: OB sort
# ══════════════════════════════════════════════════════════════════════════════
def ob_sort_key(name):
    try:
        return int(name.rsplit("_", 1)[-1])
    except ValueError:
        return name


def get_aligned_frames(ob_list, test_mode=True, test_n=5):
    """Return list of (ob_name, Path) for aligned frames in ob_list."""
    frames = []
    for ob_name in sorted(ob_list, key=ob_sort_key):
        ob_dir = ALIGNED_DIR / ob_name
        if not ob_dir.exists():
            print(f"  WARNING: {ob_dir} not found — skipping")
            continue
        fps = sorted(ob_dir.glob("HAWKI.*_1_cal_aligned.fits"))
        if test_mode:
            fps = fps[:test_n]
        frames.extend((ob_name, fp) for fp in fps)
    return frames


# ══════════════════════════════════════════════════════════════════════════════
# 1. 2MASS CATALOGUE  (queried once, cached)
# ══════════════════════════════════════════════════════════════════════════════
_cat_cache = None

def query_2mass(ra_deg, dec_deg, radius_arcmin=7.0):
    """Query 2MASS Ks via Vizier. Returns Astropy Table (cached after first call)."""
    global _cat_cache
    if _cat_cache is not None:
        return _cat_cache

    v = Vizier(columns=["RAJ2000", "DEJ2000", "Kmag", "e_Kmag"], row_limit=-1)
    coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
    result = v.query_region(coords, radius=radius_arcmin * u.arcmin,
                            catalog="II/246/out")
    if not result:
        print("  WARNING: 2MASS query returned no results — ZP will default to 1.0")
        return None

    cat = result[0]
    # Keep only rows with valid Kmag
    valid = np.isfinite(np.array(cat["Kmag"], dtype=float))
    cat   = cat[valid]
    _cat_cache = cat
    ks = np.array(cat["Kmag"], dtype=float)
    print(f"  2MASS catalogue: {len(cat)} stars  "
          f"Ks = {ks.min():.1f}–{ks.max():.1f} mag")
    return cat


# ══════════════════════════════════════════════════════════════════════════════
# 2. SOURCE EXTRACTION
# ══════════════════════════════════════════════════════════════════════════════
def extract_sources(data, detect_sigma=5.0):
    """
    sep source extraction on sky-subtracted image.
    Returns (objects recarray, background rms scalar).
    NaN pixels are masked before extraction.
    """
    d    = data.astype(np.float64)
    mask = ~np.isfinite(d)
    d[mask] = 0.0

    bkg = sep.Background(d, mask=mask, bw=64, bh=64, fw=3, fh=3)
    rms = float(np.nanmedian(bkg.rms()))

    d_sub = d - bkg.back()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        objects = sep.extract(d_sub, detect_sigma, err=bkg.rms(), mask=mask,
                              minarea=5, clean=True, clean_param=1.0,
                              deblend_nthresh=32, deblend_cont=0.005)
    return objects, rms


# ══════════════════════════════════════════════════════════════════════════════
# 3. PSF ESTIMATION — 2D MOFFAT FIT ON ISOLATED STARS
# ══════════════════════════════════════════════════════════════════════════════
def fit_psf_moffat(data, objects, stamp_size=51, max_stars=30,
                   isolation_fwhm=6):
    """
    Fit 2D Moffat to isolated bright stars.
    Returns (psf_stamp normalised to sum=1, fwhm_x_px, fwhm_y_px).
    Returns (None, nan, nan) if fewer than 3 usable stars found.
    Moffat f(r) = A*(1 + r²/γ²)^(-α) captures PSF wings better than Gaussian.
    """
    if len(objects) == 0:
        return None, np.nan, np.nan

    H, W  = data.shape
    half  = stamp_size // 2

    # Sort bright→faint
    order = np.argsort(objects["peak"])[::-1]
    objs  = objects[order]

    xs = objs["x"].astype(float)
    ys = objs["y"].astype(float)

    # Initial FWHM guess from second-moment semi-axes
    fwhm_est = float(np.nanmedian(2.355 * np.sqrt(0.5 * (objs["a"]**2 + objs["b"]**2))))
    fwhm_est = max(fwhm_est, 2.0)
    min_sep  = isolation_fwhm * fwhm_est

    # Isolation check: exclude any star with a neighbour within min_sep
    isolated = np.ones(len(objs), dtype=bool)
    for i in range(len(objs)):
        dx   = xs - xs[i]
        dy   = ys - ys[i]
        dist = np.sqrt(dx**2 + dy**2)
        dist[i] = np.inf  # ignore self
        if dist.min() < min_sep:
            isolated[i] = False

    objs = objs[isolated]

    # Exclude edge stars
    margin = half + 2
    ok     = ((objs["x"] > margin) & (objs["x"] < W - margin) &
              (objs["y"] > margin) & (objs["y"] < H - margin))
    objs   = objs[ok][:max_stars]

    if len(objs) < 3:
        return None, np.nan, np.nan

    fitter = fitting.LevMarLSQFitter()
    yy, xx = np.mgrid[0:stamp_size, 0:stamp_size]
    stamps, fwhms_x, fwhms_y = [], [], []

    # Initial gamma from fwhm_est assuming alpha=3.5 (typical seeing)
    alpha_init = 3.5
    gamma_init = fwhm_est / (2.0 * np.sqrt(2.0**(1.0 / alpha_init) - 1.0))

    for obj in objs:
        cx = int(round(obj["x"]))
        cy = int(round(obj["y"]))
        cutout = data[cy - half: cy + half + 1,
                      cx - half: cx + half + 1].copy()
        if cutout.shape != (stamp_size, stamp_size):
            continue
        if not np.all(np.isfinite(cutout)):
            continue
        cutout -= np.median(cutout)
        peak = cutout.max()
        if peak <= 0:
            continue

        m_init = models.Moffat2D(
            amplitude=peak,
            x_0=half, y_0=half,
            gamma=gamma_init,
            alpha=alpha_init,
            bounds={
                "x_0":   (half - 3, half + 3),
                "y_0":   (half - 3, half + 3),
                "gamma": (0.3, stamp_size / 2.0),
                "alpha": (1.0, 8.0),
            }
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = fitter(m_init, xx, yy, cutout)

        fwhm_px = m.fwhm  # 2*gamma*sqrt(2^(1/alpha) - 1)
        if not (1.0 < fwhm_px < stamp_size / 2):
            continue

        # Force PSF centre to (half, half) before evaluating the stamp.
        # The fit may land up to ±3 px off-centre; evaluating there embeds a
        # sub-pixel phase ramp in P̂ that shifts the PSF convolution and
        # prevents exact source cancellation in the ZOGY difference.
        m.x_0 = half
        m.y_0 = half
        stamp = m(xx, yy)
        s     = stamp.sum()
        if s <= 0:
            continue
        stamps.append(stamp / s)
        fwhms_x.append(fwhm_px)
        fwhms_y.append(fwhm_px)

    if len(stamps) < 3:
        return None, np.nan, np.nan

    psf = np.median(np.array(stamps), axis=0)
    psf /= psf.sum()
    return psf, float(np.median(fwhms_x)), float(np.median(fwhms_y))


# ══════════════════════════════════════════════════════════════════════════════
# 4. PSF → FFT  (origin-centred convention)
# ══════════════════════════════════════════════════════════════════════════════
def psf_to_hat(psf_stamp, shape):
    """
    Embed PSF stamp in image-sized array with its centre at pixel (0,0),
    then return FFT2.  Follows numpy FFT convention where DC is at [0,0].
    """
    h, w     = psf_stamp.shape
    psf_full = np.zeros(shape, dtype=np.float64)
    psf_full[:h, :w] = psf_stamp
    psf_full = np.roll(psf_full, -(h // 2), axis=0)
    psf_full = np.roll(psf_full, -(w // 2), axis=1)
    return fft2(psf_full)


# ══════════════════════════════════════════════════════════════════════════════
# 5. FLUX ZERO POINT FROM 2MASS Ks
# ══════════════════════════════════════════════════════════════════════════════
def estimate_zeropoint(data, objects, cat, ref_wcs,
                       aper_r=8, ks_min=9.5, ks_max=14.5):
    """
    Match sep-detected sources to 2MASS Ks catalogue using the reference
    frame WCS (common to all aligned frames).

    Returns F_j = median(f_inst / f_cat) where f_cat = 10^(-0.4 * Ks).
    Returns 1.0 if fewer than 3 clean matches are found.
    """
    if cat is None or ref_wcs is None or len(objects) < 5:
        return 1.0

    # Magnitude cut
    cat_ks  = np.array(cat["Kmag"],    dtype=float)
    cat_ra  = np.array(cat["RAJ2000"], dtype=float)
    cat_dec = np.array(cat["DEJ2000"], dtype=float)
    in_mag  = (cat_ks >= ks_min) & (cat_ks <= ks_max)
    cat_ks  = cat_ks[in_mag]
    cat_ra  = cat_ra[in_mag]
    cat_dec = cat_dec[in_mag]
    if len(cat_ks) < 3:
        return 1.0

    cat_sky = SkyCoord(ra=cat_ra * u.deg, dec=cat_dec * u.deg)

    # Project detected sources through the reference WCS
    src_x = objects["x"].astype(float)
    src_y = objects["y"].astype(float)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            radec = ref_wcs.all_pix2world(np.column_stack([src_x, src_y]), 0)
    except Exception:
        return 1.0
    src_sky = SkyCoord(ra=radec[:, 0] * u.deg, dec=radec[:, 1] * u.deg)

    # Match: catalogue → sources
    idx, sep2d, _ = cat_sky.match_to_catalog_sky(src_sky)
    close = sep2d < 2.0 * u.arcsec

    if close.sum() < 3:
        return 1.0

    matched_src = objects[idx[close]]
    matched_ks  = cat_ks[close]
    f_cat_arr   = 10.0 ** (-0.4 * matched_ks)

    H, W  = data.shape
    d     = data.astype(np.float64)
    d     = np.where(np.isfinite(d), d, 0.0)
    ratios = []

    # Sky annulus: inner = aper_r + 3 px, outer = aper_r + 12 px
    r_sky_in  = aper_r + 3
    r_sky_out = aper_r + 12

    for src, f_cat in zip(matched_src, f_cat_arr):
        cx = int(round(src["x"]))
        cy = int(round(src["y"]))
        r  = aper_r
        r_o = r_sky_out

        # Cutout large enough to contain sky annulus
        y0 = max(0, cy - r_o);  y1 = min(H, cy + r_o + 1)
        x0 = max(0, cx - r_o);  x1 = min(W, cx + r_o + 1)
        if y1 - y0 < 2 or x1 - x0 < 2:
            continue
        cut = d[y0:y1, x0:x1]
        yy  = np.arange(y0, y1) - cy
        xx  = np.arange(x0, x1) - cx
        r2  = xx[np.newaxis, :]**2 + yy[:, np.newaxis]**2

        aper_mask = r2 <= r**2
        sky_mask  = (r2 >= r_sky_in**2) & (r2 <= r_sky_out**2)

        if sky_mask.sum() < 10:
            continue

        # Sigma-clipped sky background per pixel
        sky_pix = cut[sky_mask]
        _, sky_med, _ = sigma_clipped_stats(sky_pix, sigma=3.0, maxiters=5)
        sky_level = sky_med

        n_aper = float(aper_mask.sum())
        f_inst = float(cut[aper_mask].sum()) - sky_level * n_aper
        if f_inst > 0 and f_cat > 0:
            ratios.append(f_inst / f_cat)

    if len(ratios) < 3:
        return 1.0

    # Sigma-clip ratios to remove outliers
    ratios = np.array(ratios)
    med    = np.median(ratios)
    std    = np.std(ratios)
    keep   = np.abs(ratios - med) < 3.0 * std
    if keep.sum() < 3:
        return float(med)
    return float(np.median(ratios[keep]))


# ══════════════════════════════════════════════════════════════════════════════
# 6. BACKGROUND NOISE
# ══════════════════════════════════════════════════════════════════════════════
def estimate_sigma(data):
    """Sigma-clipped RMS of finite pixels."""
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return 1.0
    _, _, std = sigma_clipped_stats(finite, sigma=3.0, maxiters=5)
    return float(std) if std > 0 else 1.0


# ══════════════════════════════════════════════════════════════════════════════
# 7. SOURCE CENTROID REFINEMENT
# ══════════════════════════════════════════════════════════════════════════════
def find_source_centroid(image, cx_init, cy_init, search_radius=20):
    """
    Refine a source pixel position using SEP detection in a cutout of image.

    Extracts a box around (cx_init, cy_init), runs SEP source detection, and
    returns the centroid of the nearest detected source within search_radius px.
    Falls back to a flux-weighted centroid if SEP finds nothing.

    Returns (x, y) in full-frame 0-indexed pixel coordinates.
    """
    h, w   = image.shape
    box    = search_radius + 10
    y0 = max(0, int(cy_init) - box);  y1 = min(h, int(cy_init) + box + 1)
    x0 = max(0, int(cx_init) - box);  x1 = min(w, int(cx_init) + box + 1)
    cut = np.where(np.isfinite(image[y0:y1, x0:x1]),
                   image[y0:y1, x0:x1], 0.0).astype(np.float64)

    # Local background subtraction
    try:
        _bkg  = sep.Background(cut, bw=32, bh=32, fw=3, fh=3)
        cut_s = cut - _bkg.back()
        _rms  = _bkg.rms()
    except Exception:
        _, _med, _std = sigma_clipped_stats(cut[cut != 0], sigma=3.0)
        cut_s = cut - _med
        _rms  = np.full_like(cut, max(_std, 1.0))

    # SEP detection: look for the nearest source within search_radius
    try:
        _objs, _ = sep.extract(cut_s, thresh=3.0, err=_rms, minarea=5,
                               deblend_nthresh=32, deblend_cont=0.005)
        if len(_objs) > 0:
            _dist = np.hypot(_objs["x"] - (cx_init - x0),
                             _objs["y"] - (cy_init - y0))
            _near = np.argmin(_dist)
            if _dist[_near] <= search_radius:
                return (float(_objs["x"][_near] + x0),
                        float(_objs["y"][_near] + y0))
    except Exception:
        pass

    # Flux-weighted centroid fallback within inner half of search box
    _inner = max(5, search_radius // 2)
    _yy, _xx = np.mgrid[y0:y1, x0:x1]
    _mask = ((_yy - cy_init)**2 + (_xx - cx_init)**2) <= _inner**2
    _flux = np.clip(cut_s * _mask, 0, None)
    _tot  = _flux.sum()
    if _tot > 0:
        return (float((_xx * _flux).sum() / _tot),
                float((_yy * _flux).sum() / _tot))

    return float(cx_init), float(cy_init)


# ══════════════════════════════════════════════════════════════════════════════
# 8. PER-FRAME CHARACTERISATION
# ══════════════════════════════════════════════════════════════════════════════
def characterise_frame(data, cat, ref_wcs):
    """
    For one aligned frame: extract sources, fit PSF, estimate σ and F.
    Returns (psf_stamp, P_hat, F_j, sigma_j, fwhm_px, n_psf_stars).
    Falls back to analytic Gaussian PSF if stellar fitting fails.
    """
    objects, _  = extract_sources(data, detect_sigma=DETECT_SIGMA)
    sigma_j     = estimate_sigma(data)

    psf_stamp, fwhm_x, fwhm_y = fit_psf_moffat(
        data, objects,
        stamp_size=PSF_STAMP_SIZE,
        max_stars=MAX_PSF_STARS,
        isolation_fwhm=ISOLATION_FWHM,
    )

    if psf_stamp is None:
        # Analytic fallback: symmetric Moffat at 4 px FWHM, alpha=3.5
        fwhm_px   = 4.0
        alpha_fb  = 3.5
        gamma_fb  = fwhm_px / (2.0 * np.sqrt(2.0**(1.0 / alpha_fb) - 1.0))
        half      = PSF_STAMP_SIZE // 2
        yy, xx    = np.mgrid[-half:half + 1, -half:half + 1]
        psf_stamp = (1.0 + (xx**2 + yy**2) / gamma_fb**2)**(-alpha_fb)
        psf_stamp /= psf_stamp.sum()
        fwhm_x = fwhm_y = fwhm_px
        n_stars = 0
    else:
        n_stars = MAX_PSF_STARS

    fwhm_px = 0.5 * (fwhm_x + fwhm_y)
    P_hat   = psf_to_hat(psf_stamp, data.shape)
    F_j     = estimate_zeropoint(data, objects, cat, ref_wcs,
                                 aper_r=APER_RADIUS_PX,
                                 ks_min=TWOMASS_KS_MIN,
                                 ks_max=TWOMASS_KS_MAX)
    return psf_stamp, P_hat, F_j, sigma_j, fwhm_px, n_stars


# ══════════════════════════════════════════════════════════════════════════════
# 8. PHASE 1 — PROPER COADDITION REFERENCE  (ZOGY Eq. 22–24)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("PHASE 1 — Proper coaddition reference  (ZOGY Eq. 22–24)")
print("=" * 72)

# 2MASS catalogue (cached)
cat_2mass = query_2mass(config.TARGET_RA, config.TARGET_DEC,
                        radius_arcmin=TWOMASS_RADIUS)

# Reference-frame WCS: use OB1 frame 0 (defines the aligned pixel grid)
_ref_wcs = None
_ref_ob1_dir = ALIGNED_DIR / "GX339_Ks_Imaging_1"
_ref_files   = sorted(_ref_ob1_dir.glob("HAWKI.*_1_cal_aligned.fits"))
if _ref_files:
    with fits.open(_ref_files[0]) as _h:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _ref_wcs = WCS(_h[0].header)
    print(f"  Reference WCS from: {_ref_files[0].name}")
else:
    print("  WARNING: reference frame not found — WCS-based ZP disabled")

ref_frames = get_aligned_frames(REF_OBS, test_mode=TEST_MODE,
                                test_n=TEST_N_FRAMES)
print(f"  Reference frames  : {len(ref_frames)}")

# Determine image shape from first frame
with fits.open(ref_frames[0][1]) as hdul:
    shape    = hdul[0].data.shape
    ref_hdr0 = hdul[0].header.copy()

# Image-space accumulators for inverse-variance weighted mean reference.
# The Fourier-domain proper coaddition (Eq. 22) produces R in physical-flux
# units (ADU / F_j ~ 1e-8), not ADU.  Plugging that directly into Eq. 13
# breaks the unit balance between the two numerator terms — the static scene
# no longer cancels, causing D_std >> 1.  A 1/σ² weighted mean keeps R in
# ADU so cancellation is exact; S/N improvement is equivalent to first order.
#
# W_sum is a PER-PIXEL weight array (not a scalar).  Aligned frames have NaN
# at dither-edge footprints; if NaN→0 before accumulation and divided by a
# scalar total_wt, those pixels get R ≈ 0 instead of the true background.
# The ZOGY subtraction then fails to cancel the science frame's sky at those
# edges, flooding D with Fourier-ringing artefacts and inflating D_std >> 1.
# Per-pixel W_sum correctly averages only the frames that have valid data at
# each spatial location, eliminating the edge artefact.
R_sum      = np.zeros(shape, dtype=np.float64)      # Σ wj * Rj   (ADU)
W_sum      = np.zeros(shape, dtype=np.float64)      # Σ wj * valid_mask  (per-pixel)
PSF_sum    = np.zeros((PSF_STAMP_SIZE, PSF_STAMP_SIZE))  # Σ wj * Pj
total_wt   = 0.0                                    # Σ wj  (scalar, for PSF only)
Fj_wt_sum  = 0.0                                    # Σ wj * Fj  → Fr_eff

frame_stats = []   # per-frame quality records
t0 = time.time()

for k, (ob_name, fpath) in enumerate(ref_frames):
    with fits.open(fpath) as hdul:
        data   = hdul[0].data.astype(np.float64)

    psf_j, P_hat_j, F_j, sigma_j, fwhm_px, n_stars = \
        characterise_frame(data, cat_2mass, _ref_wcs)

    # Replace NaN (alignment footprint) with 0 before accumulating.
    valid_mask = np.isfinite(data)
    d_clean    = np.where(valid_mask, data, 0.0)

    # Accumulate F_j/σ_j²-weighted mean (image space, ADU).
    # W_sum accumulates weight only where valid_mask is True, so the per-pixel
    # division R_sum/W_sum gives the correct weighted mean at every pixel
    # regardless of how many frames have valid data there (dither footprint).
    w          = F_j / sigma_j**2
    R_sum     += w * d_clean
    W_sum     += w * valid_mask          # per-pixel weight
    PSF_sum   += w * psf_j
    total_wt  += w                       # scalar — used for PSF & F_r only
    Fj_wt_sum += w * F_j

    frame_stats.append({
        "ob": ob_name, "file": fpath.name,
        "F_j": F_j, "sigma_j": sigma_j,
        "fwhm_px": fwhm_px, "n_stars": n_stars,
    })

    elapsed = time.time() - t0
    rate    = (k + 1) / elapsed
    eta     = (len(ref_frames) - k - 1) / rate if rate > 0 else 0
    print(f"  [{k+1:>3}/{len(ref_frames)}] {ob_name}/{fpath.name}  "
          f"FWHM={fwhm_px:.2f}px  F={F_j:.3e}  σ={sigma_j:.2f}  "
          f"ETA={eta:.0f}s")

# ── Build R: per-pixel inverse-variance weighted mean (ADU) ──────────────────
print("\n  Computing reference image R ...")
R         = R_sum / (W_sum + EPSILON)              # per-pixel weighted mean

# Mark pixels with zero total weight as NaN (no frame contributed there)
R[W_sum < EPSILON] = np.nan

# ── Reference PSF: weighted mean of individual PSF stamps ────────────────────
half      = PSF_STAMP_SIZE // 2
P_R_stamp = PSF_sum / (total_wt + EPSILON)
P_R_stamp /= (P_R_stamp.sum() + EPSILON)

# ── Reference flux zero point: weighted mean F_j ─────────────────────────────
F_r_scalar = Fj_wt_sum / (total_wt + EPSILON)     # representative Fr (ADU/flux)

# ── Empirical reference noise from the central region of R ───────────────────
# The edge of R has fewer contributing frames (smaller W_sum) and is noisier.
# Estimating sigma_r from the central 512×512 px (where all frames contribute)
# avoids the inflated edge RMS and gives the correct σ_r for the ZOGY formula.
_cy, _cx  = shape[0] // 2, shape[1] // 2
_half_reg = 256                                    # 512×512 central box
R_central = R[_cy - _half_reg : _cy + _half_reg,
              _cx - _half_reg : _cx + _half_reg]
sigma_r   = estimate_sigma(R_central)

# ── P_R_hat for use in Phase 2 ───────────────────────────────────────────────
P_R_hat = psf_to_hat(P_R_stamp, shape)

# ── Save reference FITS + PSF stamp ──────────────────────────────────────────
hdr_ref = ref_hdr0.copy()
hdr_ref["REFNFRM"]  = (len(ref_frames), "frames in weighted mean stack")
hdr_ref["REFFR"]    = (float(F_r_scalar), "ZOGY F_r (weighted mean F_j)")
hdr_ref["REFSIGR"]  = (float(sigma_r), "ZOGY sigma_r (empirical bkg RMS, ADU)")
hdr_ref["REFTYPE"]  = "ZOGY_INVVAR_MEAN"
fits.PrimaryHDU(data=R.astype(np.float32), header=hdr_ref).writeto(
    REF_FITS, overwrite=True)
fits.PrimaryHDU(data=P_R_stamp.astype(np.float32)).writeto(
    REF_PSF, overwrite=True)

print(f"  R saved   : {REF_FITS.name}")
print(f"  P_R saved : {REF_PSF.name}")
print(f"  F_r       = {F_r_scalar:.4e}  (weighted mean F_j)")
print(f"  σ_r       = {sigma_r:.4f} ADU  (expect ~{float(np.median([s['sigma_j'] for s in frame_stats]))/np.sqrt(len(ref_frames)):.4f})")


# ══════════════════════════════════════════════════════════════════════════════
# 9. PHASE 2 — ZOGY SUBTRACTION  (Eq. 13–15)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("PHASE 2 — ZOGY subtraction  (Eq. 13–15)")
print("=" * 72)

# Pre-compute reference Fourier array once.
# NaN pixels (zero-weight edge footprint) are filled with the median of R so
# that they contribute no signal in the difference — avoiding Fourier ringing.
R_median  = float(np.nanmedian(R))
R_clean   = np.where(np.isfinite(R), R, R_median)
R_hat_ref = fft2(R_clean)

# Compute GX 339-4 pixel position once from the reference WCS.
# Seed position from config (visually confirmed; WCS carries ~31px residual in Y).
# Tight search radius (8px < FWHM) prevents the centroid jumping to a neighbour.
_seed_x, _seed_y = config.TARGET_PIX_SEED
print(f"  GX 339-4 pixel (seed):          ({_seed_x}, {_seed_y})")
_cx_refined, _cy_refined = find_source_centroid(R, _seed_x, _seed_y, search_radius=8)
TARGET_PIX_X = int(round(_cx_refined))
TARGET_PIX_Y = int(round(_cy_refined))
print(f"  GX 339-4 pixel (centroid in R): ({TARGET_PIX_X}, {TARGET_PIX_Y})")

# WCS still needed to project 2MASS catalogue positions onto the pixel grid.
_ref_wcs_full = WCS(fits.getheader(REF_FITS))

# ── Select photometric reference stars from 2MASS catalogue ───────────────────
# Pick bright isolated stars to: (a) verify S≈0 for non-variable sources,
# (b) provide calibration anchors for flux_diff → magnitude conversion.
N_REF_STARS   = 8          # max reference stars to track
MIN_SEP_TARGET = 80        # pixels — exclude stars too close to GX 339-4
ref_stars = []
if cat_2mass is not None:
    _ks_arr  = np.array(cat_2mass["Kmag"],    dtype=float)
    _ra_arr  = np.array(cat_2mass["RAJ2000"], dtype=float)
    _dec_arr = np.array(cat_2mass["DEJ2000"], dtype=float)
    # Convert all catalogue positions to pixel coords
    _px_all, _py_all = _ref_wcs_full.all_world2pix(_ra_arr, _dec_arr, 0)
    _on_chip = ((_px_all > 50) & (_px_all < shape[1] - 50) &
                (_py_all > 50) & (_py_all < shape[0] - 50))
    _sep_from_target = np.hypot(_px_all - TARGET_PIX_X, _py_all - TARGET_PIX_Y)
    _usable = (_on_chip &
               (_sep_from_target > MIN_SEP_TARGET) &
               (_ks_arr >= TWOMASS_KS_MIN) &
               (_ks_arr <= TWOMASS_KS_MAX))
    # Sort by brightness (brightest first) and take up to N_REF_STARS
    _idx_sorted = np.where(_usable)[0][np.argsort(_ks_arr[_usable])][:N_REF_STARS]
    for _i, _idx in enumerate(_idx_sorted):
        ref_stars.append({
            "id":   f"ref_{_i+1:02d}",
            "ra":   float(_ra_arr[_idx]),
            "dec":  float(_dec_arr[_idx]),
            "x":    float(_px_all[_idx]),
            "y":    float(_py_all[_idx]),
            "Kmag": float(_ks_arr[_idx]),
        })
    print(f"  Reference stars   : {len(ref_stars)} selected  "
          f"(Ks = {ref_stars[0]['Kmag']:.1f}–{ref_stars[-1]['Kmag']:.1f} mag)")
    # Save manifest so Stage 10 knows which stars were used
    _rs_manifest = QUALITY_DIR / "reference_stars.json"
    with open(_rs_manifest, "w") as _f:
        json.dump(ref_stars, _f, indent=2)

diff_frames = get_aligned_frames(REF_OBS, test_mode=TEST_MODE,
                                 test_n=TEST_N_FRAMES)
print(f"  Science frames: {len(diff_frames)}")

diff_stats = []

for k, (ob_name, fpath) in enumerate(diff_frames):
    with fits.open(fpath) as hdul:
        N_data = hdul[0].data.astype(np.float64)
        N_hdr  = hdul[0].header.copy()

    psf_n, P_n_hat, F_n, sigma_n, fwhm_n, _ = \
        characterise_frame(N_data, cat_2mass, _ref_wcs)

    # Photometrically normalise N to the reference calibration scale.
    # ZOGY Eq. 13 cancels static sources only when Fr and Fn represent the
    # same photometric scale. Without normalisation, Fr−Fn ≠ 0 leaves star
    # residuals ∝ S·(Fr−Fn) that overwhelm the background noise floor.
    # Scaling N by Fr/Fn brings both images onto a common scale, so the
    # source terms Fr·P̂r·N̂_norm − Fr·P̂n·R̂ cancel exactly; the noise
    # in the denominator is updated accordingly (σ_n → σ_n · Fr/Fn).
    phot_scale   = F_r_scalar / F_n        # maps N → reference flux scale
    N_scaled     = N_data * phot_scale
    sigma_n_zogy = sigma_n * phot_scale    # noise of the re-scaled frame

    N_median = float(np.nanmedian(N_scaled))
    N_clean  = np.where(np.isfinite(N_scaled), N_scaled, N_median)
    N_hat    = fft2(N_clean)

    # ── ZOGY Eq. 13  (with F_n_eff = F_r after photometric normalisation) ─────
    # D̂ = (Fr·P̂r·N̂_norm − Fr·P̂n·R̂) / sqrt(Fr²·σn_norm²·|P̂r|² + Fr²·σr²·|P̂n|²)
    num_D = F_r_scalar * P_R_hat * N_hat - F_r_scalar * P_n_hat * R_hat_ref
    den_D = np.sqrt(
        F_r_scalar**2 * sigma_n_zogy**2 * np.abs(P_R_hat)**2
        + F_r_scalar**2 * sigma_r**2    * np.abs(P_n_hat)**2
        + EPSILON
    )
    D_hat = num_D / den_D
    D     = np.real(ifft2(D_hat))

    # Mask pixels where the science frame had no valid data
    D[~np.isfinite(N_data)] = np.nan

    # ── ZOGY Eq. 15 — difference image flux zero point ────────────────────────
    F_D = (F_r_scalar * F_r_scalar) / np.sqrt(
        F_r_scalar**2 * sigma_n_zogy**2 + F_r_scalar**2 * sigma_r**2 + EPSILON
    )

    # ── ZOGY Eq. 14 — PSF of difference image ─────────────────────────────────
    # P̂_D = Fr·Fr·P̂r·P̂n / (F_D · den_D)  (both F terms = F_r after normalisation)
    P_D_hat   = (F_r_scalar * F_r_scalar * P_R_hat * P_n_hat) / (F_D * den_D + EPSILON)
    P_D_full  = np.real(ifft2(P_D_hat))
    P_D_full  = np.roll(np.roll(P_D_full, half, axis=0), half, axis=1)
    P_D_stamp = P_D_full[:PSF_STAMP_SIZE, :PSF_STAMP_SIZE].copy()
    s = P_D_stamp.sum()
    if abs(s) > EPSILON:
        P_D_stamp /= s

    # ── ZOGY Eq. 17 — Optimal S statistic (matched-filter detection image) ──────
    # S = IFFT(F_D · P̂_D* · D̂)
    # For a point source with differential flux Δf at position x0:
    #   E[S(x0)] = F_D · Δf  →  Δf = S(x0) / F_D
    # S is the optimal linear statistic for detecting point-source variability.
    S_hat    = F_D * np.conj(P_D_hat) * D_hat
    S        = np.real(ifft2(S_hat))
    S[~np.isfinite(N_data)] = np.nan

    # Extract S and differential flux at GX 339-4's fixed pixel position
    S_target  = float(S[TARGET_PIX_Y, TARGET_PIX_X])
    flux_diff = S_target / F_D   # differential flux in reference-scale ADU

    # Extract S at each reference star position (should be ≈ 0 for non-variables)
    ref_S = {}
    for rs in ref_stars:
        rx, ry = int(round(rs["x"])), int(round(rs["y"]))
        if 0 <= rx < shape[1] and 0 <= ry < shape[0] and np.isfinite(S[ry, rx]):
            ref_S[rs["id"]] = float(S[ry, rx])
        else:
            ref_S[rs["id"]] = float("nan")

    # Parse MJD from filename: HAWKI.2025-05-17T05_45_20.232_1_cal_aligned.fits
    _fname_ts  = fpath.stem.replace("HAWKI.", "").split("_1_")[0]  # "2025-05-17T05_45_20.232"
    _date, _t  = _fname_ts.split("T")
    _ts_isot   = f"{_date}T{_t.replace('_', ':')}"                # "2025-05-17T05:45:20.232"
    mjd        = float(Time(_ts_isot, format="isot", scale="utc").mjd)

    # ── Save D and P_D ────────────────────────────────────────────────────────
    out_ob  = DIFF_DIR / ob_name
    out_ob.mkdir(parents=True, exist_ok=True)
    d_name  = fpath.name.replace("_cal_aligned.fits", "_diff.fits")
    pd_name = fpath.name.replace("_cal_aligned.fits", "_diff_psf.fits")

    hdr_d = N_hdr.copy()
    hdr_d["DIFFTYPE"] = ("ZOGY_D",   "proper difference image")
    hdr_d["ZOGYFD"]   = (float(F_D),          "ZOGY F_D")
    hdr_d["ZOGYFR"]   = (float(F_r_scalar),   "ZOGY F_r")
    hdr_d["ZOGYFN"]   = (float(F_n),          "ZOGY F_n")
    hdr_d["ZOGYSIGR"] = (float(sigma_r),      "ZOGY sigma_r")
    hdr_d["ZOGYSIGN"] = (float(sigma_n),      "ZOGY sigma_n")
    hdr_d["FWHMN_PX"] = (float(fwhm_n),       "science FWHM (px)")

    fits.PrimaryHDU(data=D.astype(np.float32),       header=hdr_d).writeto(
        out_ob / d_name,  overwrite=True)
    fits.PrimaryHDU(data=P_D_stamp.astype(np.float32)).writeto(
        out_ob / pd_name, overwrite=True)

    # Use sigma-clipped stats: in a crowded galactic field, bright star residuals
    # dominate the unclipped std and give a misleading D_std >> 1.  After 3σ
    # clipping the background noise in D is close to the ZOGY theoretical ~1.0.
    D_fin = D[np.isfinite(D)]
    d_mean_clip, d_med_clip, d_std_clip = sigma_clipped_stats(D_fin, sigma=3.0, maxiters=5)
    d_std_raw  = float(np.nanstd(D))
    row = {
        "ob": ob_name, "file": fpath.name,
        "mjd": mjd,
        "D_mean": float(d_mean_clip), "D_std": float(d_std_clip),
        "D_std_raw": d_std_raw,
        "S_target": S_target, "flux_diff": flux_diff,
        "F_n": F_n, "F_D": float(F_D), "sigma_n": sigma_n, "fwhm_n": fwhm_n,
    }
    row.update(ref_S)   # adds S_ref_01, S_ref_02, … columns
    diff_stats.append(row)

    elapsed = time.time() - t0
    rate    = (k + 1) / elapsed if elapsed > 0 else 1
    eta     = (len(diff_frames) - k - 1) / rate
    print(f"  [{k+1:>3}/{len(diff_frames)}] {ob_name}/{d_name}  "
          f"D(std_clip)={d_std_clip:.3f}  D(std_raw)={d_std_raw:.3f}  ETA={eta:.0f}s")

print(f"\n  Difference images saved to: {DIFF_DIR}")


# ══════════════════════════════════════════════════════════════════════════════
# 10. PHASE 3 — REFERENCE QUALITY METRICS
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("PHASE 3 — Quality metrics")
print("=" * 72)

# ── CSVs ─────────────────────────────────────────────────────────────────────
stats_csv = QUALITY_DIR / "frame_stats.csv"
with open(stats_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["ob","file","F_j","sigma_j","fwhm_px","n_stars"])
    w.writeheader(); w.writerows(frame_stats)

_ref_ids = [rs["id"] for rs in ref_stars]

diff_csv = QUALITY_DIR / "diff_stats.csv"
with open(diff_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=[
        "ob","file","mjd",
        "D_mean","D_std","D_std_raw",
        "S_target","flux_diff","F_n","F_D","sigma_n","fwhm_n",
    ] + _ref_ids)
    w.writeheader(); w.writerows(diff_stats)

# ── Lightcurve CSV — one row per frame, ordered by MJD ───────────────────────
lc_csv = DIFF_DIR / "lightcurve_raw.csv"
lc_cols = ["mjd","ob","file","S_target","flux_diff","F_D","fwhm_n"] + _ref_ids
lc_rows = sorted(diff_stats, key=lambda r: r["mjd"])
with open(lc_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=lc_cols)
    w.writeheader()
    for r in lc_rows:
        w.writerow({k: r.get(k, float("nan")) for k in lc_cols})
print(f"  Lightcurve CSV  : {lc_csv}")

print(f"  Frame stats : {stats_csv.name}")
print(f"  Diff stats  : {diff_csv.name}")

# Convenience arrays
fwhms  = [s["fwhm_px"] for s in frame_stats if np.isfinite(s["fwhm_px"])]
Fjs    = [s["F_j"]     for s in frame_stats]
D_stds = [s["D_std"]   for s in diff_stats  if np.isfinite(s["D_std"])]
D_means= [s["D_mean"]  for s in diff_stats  if np.isfinite(s["D_mean"])]

# ── Plot 1: FWHM histogram ────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(fwhms, bins=25, color="steelblue", edgecolor="white", alpha=0.85)
med_fwhm = np.median(fwhms)
ax.axvline(med_fwhm, color="red", ls="--",
           label=f"Median {med_fwhm:.2f} px = {med_fwhm*PIXSCALE:.2f}\"")
ax.set_xlabel("FWHM (pixels)")
ax.set_ylabel("Frames")
ax.set_title("Per-frame PSF FWHM  (OBs 1–10,  quiescent)")
ax.legend()
fig.tight_layout()
fig.savefig(QUALITY_DIR / "q01_fwhm_histogram.png", dpi=150)
plt.close(fig)

# ── Plot 2: FWHM timeline (seeing variation) ──────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(fwhms, "o", ms=3, color="steelblue", alpha=0.7)
ax.axhline(med_fwhm, color="red", ls="--", alpha=0.8,
           label=f"Median {med_fwhm:.2f} px")
# OB boundaries
ob_labels = [s["ob"] for s in frame_stats]
boundaries = [i for i in range(1, len(ob_labels))
              if ob_labels[i] != ob_labels[i-1]]
for b in boundaries:
    ax.axvline(b - 0.5, color="grey", lw=0.7, alpha=0.5)
ax.set_xlabel("Frame index")
ax.set_ylabel("FWHM (pixels)")
ax.set_title("Seeing timeline — quiescent frames")
ax.legend()
fig.tight_layout()
fig.savefig(QUALITY_DIR / "q02_fwhm_timeline.png", dpi=150)
plt.close(fig)

# ── Plot 3: Flux zero-point timeline ─────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(Fjs, "o", ms=3, color="darkorange", alpha=0.8)
med_F = np.median(Fjs)
ax.axhline(med_F, color="red", ls="--", label=f"Median {med_F:.3e}")
for b in boundaries:
    ax.axvline(b - 0.5, color="grey", lw=0.7, alpha=0.5)
ax.set_xlabel("Frame index")
ax.set_ylabel("F_j  (ADU / 2MASS linear flux)")
ax.set_title("Flux zero point per frame  (2MASS Ks matching)")
ax.legend()
fig.tight_layout()
fig.savefig(QUALITY_DIR / "q03_zeropoint_timeline.png", dpi=150)
plt.close(fig)

# ── Plot 4: Reference PSF profile ─────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

im = axes[0].imshow(P_R_stamp, origin="lower", cmap="inferno",
                    interpolation="nearest")
axes[0].set_title("Reference PSF  $P_R$")
plt.colorbar(im, ax=axes[0])

# Radial profile from stamp centre
ctr  = PSF_STAMP_SIZE // 2
rads = np.arange(ctr + 1)
profile = [float(P_R_stamp[ctr, ctr + r]) for r in rads]
arcsec  = rads * PIXSCALE

# Estimate FWHM of P_R by half-max interpolation
peak     = P_R_stamp.max()
half_max = peak / 2.0
fwhm_R_px = np.nan
for r in range(1, len(profile)):
    if profile[r] < half_max:
        frac = (profile[r-1] - half_max) / (profile[r-1] - profile[r])
        fwhm_R_px = 2.0 * (r - 1 + frac)
        break

axes[1].semilogy(arcsec, profile, "o-", ms=4, color="steelblue")
if not np.isnan(fwhm_R_px):
    axes[1].axvline(fwhm_R_px / 2 * PIXSCALE, color="red", ls="--",
                    label=f"FWHM/2 = {fwhm_R_px/2*PIXSCALE:.2f}\"")
    axes[1].legend(fontsize=8)
axes[1].set_xlabel("Radius (arcsec)")
axes[1].set_ylabel("Normalised amplitude")
axes[1].set_title("Radial profile of $P_R$  (log)")
axes[1].grid(True, alpha=0.3)

axes[2].plot(P_R_stamp[ctr, :], label="Row")
axes[2].plot(P_R_stamp[:, ctr], label="Col")
axes[2].set_xlabel("Pixel from centre")
axes[2].set_ylabel("Amplitude")
axes[2].set_title("PSF cross-sections")
axes[2].legend()
axes[2].grid(True, alpha=0.3)

fwhm_str = (f"{fwhm_R_px:.2f} px = {fwhm_R_px*PIXSCALE:.2f}\""
            if not np.isnan(fwhm_R_px) else "N/A")
fig.suptitle(f"Reference PSF  $P_R$   FWHM ≈ {fwhm_str}  "
             f"({len(ref_frames)} frames coadded)", fontsize=10)
fig.tight_layout()
fig.savefig(QUALITY_DIR / "q04_reference_psf.png", dpi=150)
plt.close(fig)

# ── Plot 5: D-image pixel statistics ─────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].hist(D_stds, bins=25, color="steelblue", edgecolor="white", alpha=0.85)
axes[0].axvline(1.0, color="red", ls="--", lw=1.5, label="Expected = 1.0")
med_std = np.median(D_stds)
axes[0].axvline(med_std, color="orange", ls=":",
                label=f"Median = {med_std:.3f}")
axes[0].set_xlabel("3σ-clipped σ(D)  [should ≈ 1 for correct ZOGY]")
axes[0].set_ylabel("Frames")
axes[0].set_title("Difference image background σ  (sigma-clipped)")
axes[0].legend()

axes[1].hist(D_means, bins=25, color="coral", edgecolor="white", alpha=0.85)
axes[1].axvline(0.0, color="red", ls="--", lw=1.5, label="Expected = 0")
med_mean = np.median(D_means)
axes[1].axvline(med_mean, color="orange", ls=":",
                label=f"Median = {med_mean:.4f}")
axes[1].set_xlabel("mean(D)  [should ≈ 0]")
axes[1].set_ylabel("Frames")
axes[1].set_title("Difference image mean")
axes[1].legend()

fig.suptitle("ZOGY D-image pixel statistics  (diagnostic for subtraction quality)",
             fontsize=10)
fig.tight_layout()
fig.savefig(QUALITY_DIR / "q05_diff_pixel_stats.png", dpi=150)
plt.close(fig)

# ── Numerical summary ─────────────────────────────────────────────────────────
n_frames   = len(ref_frames)
med_sigma  = float(np.median([s["sigma_j"] for s in frame_stats]))
exp_sigma_r = med_sigma / np.sqrt(n_frames)

# Approximate 5σ depth in the reference
n_aper     = np.pi * APER_RADIUS_PX**2
depth_adu  = 5.0 * sigma_r * np.sqrt(n_aper)
depth_mag  = (-2.5 * np.log10(depth_adu / med_F)
              if (med_F > 0 and depth_adu > 0) else np.nan)

print(f"""
  ┌── Reference quality ─────────────────────────────────────────┐
  │  Frames in coadd         : {n_frames}
  │  Median per-frame σ      : {med_sigma:.4f} ADU
  │  Expected σ_r (σ/√N)    : {exp_sigma_r:.4f} ADU
  │  Measured σ_r            : {sigma_r:.4f} ADU
  │  σ_r ratio (meas/exp)   : {sigma_r/exp_sigma_r:.3f}   (1.0 = optimal)
  │  Median seeing FWHM      : {med_fwhm:.2f} px = {med_fwhm*PIXSCALE:.2f} arcsec
  │  Reference PSF FWHM      : {fwhm_R_px:.2f} px = {fwhm_R_px*PIXSCALE:.2f} arcsec
  │
  │  Median D std (3σ-clipped): {med_std:.3f}   (expect ~1.0; crowd-field OK if 0.7–1.3)
  │  Median D mean (clipped) : {med_mean:.4f}  (expect 0.0; subtract locally for phot)
  │
  │  5σ depth (r={APER_RADIUS_PX} px aper)  : {depth_adu:.1f} ADU ≈ Ks {depth_mag:.2f} mag
  └──────────────────────────────────────────────────────────────┘""")

# ── Summary JSON ──────────────────────────────────────────────────────────────
summary = {
    "n_ref_frames":        n_frames,
    "n_diff_frames":       len(diff_frames),
    "F_r":                 float(F_r_scalar),
    "sigma_r_measured":    float(sigma_r),
    "sigma_r_expected":    float(exp_sigma_r),
    "sigma_r_ratio":       float(sigma_r / exp_sigma_r),
    "fwhm_median_px":      float(med_fwhm),
    "fwhm_median_arcsec":  float(med_fwhm * PIXSCALE),
    "fwhm_R_px":           float(fwhm_R_px) if not np.isnan(fwhm_R_px) else None,
    "fwhm_R_arcsec":       float(fwhm_R_px * PIXSCALE) if not np.isnan(fwhm_R_px) else None,
    "D_std_median_clipped": float(med_std),
    "D_std_median_raw":    float(np.median([s["D_std_raw"] for s in diff_stats])),
    "D_mean_median":       float(med_mean),
    "depth_5sigma_adu":    float(depth_adu),
    "depth_5sigma_Ks_mag": float(depth_mag) if not np.isnan(depth_mag) else None,
    "ref_fits":            str(REF_FITS),
    "ref_psf_fits":        str(REF_PSF),
    "quality_dir":         str(QUALITY_DIR),
}
with open(QUALITY_DIR / "summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print(f"\n  Quality plots + CSV : {QUALITY_DIR}")
print(f"  Summary JSON        : {QUALITY_DIR / 'summary.json'}")

print("\n" + "=" * 72)
print("ZOGY PIPELINE COMPLETE")
print(f"  Reference image  : {REF_FITS}")
print(f"  Difference images: {DIFF_DIR}")
print(f"  Quality outputs  : {QUALITY_DIR}")
print("=" * 72)

# ══════════════════════════════════════════════════════════════════════════════
# EMAIL NOTIFICATION  (full run only)
# Set GMAIL_USER and GMAIL_APP_PW environment variables to enable.
# ══════════════════════════════════════════════════════════════════════════════
if not TEST_MODE:
    import os, smtplib
    from email.message import EmailMessage
    from datetime import datetime

    gmail_user = os.environ.get("GMAIL_USER")
    gmail_pw   = os.environ.get("GMAIL_APP_PW")

    if gmail_user and gmail_pw:
        try:
            elapsed_total = time.time() - t0
            hours, rem    = divmod(int(elapsed_total), 3600)
            mins, secs    = divmod(rem, 60)

            msg = EmailMessage()
            msg["Subject"] = "GX 339-4 Pipeline — ZOGY complete"
            msg["From"]    = gmail_user
            msg["To"]      = gmail_user
            msg.set_content(
                f"ZOGY pipeline finished at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
                f"  Reference frames    : {len(ref_frames)}\n"
                f"  Difference frames   : {len(diff_frames)}\n"
                f"  Runtime             : {hours}h {mins}m {secs}s\n\n"
                f"Quality summary:\n"
                f"  σ_r (measured)      : {sigma_r:.4f} ADU\n"
                f"  σ_r (expected σ/√N) : {exp_sigma_r:.4f} ADU\n"
                f"  σ_r ratio           : {sigma_r/exp_sigma_r:.3f}  (1.0 = optimal)\n"
                f"  Median FWHM         : {med_fwhm:.2f} px = {med_fwhm*PIXSCALE:.2f}\"\n"
                f"  Median D std        : {med_std:.3f}  (expect 1.0)\n"
                f"  Median D mean       : {med_mean:.4f}  (expect 0.0)\n\n"
                f"Outputs:\n"
                f"  Reference image     : {REF_FITS}\n"
                f"  Difference images   : {DIFF_DIR}\n"
                f"  Quality plots       : {QUALITY_DIR}\n"
            )
            with smtplib.SMTP_SSL("smtp.gmail.com", 465) as smtp:
                smtp.login(gmail_user, gmail_pw)
                smtp.send_message(msg)
            print(f"\nEmail notification sent to {gmail_user}")
        except Exception as e:
            print(f"\nEmail notification failed: {e}")
    else:
        print("\n(No email sent — set GMAIL_USER and GMAIL_APP_PW env vars to enable)")
print("=" * 72)
