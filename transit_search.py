#!/usr/bin/env python3
"""
transit_search.py
=================

Search photometric light curves for periodic transit signals, vet the
candidates, and estimate planet sizes.

Pipeline for each light curve file:
    1. Load time and relative flux (e.g. AstroImageJ .dat output)
    2. Remove invalid (NaN, zero or negative) points
    3. Normalize and detrend with a low-order polynomial,
       then clip upward outliers (cosmic rays, bad frames)
    4. Optionally smooth with a Savitzky-Golay filter
    5. Run a Box Least Squares (BLS) search over periods and durations
    6. Vet the best signal: odd/even depths, secondary eclipse, number
       of transits observed
    7. Mask the signal and repeat to look for additional planets
    8. Save diagnostic plots and a summary table

Author: Javed Iqbal (https://github.com/Jiqbal-astro)
"""

import argparse
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # save figures without needing a display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.stats import sigma_clip
from astropy.timeseries import BoxLeastSquares
from scipy.signal import savgol_filter

R_SUN_IN_R_EARTH = 109.1
R_SUN_IN_R_JUP = 9.731


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load_light_curves(directory, extension=".dat"):
    """Read every file with the given extension in `directory`.

    Files are expected to be whitespace-separated, with time in the first
    column (e.g. BJD_TDB) and relative flux in the second column
    (e.g. rel_flux_T1). Lines starting with '#' are ignored.
    """
    light_curves = {}
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(extension):
            path = os.path.join(directory, filename)
            light_curves[filename] = pd.read_csv(path, sep=r"\s+", comment="#")
    return light_curves


# --------------------------------------------------------------------------
# Cleaning and detrending
# --------------------------------------------------------------------------
def preprocess_light_curve(data, poly_degree=3, sigma=3.0,
                           smooth=False, window_length=51, polyorder=2):
    """Clean, normalize, detrend and (optionally) smooth one light curve.

    Returns (time, flux) arrays, or (None, None) if too few points remain.
    """
    time = data.iloc[:, 0].to_numpy(dtype=float)
    flux = data.iloc[:, 1].to_numpy(dtype=float)

    # Drop NaNs, infinities and non-physical (zero/negative) flux values
    good = np.isfinite(time) & np.isfinite(flux) & (flux > 0)
    time, flux = time[good], flux[good]
    order = np.argsort(time)
    time, flux = time[order], flux[order]

    if len(time) < max(5, poly_degree + 2):
        return None, None

    # Normalize to a median of 1
    flux = flux / np.median(flux)

    # Detrend: divide out a low-order polynomial in time.
    # Time is centred first to keep the polynomial fit numerically stable.
    t0 = time.mean()
    coeffs = np.polyfit(time - t0, flux, deg=poly_degree)
    flux = flux / np.polyval(coeffs, time - t0)

    # Clip only UPWARD outliers (cosmic rays, bad frames). Clipping low
    # points too would remove the transits we are trying to find.
    clipped = sigma_clip(flux, sigma_upper=sigma, sigma_lower=np.inf)
    keep = ~np.ma.getmaskarray(clipped)
    time, flux = time[keep], flux[keep]

    # Optional Savitzky-Golay smoothing
    if smooth:
        if len(flux) < window_length:
            return None, None
        flux = savgol_filter(flux, window_length=window_length, polyorder=polyorder)

    return time, flux


def robust_scatter(flux):
    """Standard deviation estimated from the median absolute deviation."""
    return 1.4826 * np.median(np.abs(flux - np.median(flux)))


# --------------------------------------------------------------------------
# Transit search and vetting
# --------------------------------------------------------------------------
def search_transit(time, flux, periods, durations):
    """Run a BLS search and return the best signal, the periodogram and the model."""
    dy = np.full_like(flux, robust_scatter(flux))
    bls = BoxLeastSquares(time, flux, dy=dy)
    result = bls.power(periods, durations)

    i = int(np.argmax(result.power))
    best = {
        "period_days": float(result.period[i]),
        "t0": float(result.transit_time[i]),
        "duration_days": float(result.duration[i]),
        "depth": float(result.depth[i]),
        "depth_err": float(result.depth_err[i]),
        "depth_snr": float(result.depth_snr[i]),
        "power": float(result.power[i]),
    }
    model = bls.model(time, best["period_days"], best["duration_days"], best["t0"])
    stats = bls.compute_stats(best["period_days"], best["duration_days"], best["t0"])
    return bls, best, result, model, stats


def _significance(a, b):
    """Difference between two (value, error) pairs in units of sigma."""
    err = np.hypot(a[1], b[1])
    return float(abs(a[0] - b[0]) / err) if err > 0 else 0.0


def vet_candidate(best, stats, snr_threshold, flag_sigma=3.0):
    """Basic false-positive checks. Returns a dict of metrics and flags."""
    n_transits = int(np.sum(np.asarray(stats["per_transit_count"]) > 0))
    odd_even_sigma = _significance(stats["depth_odd"], stats["depth_even"])

    sec_depth, sec_err = stats["depth_phased"]
    secondary_sigma = float(sec_depth / sec_err) if sec_err > 0 else 0.0

    flags = []
    if best["depth_snr"] < snr_threshold:
        flags.append("low_snr")
    if n_transits < 2:
        flags.append("single_transit")
    if odd_even_sigma > flag_sigma:
        flags.append("odd_even_mismatch")
    if secondary_sigma > flag_sigma:
        flags.append("secondary_eclipse")

    return {
        "n_transits": n_transits,
        "depth_odd": float(stats["depth_odd"][0]),
        "depth_even": float(stats["depth_even"][0]),
        "odd_even_sigma": odd_even_sigma,
        "secondary_depth": float(sec_depth),
        "secondary_sigma": secondary_sigma,
        "flags": ";".join(flags) if flags else "pass",
    }


def planet_size(depth, rstar=None):
    """Radius ratio from the transit depth, and planet radius if R* is given."""
    ratio = float(np.sqrt(depth)) if depth > 0 else np.nan
    out = {"rp_rs": ratio}
    if rstar is not None:
        out["rp_earth"] = ratio * rstar * R_SUN_IN_R_EARTH
        out["rp_jup"] = ratio * rstar * R_SUN_IN_R_JUP
    return out


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------
def phase_fold(time, period, t0):
    """Time from mid-transit (days), folded on the period."""
    return (time - t0 + 0.5 * period) % period - 0.5 * period


def binned(x, y, n_bins):
    """Median-bin y in n_bins equal-width bins of x."""
    edges = np.linspace(x.min(), x.max(), n_bins + 1)
    idx = np.digitize(x, edges) - 1
    centres, values = [], []
    for k in range(n_bins):
        sel = idx == k
        if np.count_nonzero(sel) >= 3:
            centres.append(0.5 * (edges[k] + edges[k + 1]))
            values.append(np.median(y[sel]))
    return np.array(centres), np.array(values)


def plot_candidate(time, flux, model, result, best, vet, title, out_path):
    """Save a 2x2 diagnostic figure for one candidate."""
    P, t0, dur = best["period_days"], best["t0"], best["duration_days"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    ax_lc, ax_pg, ax_fold, ax_oe = axes.ravel()

    # Light curve with model
    offset = np.floor(time.min())
    ax_lc.plot(time - offset, flux, ".", ms=2, color="tab:blue", alpha=0.6, label="Flux")
    ax_lc.plot(time - offset, model, "-", lw=1.2, color="tab:red", label="BLS model")
    ax_lc.set_xlabel(f"Time (BJD_TDB − {offset:.0f})")
    ax_lc.set_ylabel("Normalized flux")
    ax_lc.set_title("Light curve")
    ax_lc.legend(loc="lower right", fontsize=8)

    # Periodogram
    ax_pg.plot(result.period, result.power, color="k", lw=0.7)
    ax_pg.axvline(P, color="tab:red", ls="--", lw=1, label=f"P = {P:.4f} d")
    for harmonic in (0.5, 2.0):
        if result.period.min() <= P * harmonic <= result.period.max():
            ax_pg.axvline(P * harmonic, color="tab:red", ls=":", lw=0.8, alpha=0.6)
    ax_pg.set_xlabel("Period (days)")
    ax_pg.set_ylabel("BLS power")
    ax_pg.set_title("BLS periodogram (dotted: P/2, 2P)")
    ax_pg.legend(fontsize=8)

    # Phase-folded transit, zoomed in
    phase = phase_fold(time, P, t0) * 24.0  # hours
    window = 3 * dur * 24.0
    sel = np.abs(phase) < window
    ax_fold.plot(phase[sel], flux[sel], ".", ms=2, color="0.6", label="Data")
    if np.count_nonzero(sel) > 30:
        bx, by = binned(phase[sel], flux[sel], 40)
        ax_fold.plot(bx, by, "o", ms=4, color="tab:blue", label="Binned")
    order = np.argsort(phase[sel])
    ax_fold.plot(phase[sel][order], model[sel][order], "-", color="tab:red", lw=1.5, label="Model")
    ax_fold.set_xlabel("Hours from mid-transit")
    ax_fold.set_ylabel("Normalized flux")
    ax_fold.set_title(f"Phase-folded, depth = {best['depth'] * 1e3:.2f} ppt, SNR = {best['depth_snr']:.1f}")
    ax_fold.legend(fontsize=8)

    # Odd vs even transits
    epoch = np.round((time - t0) / P).astype(int)
    for parity, colour, label in ((1, "tab:orange", "Odd"), (0, "tab:green", "Even")):
        s = sel & (epoch % 2 == parity)
        if np.count_nonzero(s) > 30:
            bx, by = binned(phase[s], flux[s], 30)
            ax_oe.plot(bx, by, "o-", ms=3, lw=1, color=colour, label=label)
    ax_oe.set_xlabel("Hours from mid-transit")
    ax_oe.set_ylabel("Normalized flux")
    ax_oe.set_title(f"Odd vs even, difference = {vet['odd_even_sigma']:.1f}σ")
    ax_oe.legend(fontsize=8)

    for ax in axes.ravel():
        ax.grid(alpha=0.3)
    fig.suptitle(f"{title}  |  vetting: {vet['flags']}", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def analyse_light_curve(name, time, flux, args, out_dir):
    """Iteratively search one light curve for up to args.max_planets signals."""
    periods = np.linspace(args.min_period, args.max_period, args.n_periods)
    durations = np.array([d for d in args.durations if d < args.min_period])
    if durations.size == 0:
        raise ValueError("All trial durations must be shorter than --min-period")

    rows = []
    remaining = np.ones_like(time, dtype=bool)
    for k in range(1, args.max_planets + 1):
        t, f = time[remaining], flux[remaining]
        if len(t) < 50:
            break

        bls, best, result, model, stats = search_transit(t, f, periods, durations)
        vet = vet_candidate(best, stats, args.snr_threshold)
        if best["depth_snr"] < args.snr_threshold or best["depth"] <= 0:
            if k == 1:
                print(f"  no signal above SNR {args.snr_threshold} "
                      f"(best: P = {best['period_days']:.4f} d, SNR = {best['depth_snr']:.1f})")
            break

        size = planet_size(best["depth"], args.rstar)
        label = f"{Path(name).stem}_c{k}"
        msg = (f"  candidate {k}: P = {best['period_days']:.4f} d, "
               f"depth = {best['depth'] * 1e3:.2f} ppt, SNR = {best['depth_snr']:.1f}, "
               f"Rp/R* = {size['rp_rs']:.3f}")
        if "rp_earth" in size:
            msg += f", Rp = {size['rp_earth']:.2f} R_earth"
        print(msg + f"  [{vet['flags']}]")

        plot_candidate(t, f, model, result, best, vet,
                       title=f"{name}  candidate {k}",
                       out_path=out_dir / f"{label}.png")
        rows.append({"file": name, "candidate": k, **best, **size, **vet})

        # Mask this signal's transits before searching for the next one
        in_transit = bls.transit_mask(time, best["period_days"],
                                      args.mask_factor * best["duration_days"], best["t0"])
        remaining &= ~in_transit

    return rows


def process_directory(args):
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    light_curves = load_light_curves(args.input, args.extension)
    if not light_curves:
        print(f"No '{args.extension}' files found in {args.input}")
        return

    rows = []
    for name, data in light_curves.items():
        print(f"Processing {name} ...")
        time, flux = preprocess_light_curve(
            data,
            poly_degree=args.poly_degree,
            sigma=args.sigma,
            smooth=args.smooth,
            window_length=args.window,
        )
        if time is None:
            print("  skipped: not enough valid data points")
            continue
        rows.extend(analyse_light_curve(name, time, flux, args, out_dir))

    if rows:
        summary = pd.DataFrame(rows)
        summary.to_csv(out_dir / "bls_summary.csv", index=False)
        print(f"\n{len(rows)} candidate(s) saved to {out_dir}/")
    else:
        print("\nNo candidates found.")


def parse_args():
    p = argparse.ArgumentParser(
        description="Search light curves for transit signals with Box Least Squares, "
                    "vet the candidates and estimate planet sizes.")
    p.add_argument("input", help="directory containing light curve files")
    p.add_argument("-o", "--output", default="results", help="output directory (default: results)")
    p.add_argument("--extension", default=".dat", help="light curve file extension (default: .dat)")

    s = p.add_argument_group("search")
    s.add_argument("--min-period", type=float, default=0.5, help="minimum trial period, days (default: 0.5)")
    s.add_argument("--max-period", type=float, default=10.0, help="maximum trial period, days (default: 10)")
    s.add_argument("--n-periods", type=int, default=10000, help="number of trial periods (default: 10000)")
    s.add_argument("--durations", type=float, nargs="+", default=[0.04, 0.06, 0.08, 0.12, 0.16, 0.2],
                   help="trial transit durations, days (default: 0.04 0.06 0.08 0.12 0.16 0.2)")
    s.add_argument("--max-planets", type=int, default=3, help="maximum signals to search for per file (default: 3)")
    s.add_argument("--snr-threshold", type=float, default=7.0, help="minimum depth SNR to report a candidate (default: 7)")
    s.add_argument("--mask-factor", type=float, default=1.5,
                   help="mask width, in transit durations, when removing a found signal (default: 1.5)")
    s.add_argument("--rstar", type=float, default=None, help="stellar radius in solar radii, to report planet radius")

    c = p.add_argument_group("cleaning")
    c.add_argument("--poly-degree", type=int, default=3, help="detrending polynomial degree (default: 3)")
    c.add_argument("--sigma", type=float, default=3.0, help="upper sigma-clipping threshold (default: 3)")
    c.add_argument("--smooth", action="store_true", help="apply Savitzky-Golay smoothing before the search")
    c.add_argument("--window", type=int, default=51, help="Savitzky-Golay window length, odd (default: 51)")
    return p.parse_args()


if __name__ == "__main__":
    process_directory(parse_args())
