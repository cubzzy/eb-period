"""Period finding on ZTF light curve CSVs (IRSA Time Series Tool or API downloads).

Usage:
    python3 ztf_period.py                      # every CSV in light-data/ next to this script
    python3 ztf_period.py light-data/B         # just one spectral class
    python3 ztf_period.py path/to/folder       # every CSV in that folder
    python3 ztf_period.py path/to/file.csv     # just that file
    options: --outdir some/folder

Each CSV must hold exactly one star (one oid); a file with several oids is skipped. Every
entry in METHODS is run on the same frequency grid so the results are directly comparable.

Subfolders count as spectral classes: light-data/B/*.csv, light-data/A/*.csv, ... give
spectral_class B, A, ... in results.csv (CSVs sitting directly in light-data get '').

For each star the periodogram peak P is compared with 2P (eclipsing binaries often peak
at half the true period): the fold with the tighter scatter wins, and both folds are drawn
in the PNG with the chosen one outlined in green. fold_ratio in results.csv is the 2P/P
scatter; values near 1 mean "no evidence for doubling", small values mean "clearly doubled".

Writes to output/<name of the folder you ran>/ (next to this script), e.g. output/B/ for
light-data/B, so each class keeps its own results:
    results.csv   one row per star: period + uncertainty (PRIMARY_METHOD) and every method's
    skipped.csv   stars that were skipped or failed, with the reason
    plots/        one PNG per star, named <csv name>_<oid>.png (the oid is also in results.csv)
"""
from astropy.timeseries import LombScargle
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
import argparse


matplotlib.use("Agg")  # no window popups

HERE = Path(__file__).resolve().parent

TIME_COL = "hjd"
PERIOD_RANGE = (0.2, 20.0)  # days, shared by every method
OVERSAMPLE = 10  # frequency grid spacing = 1 / (OVERSAMPLE * baseline)
REFINE_FACTOR = 100  # second p ss: grid this many times finer, +/- REFINE_STEPS coarse steps around the peak
REFINE_STEPS = 1
N_BOOT = 100  # bootstrap resamples for the period uncertainty
BOOT_STEPS = 5  # search window for each resample, in coarse steps (wider than REFINE_STEPS so scatter isn't clipped)

MIN_EPOCHS = 50  # skip stars with fewer good epochs
MIN_BASELINE_D = 100.0  # skip stars observed over less than this many days
PRIMARY_METHOD = "LS 2-term"  # its period and sigma go in the period_d / sigma_d columns
FOLD_BINS = 50  # phase bins when folding at P; folding at 2P uses twice as many (same time resolution)
DOUBLE_RATIO = 0.8  # use 2P only if its fold scatter is below this fraction of P's

# One entry per run: add, remove or edit entries to compare settings.
# kind="ls":  nterms = number of harmonics (2 helps eclipsing binaries)

METHODS = {
    "LS 1-term": dict(kind="ls", nterms=1),
    "LS 2-term": dict(kind="ls", nterms=2),
}


def load_epochs(path):
    """cuts out catflags != 0"""
    df = pd.read_csv(path, encoding="utf-8-sig")  # utf-8-sig tolerates an Excel BOM
    good = df[df["catflags"] == 0].dropna(subset=[TIME_COL, "mag", "magerr"]) # catflags == 32768 denotes a bad observation
    return good[good["magerr"] > 0]  # a zero error makes the weighted periodogram divide by zero


def frequency_grid(t):
    pmin, pmax = PERIOD_RANGE
    return np.arange(1 / pmax, 1 / pmin, 1 / (OVERSAMPLE * np.ptp(t)))


def power_spectrum(t, y, dy, freq, cfg):
    cfg = dict(cfg)
    kind = cfg.pop("kind")
    if kind == "ls":
        return LombScargle(t, y, dy, nterms=cfg.get("nterms", 1)).power(freq)
    raise ValueError(f"unknown method kind: {kind}")


def refine_peak(t, y, dy, freq, cfg, period, steps=REFINE_STEPS):
    """Re-search a narrow, much finer frequency window around a coarse-grid peak."""
    step = freq[1] - freq[0]
    centre = 1 / period
    fine = np.arange(centre - steps * step, centre + steps * step, step / REFINE_FACTOR)
    return 1 / fine[np.argmax(power_spectrum(t, y, dy, fine, cfg))]


def period_uncertainty(t, y, dy, freq, cfg, coarse, multiplier, rng):
    """Bootstrap: resample epochs with replacement, redo the fine search each time, and
    return the std of the final periods. The P-vs-2P choice is made once on the full data
    and applied to every resample, so a flip between P and 2P can't inflate the scatter."""
    n = len(t)
    finals = []
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        finals.append(multiplier * refine_peak(t[i], y[i], dy[i], freq, cfg, coarse, BOOT_STEPS))
    return np.std(finals, ddof=1)


def fold_variance(period, t, y, nbins):
    """Phase-dispersion measure: variance of the points inside phase bins, relative to the
    overall variance. Small means the fold at this period is tight."""
    bins = np.minimum(((t / period) % 1 * nbins).astype(int), nbins - 1)
    n = np.bincount(bins, minlength=nbins)
    s = np.bincount(bins, weights=y, minlength=nbins)
    ss = np.bincount(bins, weights=y * y, minlength=nbins)
    used = n > 0
    within = (ss - s**2 / np.maximum(n, 1))[used].sum() / max(len(y) - used.sum(), 1)
    return within / np.var(y, ddof=1)


def choose_multiplier(period, t, y):
    """Return (1 or 2, ratio). Eclipsing binaries often peak at half the true period because
    the two eclipses look alike, so fold at P and at 2P and take 2P only if it is clearly tighter."""
    ratio = fold_variance(2 * period, t, y, 2 * FOLD_BINS) / fold_variance(period, t, y, FOLD_BINS)
    return (2 if ratio < DOUBLE_RATIO else 1), ratio


def analyze_star(lc, rng):
    """Run every method on one light curve. Raises ValueError if it isn't usable."""
    t, y, dy = (lc[c].to_numpy() for c in (TIME_COL, "mag", "magerr"))
    if len(t) < MIN_EPOCHS:
        raise ValueError(f"only {len(t)} good epochs (need {MIN_EPOCHS})")
    if np.ptp(t) < MIN_BASELINE_D:
        raise ValueError(f"baseline only {np.ptp(t):.0f} d (need {MIN_BASELINE_D:.0f})")
    freq = frequency_grid(t)
    results = []
    for name, cfg in METHODS.items():
        power = power_spectrum(t, y, dy, freq, cfg)
        if not np.isfinite(power).all():  # argmax of NaNs would quietly return a bogus period
            raise ValueError(f"{name}: periodogram has non-finite values")
        coarse = 1 / freq[np.argmax(power)]
        fine = refine_peak(t, y, dy, freq, cfg, coarse)
        multiplier, ratio = choose_multiplier(fine, t, y)
        sigma = period_uncertainty(t, y, dy, freq, cfg, coarse, multiplier, rng)
        results.append({"name": name, "power": power, "fine": fine, "multiplier": multiplier, "ratio": ratio, "period": multiplier * fine, "sigma": sigma})
    return t, y, dy, freq, results


def result_row(source, spectral_class, lc, t, results):
    """One results.csv row: star info, the primary method's period, then every method's."""
    by_name = {r["name"]: r for r in results}
    primary = by_name[PRIMARY_METHOD]
    row = {
        "source_file": source,
        "spectral_class": spectral_class,
        "oid": lc["oid"].iloc[0],
        "ra": lc["ra"].median(),
        "dec": lc["dec"].median(),
        "filtercode": lc["filtercode"].iloc[0],
        "n_epochs": len(lc),
        "baseline_d": np.ptp(t),
        "period_d": primary["period"],
        "sigma_d": primary["sigma"],
        "doubled": primary["multiplier"] == 2,  # True if the periodogram peak was doubled
        "fold_ratio": primary["ratio"],  # 2P fold scatter / P fold scatter (below DOUBLE_RATIO = doubled)
        "peak_power": primary["power"].max(),
    }
    for r in results:
        row[f"{r['name']} period_d"] = r["period"]
        row[f"{r['name']} sigma_d"] = r["sigma"]
    return row


def fold_panel(ax, t, y, dy, period, label, chosen):
    phase = (t / period) % 1
    for cycle in (0, 1):  # two cycles so eclipses at the phase edge are visible
        ax.errorbar(phase + cycle, y, dy, fmt=".", ms=3, alpha=0.5, color="C0", elinewidth=0.5)
    ax.invert_yaxis()
    ax.set(xlabel="phase (2 cycles)", ylabel="mag", title=f"{label} = {period:.5f} d" + ("  [CHOSEN]" if chosen else ""))
    if chosen:
        for spine in ax.spines.values():
            spine.set_edgecolor("green")
            spine.set_linewidth(2.5)


def plot_results(out_path, title, t, y, dy, freq, results):
    """One row per method: periodogram, then the light curve folded at P and at 2P (chosen one outlined)."""
    fig, axes = plt.subplots(len(results), 3, figsize=(17, 3.2 * len(results)), squeeze=False)
    for (ax_pow, ax_p, ax_2p), r in zip(axes, results):
        ax_pow.plot(1 / freq, r["power"], lw=0.6)
        ax_pow.axvline(r["fine"], color="r", ls="--", lw=0.8, alpha=0.5, zorder=0)  # behind the data so it doesn't hide the peak
        ax_pow.set(xscale="log", xlabel="period (d)", ylabel="power", title=f"{r['name']}: peak at {r['fine']:.5f} d")
        fold_panel(ax_p, t, y, dy, r["fine"], "Fold at P", r["multiplier"] == 1)
        fold_panel(ax_2p, t, y, dy, 2 * r["fine"], f"Fold at 2P (scatter ratio {r['ratio']:.2f})", r["multiplier"] == 2)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)  # a batch run makes hundreds of figures; don't keep them in memory


def csv_files(target):
    """The CSVs to process: every .csv under a folder (subfolders included), or one file."""
    target = Path(target)
    return sorted(target.rglob("*.csv")) if target.is_dir() else [target]


def log_skip(skipped, path, oid, reason):
    skipped.append({"source_file": path.stem, "oid": oid, "reason": reason})
    print(f"{path.name}: SKIPPED ({reason})")


def load_star(path, skipped):
    """The light curve in one CSV (one oid per file), sorted by time. A file that can't be read,
    has no good epochs, or holds more than one oid is logged to skipped and gives None, so it
    can't stop a long run, vanish silently, or mix two light curves in one fit."""
    try:
        good = load_epochs(path)
    except Exception as err:
        log_skip(skipped, path, "", f"unreadable file: {type(err).__name__}: {err}")
        return None
    if good.empty:
        log_skip(skipped, path, "", "no epochs with catflags == 0")
    elif good["oid"].nunique() > 1:
        log_skip(skipped, path, "", f"{good['oid'].nunique()} oids in one file (expected 1)")
    else:
        return good.sort_values(TIME_COL)
    return None


def spectral_class_of(path):
    """The class folder a CSV sits in: light-data/B/x.csv -> 'B'. A CSV directly in light-data
    gives ''; a CSV outside light-data gives the name of its own folder."""
    path = path.resolve()
    try:
        parts = path.relative_to(HERE / "light-data").parts
    except ValueError:
        return path.parent.name
    return parts[0] if len(parts) > 1 else ""


def save_tables(outdir, rows, skipped):
    pd.DataFrame(rows).to_csv(outdir / "results.csv", index=False)
    pd.DataFrame(skipped, columns=["source_file", "oid", "reason"]).to_csv(outdir / "skipped.csv", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default=HERE / "light-data", help="a CSV file or a folder of CSVs")
    ap.add_argument("--outdir", type=Path, help="default: output/<name of the folder or file you ran>, so runs don't overwrite each other")
    args = ap.parse_args()

    target = Path(args.target).resolve()
    if args.outdir is None:
        args.outdir = HERE / "output" / (target.name if target.is_dir() else target.stem)
    plot_dir = args.outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    files = csv_files(args.target)
    if not files:
        raise SystemExit(f"No CSV files found in {args.target}")

    rows, skipped = [], []
    for path in files:
        lc = load_star(path, skipped)
        if lc is not None:
            oid = lc["oid"].iloc[0]
            try:
                t, y, dy, freq, results = analyze_star(lc, np.random.default_rng(int(oid)))  # seeded by oid: same star, same sigma, any run order
            except Exception as err:  # keep going: one bad star shouldn't stop a long run
                log_skip(skipped, path, oid, f"{type(err).__name__}: {err}")
            else:
                rows.append(result_row(path.stem, spectral_class_of(path), lc, t, results))
                plot_results(plot_dir / f"{path.stem}_{oid}.png", f"{path.stem}  oid {oid} ({lc['filtercode'].iloc[0]}, {len(lc)} epochs)", t, y, dy, freq, results)
                print(f"{path.name}: oid {oid}: P = {rows[-1]['period_d']:.6f} +/- {rows[-1]['sigma_d']:.6f} d" + ("  (doubled)" if rows[-1]["doubled"] else ""))
        save_tables(args.outdir, rows, skipped)  # rewritten after each file so a crash loses little

    print(f"\n{len(rows)} stars done, {len(skipped)} skipped. Results in {args.outdir}")


if __name__ == "__main__":
    main()