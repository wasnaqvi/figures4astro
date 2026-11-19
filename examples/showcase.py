"""Build the README gallery for figures4astro.

    python examples/showcase.py [--out gallery] [--formats png]

Data:
  * mass-radius: NASA Exoplanet Archive (pscomppars) + Zeng et al. 2016/2019
    composition tables, downloaded once into examples/data/ (not committed).
  * everything else is synthetic but physically shaped: limb-darkened transits
    integrated over the stellar disk, an analytic isothermal transmission model (Heng & Kitzmann 2017
    form) with molecular bands at their real wavelengths, a cosine phase curve
    with a geometric eclipse, Keplerian RVs with an astropy Lomb-Scargle periodogram.

Needs numpy, scipy, matplotlib, astropy, corner.
"""
from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "figures4astro" / "scripts"))
import exoplot as xp  # noqa: E402

rng = np.random.default_rng(2026)
DATA = HERE / "data"

ARCHIVE_URL = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
               "select+pl_name,pl_bmasse,pl_bmasseerr1,pl_bmasseerr2,pl_rade,pl_radeerr1,"
               "pl_radeerr2,pl_eqt,pl_bmassprov+from+pscomppars+where+pl_bmasse+is+not+null"
               "+and+pl_rade+is+not+null+and+pl_rade<4.5+and+pl_bmasse<40&format=csv")
ZENG_URL = "https://lweb.cfa.harvard.edu/~lzeng/tables/{}"
ZENG = {"Pure iron": ("massradiusFe.txt", "#7f7f7f", "-"),
        "Earth-like": ("massradiusEarthlikeRocky.txt", "#8c564b", "-"),
        "50% water": ("massradius_50percentH2O_300K_1mbar.txt", "#1f78b4", "--"),
        "100% water": ("massradius_100percentH2O_300K_1mbar.txt", "#6baed6", ":")}


def fetch(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        urllib.request.urlretrieve(url, path)
    return path


# ---------------------------------------------------------------------------
# physics helpers
# ---------------------------------------------------------------------------
BANDS = {  # band centres (um), relative strengths, log-widths
    "H2O": [(0.94, 0.3, 0.035), (1.14, 0.6, 0.04), (1.38, 1.0, 0.05), (1.87, 1.0, 0.06),
            (2.69, 1.6, 0.06)],
    "CH4": [(1.16, 0.4, 0.03), (1.67, 0.9, 0.035), (2.32, 1.2, 0.04), (3.31, 2.2, 0.035)],
    "CO2": [(2.01, 0.6, 0.02), (2.69, 0.8, 0.02), (4.27, 2.8, 0.022)],
}


def line_forest(wl, seed, depth=0.9, scale=0.004):
    """Correlated multiplicative 'line forest' so bands look like real cross-sections."""
    r = np.random.default_rng(seed)
    z = gaussian_filter1d(r.normal(0, 1, wl.size), sigma=max(1.0, scale / np.diff(np.log(wl)).mean()))
    return np.exp(depth * z / z.std())


def cross_section(wl, species):
    s = np.zeros_like(wl)
    for c, a, w in BANDS[species]:
        s += a * np.exp(-0.5 * ((np.log(wl) - np.log(c)) / w) ** 2)
    return 1e-3 * (s * line_forest(wl, sum(map(ord, species))) + 1e-7)


def transmission(wl, vmr, rp=2.6, rs=0.45, H_km=110.0, haze=8000.0, cloud=None):
    """Analytic isothermal transit depth (ppm) and per-species contributions."""
    R_earth, R_sun = 6371.0, 695700.0
    rp_km, rs_km = rp * R_earth, rs * R_sun
    kappa = {sp: x * cross_section(wl, sp) for sp, x in vmr.items()}
    ray = 1e-9 * haze * (wl / 1.0) ** -4
    tot = sum(kappa.values()) + ray
    def depth(k):
        z = H_km * np.log(np.maximum(k, 1e-14) / 1e-9)
        if cloud is not None:
            z = np.maximum(z, cloud * H_km)
        return 1e6 * ((rp_km + z) / rs_km) ** 2
    cont = depth(ray + 1e-12)
    parts = {sp: depth(ray + k) for sp, k in kappa.items()}
    return depth(tot), cont, parts


def _overlap_angle(r, z, p):
    """Half-angle of the arc of a ring of radius r (star) covered by a disk of radius p at distance z."""
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (r ** 2 + z ** 2 - p ** 2) / (2 * r * z)
    ang = np.arccos(np.clip(c, -1, 1))
    ang = np.where(r + z <= p, np.pi, ang)          # ring fully inside the planet disk
    ang = np.where((np.abs(r - z) >= p) & (r + z > p), 0.0, ang)
    return ang


def transit_flux(z, p, u1, u2, nr=1500):
    """Quadratic limb-darkened transit by direct integration over stellar radius."""
    r = (np.arange(nr) + 0.5) / nr
    mu = np.sqrt(1 - r ** 2)
    I = 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2
    total = np.sum(I * 2 * np.pi * r) / nr
    z = np.atleast_1d(z)[:, None]
    blocked = np.sum(I * 2 * _overlap_angle(r[None, :], z, p) * r, axis=1) / nr
    return 1 - blocked / total


def disk_overlap(z, p):
    """Fraction of the planet disk (radius p) hidden behind the star (radius 1): eclipse."""
    z = np.abs(np.atleast_1d(z))
    out = np.zeros_like(z, float)
    full = z <= 1 - p
    out[full] = 1.0
    part = (z > 1 - p) & (z < 1 + p)
    zz = z[part]
    k0 = np.arccos((p ** 2 + zz ** 2 - 1) / (2 * p * zz))
    k1 = np.arccos((1 - p ** 2 + zz ** 2) / (2 * zz))
    area = p ** 2 * k0 + k1 - 0.5 * np.sqrt(np.clip(4 * zz ** 2 - (1 + zz ** 2 - p ** 2) ** 2, 0, None))
    out[part] = area / (np.pi * p ** 2)
    return out


def sky_separation(phase, a, inc_deg):
    """Projected star-planet separation (stellar radii) for a circular orbit; phase 0 = transit."""
    ph = 2 * np.pi * phase
    return a * np.sqrt(np.sin(ph) ** 2 + (np.cos(np.radians(inc_deg)) * np.cos(ph)) ** 2)


def bin_spectrum(wl, y, edges):
    idx = np.digitize(wl, edges) - 1
    return np.array([y[idx == k].mean() for k in range(len(edges) - 1)])


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def fig_transmission(out, fmt):
    wl = np.geomspace(0.58, 5.3, 12000)
    vmr = {"H2O": 1.0, "CH4": 2.0, "CO2": 0.6}
    model, cont, parts = transmission(wl, vmr)
    off = 2900 - np.median(model)
    model, cont = model + off, cont + off
    parts = {k: v + off for k, v in parts.items()}
    # NIRISS SOSS 0.6-2.8 um at R~35, NIRSpec G395H NRS1/NRS2 at R~45, real detector gap
    e1 = xp.constant_r_edges(0.62, 2.80, 35)
    e2 = np.concatenate([xp.constant_r_edges(2.87, 3.71, 45), xp.constant_r_edges(3.82, 5.17, 45)])
    rows = []
    for edges, sig, name in ((e1, 22, "NIRISS SOSS"), (e2, 30, "NIRSpec G395H")):
        keep = np.ones(len(edges) - 1, bool)
        if name.startswith("NIRSpec"):
            c = np.sqrt(edges[:-1] * edges[1:]); keep = (c < 3.71) | (c > 3.82)
        wb = np.sqrt(edges[:-1] * edges[1:])[keep]
        hw = np.vstack([wb - edges[:-1][keep], edges[1:][keep] - wb])
        tru = bin_spectrum(wl, model, edges)[keep]
        err = sig * (1 + 0.6 * (wb > 4.5)) * rng.uniform(0.85, 1.15, wb.size)
        rows.append((name, wb, hw, tru + rng.normal(0, err), err))
    sm = gaussian_filter1d(model, 30)
    sig1, sig2 = 14 + 0.004 * (sm - sm.min()), 30 + 0.008 * (sm - sm.min())

    fig, (ax, axc) = plt.subplots(2, 1, figsize=xp.figsize("aas", "double", 0.56), sharex=True,
                                  gridspec_kw=dict(hspace=0, height_ratios=(1.55, 1)))
    xp.model_band(ax, wl, sm, sm - sig1, sm + sig1, sm - sig2, sm + sig2, color="#2c7bb6")
    ax.plot(wl, gaussian_filter1d(cont, 30) + 0 * wl, color="0.55", ls=":", lw=0.9)
    for (name, wb, hw, y, e), mk in zip(rows, ("o", "D")):
        xp.data_points(ax, wb, y, e, xerr=hw, fmt=mk, ms=3.6 if mk == "o" else 3.1, label=name)
    xp.log_wavelength_axis(ax)
    ax.set_ylabel(xp.label("Transit depth", "ppm"))
    handles = [Line2D([], [], color="#2c7bb6", lw=1.2, label=r"Retrieved model, 1$\sigma$/2$\sigma$"),
               Line2D([], [], color="0.55", ls=":", lw=0.9, label="Haze only"),
               *ax.get_legend_handles_labels()[0][-2:]]
    ax.legend(handles=handles, loc="upper left", ncol=2)
    ax.set_ylim(sm.min() - 160, sm.max() + 260)
    xp.panel_label(ax, "a", loc="upper right")
    # contributions, filled in fixed species colours
    base = gaussian_filter1d(cont, 30)
    axc.fill_between(wl, base.min() - 400, base, color="0.82", lw=0)
    for sp in ("H2O", "CH4", "CO2"):
        p = gaussian_filter1d(parts[sp], 30)
        axc.fill_between(wl, base, np.maximum(p, base), color=xp.MOLECULE_COLORS[sp], alpha=0.55, lw=0,
                         label=xp.mol(sp))
    axc.plot(wl, sm, color="k", lw=0.8)
    axc.legend(loc="upper left", ncol=3)
    axc.set_ylim(base.min() - 40, sm.max() + 120)
    xp.log_wavelength_axis(axc)
    axc.set_ylabel(xp.label("Depth", "ppm"))
    axc.set_xlabel(xp.label("Wavelength", r"$\mu$m"))
    xp.panel_label(axc, "b", loc="upper right")
    return xp.finalize(fig, out / "transmission_spectrum", formats=fmt)


def fig_light_curves(out, fmt):
    t = np.linspace(-2.6, 2.6, 520) / 24.0
    per, a_rs, inc = 3.2, 14.0, 89.3
    z = sky_separation(t / per, a_rs, inc)
    edges = np.geomspace(0.85, 5.1, 8)
    wc = np.sqrt(edges[:-1] * edges[1:])
    fig = plt.figure(figsize=xp.figsize("aas", "double", 0.62))
    gs = fig.add_gridspec(1, 3, width_ratios=(1.25, 1.0, 0.45), wspace=0.17)
    a0, a1, a2 = (fig.add_subplot(gs[k]) for k in range(3))
    cmap = matplotlib.colormaps[xp.CMAPS["wavelength"]]
    off, roff = 0.0105, 1600.0
    th = t * 24
    for k, w in enumerate(wc):
        col = cmap(k / (len(wc) - 1))
        rp = 0.118 + 0.0016 * np.sin(2.2 * np.log(w)) + 0.0012 * np.exp(-0.5 * ((w - 4.3) / 0.12) ** 2)
        m = transit_flux(z, rp, 0.42 * (1.0 / w) ** 0.6, 0.18)
        sig = 650e-6 * (1 + 0.25 * (w > 3.5))
        f = m + rng.normal(0, sig, t.size)
        xb, yb, eb = xp.bin_data(th, f, width=0.12)
        _, rb, _ = xp.bin_data(th, (f - m) * 1e6, width=0.12)
        a0.plot(th, m - k * off, color="k", lw=0.7, zorder=4)
        a0.plot(xb, yb - k * off, "o", ms=2.4, mfc=col, mec="k", mew=0.25, zorder=3)
        a0.text(-2.5, 1.0 - k * off + 0.0016, f"{edges[k]:.2f}–{edges[k + 1]:.2f} µm", fontsize=6.3,
                va="bottom")
        a1.axhline(-k * roff, color="0.55", lw=0.6, ls="--")
        a1.plot(xb, rb - k * roff, "o", ms=2.4, mfc=col, mec="k", mew=0.25)
        a2.hist((f - m) / sig, bins=22, range=(-3.8, 3.8), histtype="step", color=col, lw=0.9,
                density=True)
    g = np.linspace(-3.8, 3.8, 200)
    a2.plot(g, np.exp(-0.5 * g ** 2) / np.sqrt(2 * np.pi), "k--", lw=0.9)
    a0.set_ylabel("Relative flux + offset")
    a1.set_ylabel("Residuals + offset"); a1.set_yticks([])
    a1.yaxis.set_label_position("left")
    a2.set_yticks([]); a2.set_xlabel(r"Residual / $\sigma$")
    for a in (a0, a1):
        a.set_xlabel(xp.label("Time from mid-transit", "h"))
        a.set_xlim(-2.6, 2.6)
    a0.set_ylim(1 - (len(wc) - 1) * off - 0.019, 1.0045)
    a1.set_ylim(-(len(wc) - 1) * roff - 900, 900)
    for a, l in zip((a0, a1, a2), "abc"):
        xp.panel_label(a, l, inside=False)
    return xp.finalize(fig, out / "spectroscopic_light_curves", formats=fmt)


def fig_phase_curve_map(out, fmt):
    ph = np.linspace(-0.02, 1.02, 4200)
    z = sky_separation(ph, 3.8, 88.5)
    behind = np.cos(2 * np.pi * ph) < 0                       # planet on the far side
    ecl = 1 - np.where(behind, disk_overlap(z, 0.155), 0.0)   # visible fraction of the planet
    offset = np.radians(18)
    emis = 1350 + 1150 * np.cos(2 * np.pi * (ph - 0.5) + offset)
    refl = 260 * np.clip(np.cos(2 * np.pi * (ph - 0.5)), 0, None) ** 1.5
    tot = (emis + refl) * ecl
    transit = (np.abs(ph) < 0.045) | (np.abs(ph - 1) < 0.045)
    y = tot + rng.normal(0, 210, ph.size)
    fig = plt.figure(figsize=xp.figsize("aas", "double", 0.42))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.0, 1.05), height_ratios=(3, 1), hspace=0, wspace=0.14)
    a = fig.add_subplot(gs[0, 0]); r = fig.add_subplot(gs[1, 0], sharex=a)
    a.tick_params(labelbottom=False)
    xp.plot_phase_curve(a, r, ph, y, np.full(ph.size, 210.0), ph, tot, tot,
                        components={"Thermal": dict(y=emis * ecl, lo=(emis - 60) * ecl, hi=(emis + 60) * ecl,
                                                    color=xp.DAYSIDE),
                                    "Reflected": dict(y=refl * ecl, lo=0.75 * refl * ecl, hi=1.25 * refl * ecl,
                                                      color=xp.REFLECTED)},
                        transit_mask=transit, bin_n=60)
    a.set_xlim(0, 1); a.set_ylim(-150, 3150); r.set_ylim(-330, 330)
    a.legend(loc="upper left", fontsize=6.5)
    from matplotlib.ticker import MaxNLocator
    a.yaxis.set_major_locator(MaxNLocator(nbins=6, prune="lower"))
    xp.panel_label(a, "a", loc="upper right")
    # Mollweide temperature map, hotspot shifted east by the phase offset
    m = fig.add_subplot(gs[:, 1], projection="mollweide")
    lon = np.linspace(-np.pi, np.pi, 361); lat = np.linspace(-np.pi / 2, np.pi / 2, 181)
    L, B = np.meshgrid(lon, lat)
    mu = np.cos(L - offset) * np.cos(B)
    T = 1150 + 1400 * np.clip(mu, 0, None) ** 0.55 + 90 * np.cos(B) ** 4 * np.cos(L + 0.9)
    pm = m.pcolormesh(L, B, T, cmap=xp.CMAPS["temperature"], shading="auto", rasterized=True)
    m.contour(L, B, T, levels=[1600, 2000, 2400], colors="w", linewidths=0.45, alpha=0.7)
    for lt in (-np.pi / 2, np.pi / 2):
        m.plot([lt, lt], [-np.pi / 2, np.pi / 2], color="w", ls="--", lw=0.8)
    m.plot(offset, 0, marker="+", color="w", ms=7, mew=1.1)
    m.grid(True, color="w", lw=0.3, alpha=0.35)
    m.set_xticklabels([]); m.tick_params(labelsize=6.5, colors="0.25")
    cb = xp.colorbar(fig, pm, m, label=xp.label("Temperature", "K"), orientation="horizontal",
                     pad=0.06, fraction=0.06, aspect=30)
    xp.panel_label(m, "b", inside=False)
    return xp.finalize(fig, out / "phase_curve_map", formats=fmt)


def fig_mass_radius(out, fmt):
    import csv
    path = fetch(ARCHIVE_URL, DATA / "archive_mr.csv")
    rows = list(csv.DictReader(open(path)))
    def f(r, k):
        try:
            return float(r[k])
        except ValueError:
            return np.nan
    M = np.array([f(r, "pl_bmasse") for r in rows]); R = np.array([f(r, "pl_rade") for r in rows])
    Me = np.array([0.5 * (f(r, "pl_bmasseerr1") - f(r, "pl_bmasseerr2")) for r in rows])
    Re = np.array([0.5 * (f(r, "pl_radeerr1") - f(r, "pl_radeerr2")) for r in rows])
    Teq = np.array([f(r, "pl_eqt") for r in rows])
    names = np.array([r["pl_name"] for r in rows])
    measured = np.array([r["pl_bmassprov"] == "Mass" for r in rows])
    good = measured & (M / Me > 5) & (R / Re > 12) & np.isfinite(Teq) & (R < 4.2) & (M < 30)
    fig, ax = plt.subplots(figsize=xp.figsize("aas", "single", 0.95, width=4.6))
    ax.plot(M[measured & ~good], R[measured & ~good], ".", color="0.82", ms=2.2, mew=0, zorder=0)
    for name, (fn, col, ls) in ZENG.items():
        d = np.loadtxt(fetch(ZENG_URL.format(fn), DATA / fn))
        ax.plot(d[:, 0], d[:, 1], color=col, ls=ls, lw=1.1, label=name, zorder=1)
    ax.errorbar(M[good], R[good], Re[good], Me[good], fmt="none", ecolor="0.55", elinewidth=0.45, zorder=2)
    sc = ax.scatter(M[good], R[good], c=Teq[good], cmap=xp.CMAPS["teq"], vmin=200, vmax=2000, s=16,
                    edgecolors="k", linewidths=0.4, zorder=3)
    labels = {"GJ 1214 b": (3.2, 3.55), "K2-18 b": (16, 3.05), "TOI-270 d": (1.6, 2.55),
              "LHS 1140 b": (14, 1.25), "55 Cnc e": (17, 1.55), "TRAPPIST-1 e": (0.5, 1.35)}
    for nm, xy in labels.items():
        k = np.where(names == nm)[0]
        if k.size:
            k = k[0]
            ax.scatter(M[k], R[k], marker="*", s=95, c="gold", edgecolors="k", linewidths=0.6, zorder=5)
            ax.annotate(nm, (M[k], R[k]), xytext=xy, textcoords="data", fontsize=6.8, zorder=6,
                        ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.6", lw=0.4),
                        arrowprops=dict(arrowstyle="-", color="0.3", lw=0.6, shrinkA=0, shrinkB=4))
    ax.set_xscale("log"); ax.set_xlim(0.3, 30); ax.set_ylim(0.6, 4.2)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel(xp.label("Mass", r"$M_\oplus$")); ax.set_ylabel(xp.label("Radius", r"$R_\oplus$"))
    ax.legend(loc="upper left", title="Zeng et al. (2016, 2019)", title_fontsize=6.5, fontsize=6.5)
    xp.colorbar(fig, sc, ax, label=xp.label(r"$T_{\rm eq}$", "K"), extend="max")
    return xp.finalize(fig, out / "mass_radius", formats=fmt)


def fig_rv(out, fmt):
    from astropy.timeseries import LombScargle
    P, K, Prot = 5.36, 2.1, 37.0
    t = np.sort(np.concatenate([rng.uniform(0, 140, 70), rng.uniform(360, 520, 90)]))
    inst = np.where(rng.random(t.size) < 0.55, "NIRPS", "HARPS")
    err = np.where(inst == "NIRPS", 1.1, 1.6)
    act = 1.6 * np.sin(2 * np.pi * t / Prot + 0.4)
    kep = K * np.sin(2 * np.pi * t / P)
    rv = kep + act + rng.normal(0, err)
    fig = plt.figure(figsize=xp.figsize("aas", "double", 0.36))
    gs = fig.add_gridspec(2, 2, height_ratios=(3, 1), width_ratios=(1, 1.2), hspace=0, wspace=0.26)
    a = fig.add_subplot(gs[0, 0]); r = fig.add_subplot(gs[1, 0], sharex=a); a.tick_params(labelbottom=False)
    ph = ((t / P + 0.5) % 1) - 0.5
    mp = np.linspace(-0.5, 0.5, 400)
    draws = [(K + rng.normal(0, 0.18)) * np.sin(2 * np.pi * mp + rng.normal(0, 0.05)) for _ in range(60)]
    xp.plot_rv_phased(a, r, ph, rv - act, err, inst, mp, K * np.sin(2 * np.pi * mp),
                      resid=rv - act - kep, draws=draws, nbins=12)
    a.set_ylim(-7.5, 7.5); r.set_ylim(-6, 6)
    a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=6.5, handletextpad=0.3)
    xp.panel_label(a, "a", loc="upper right")
    gp = fig.add_gridspec(3, 2, width_ratios=(1, 1.2), hspace=0, wspace=0.26)
    freq = np.linspace(1 / 300, 1 / 1.2, 30000)
    stages = [("RVs", rv), ("Activity removed", rv - act), (r"Activity + planet b removed", rv - act - kep)]
    for k, (name, y) in enumerate(stages):
        q = fig.add_subplot(gp[k, 1])
        ls = LombScargle(t, y, err)
        pw = ls.power(freq)
        fap = ls.false_alarm_level([0.1, 0.01, 0.001])
        xp.plot_periodogram(q, 1 / freq, pw, fap={"10%": fap[0], "1%": fap[1], "0.1%": fap[2]},
                            mark={r"$P_{\rm b}$": P, r"$P_{\rm rot}$": Prot} if k == 0 else None)
        q.set_ylim(0, 0.62); q.set_yticks([0, 0.2, 0.4])
        q.text(0.985, 0.88, name, transform=q.transAxes, ha="right", va="top", fontsize=6.5)
        if k < 2:
            q.tick_params(labelbottom=False); q.set_xlabel("")
        if k != 1:
            q.set_ylabel("")
        if k == 0:
            q.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), fontsize=6.5, ncol=3,
                     handlelength=1.6, columnspacing=1.0)
            xp.panel_label(q, "b", inside=False)
    return xp.finalize(fig, out / "radial_velocity", formats=fmt)


def fig_tp_corner(out, fmt):
    import corner
    P = np.logspace(-6, 1, 250)
    lp = np.log10(P)
    T = 760 + 230 * (1 + np.tanh(1.3 * (lp + 2.2))) / 2 + 40 * np.exp(-((lp + 4.5) / 0.8) ** 2)
    fig = plt.figure(figsize=xp.figsize("aas", "double", 0.48))
    left, right = fig.subfigures(1, 2, width_ratios=(1.75, 1.3), wspace=0.04)
    a, b = left.subplots(1, 2, sharey=True, gridspec_kw=dict(width_ratios=(1, 0.75), wspace=0.08))
    w1 = 25 + 35 * np.exp(-((lp + 2) / 2.5) ** 2); w2 = 2.1 * w1
    xp.plot_tp(a, T, P, T - w1, T + w1, T - w2, T + w2, color="k", label_="Retrieved")
    xp.plot_tp(a, T + 60 * np.tanh(lp + 3) - 20, P, color=xp.DAYSIDE, ls="--", label_="RCE model")
    for x0, (sp, lo, hi) in zip((620, 650, 680), (("H2O", 3e-4, 3e-2), ("CH4", 1e-4, 1e-2), ("CO2", 3e-5, 1e-3))):
        xp.contribution_bar(a, x0, lo, hi, xp.MOLECULE_COLORS[sp])
    a.set_xlim(590, 1120); a.legend(loc="lower left", fontsize=6.5)
    xp.panel_label(a, "a", loc="upper right")
    prof = {"H2O": -2.2 - 1.8 * (lp < -3.3), "CH4": -2.0 + 0.0 * lp - 1.2 * (lp < -4.5) * (-4.5 - lp),
            "CO2": -3.3 + 0.0 * lp, "CO": -4.4 + 0.25 * np.clip(lp + 1, 0, None)}
    for sp, x in prof.items():
        b.plot(x, P, color=xp.MOLECULE_COLORS[sp], lw=1.2, label=xp.mol(sp))
    b.set_xlim(-6.5, -1.2); b.set_xlabel(r"log$_{10}$ VMR"); b.tick_params(labelleft=False)
    b.legend(loc="lower left", fontsize=6.5)
    xp.panel_label(b, "b", loc="upper right")
    # corner of the four parameters that carry the claim
    cov = np.array([[0.09, 0.05, 0.02, 2.0], [0.05, 0.06, 0.01, 1.2], [0.02, 0.01, 0.12, 0.3],
                    [2.0, 1.2, 0.3, 900]])
    s = rng.multivariate_normal([-2.2, -2.0, -3.3, 830], cov, 40000)
    axes = right.subplots(4, 4, gridspec_kw=dict(hspace=0.06, wspace=0.06))
    corner.corner(s, fig=right, labels=[r"log H$_2$O", r"log CH$_4$", r"log CO$_2$", r"$T$ (K)"],
                  color="#2c7bb6", bins=30, smooth=1.0, plot_datapoints=False, fill_contours=True,
                  levels=(0.393, 0.865), label_kwargs=dict(fontsize=6.5), max_n_ticks=3, labelpad=0.18,
                  hist_kwargs=dict(lw=0.9))
    for ax_ in axes.flat:
        ax_.tick_params(labelsize=5.5, length=2, which="both")
        ax_.minorticks_off()
    left.subplots_adjust(top=0.94, bottom=0.14)
    right.subplots_adjust(left=0.2, bottom=0.2, top=0.94, right=0.98)
    xp.panel_label(axes[0, 0], "c", inside=False)
    return xp.finalize(fig, out / "tp_profile_posteriors", formats=fmt)


def fig_flowchart(out, fmt):
    """A JWST time-series analysis, from raw ramps to atmospheric constraints."""
    fig, ax = xp.flow_canvas(size=xp.figsize("aas", "double", 0.5))
    W, H = 0.168, 0.17
    y1, y2, y3 = 0.70, 0.40, 0.11
    row1 = [("Raw ramps", "uncalibrated\nintegrations", "data"),
            ("Stage 1", "ramps to slopes;\njumps, 1/f noise", "process"),
            ("Stage 2", "flat, background,\nwavelength solution", "process"),
            ("Stage 3", "1D spectral\nextraction", "process"),
            ("Light curves", "white light +\nspectroscopic bins", "data")]
    xs1 = [0.095, 0.2925, 0.49, 0.6875, 0.885]
    b1 = [xp.flow_box(ax, x, y1, t, d, k, w=W, h=H) for x, (t, d, k) in zip(xs1, row1)]
    for a, b in zip(b1, b1[1:]):
        xp.flow_arrow(ax, a, b)
    fits = xp.flow_box(ax, 0.885, y2, "Light-curve fits", "transit + systematics\n+ red-noise GP", "fit", w=W, h=H)
    spec = xp.flow_box(ax, 0.6875, y2, "Transmission\nspectrum", None, "product", w=W, h=H)
    retr = xp.flow_box(ax, 0.49, y2, "Retrieval", "free chemistry or\nequilibrium grids", "fit", w=W, h=H)
    post = xp.flow_box(ax, 0.2925, y2, "Posteriors", "abundances, T-P,\nclouds, model odds", "product", w=W, h=H)
    sci = xp.flow_box(ax, 0.095, y2, "Interpretation", "metallicity, C/O,\nformation, climate", "product", w=W, h=H)
    xp.flow_arrow(ax, b1[-1], fits, "bottom", "top")
    for a, b in ((fits, spec), (spec, retr), (retr, post), (post, sci)):
        xp.flow_arrow(ax, a, b, "left", "right")
    # external inputs, dashed
    ext = [(0.915, "Limb darkening", "stellar models", fits), (0.745, "Ephemeris", "archive priors", fits),
           (0.57, "Opacities", "line lists, CIA", retr), (0.40, "Stellar spectrum", "spots, faculae", retr)]
    for x, t, d, tgt in ext:
        e = xp.flow_box(ax, x, y3, t, d, "external", w=0.15, h=0.13, dashed=True)
        xp.flow_arrow(ax, e, tgt, "top", "bottom", dashed=True, color="0.45")
    # feedback: diagnostics send the reduction back a step
    ax.annotate("", xy=(0.70, y1 - H / 2 - 0.004), xytext=(0.83, y2 + H / 2 + 0.004), zorder=4,
                arrowprops=dict(arrowstyle="-|>", color=xp.FLOW_KINDS["fit"][1], lw=0.9, ls="--",
                                mutation_scale=8, shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=-0.3"))
    ax.text(0.69, 0.535, "residual\nchecks", ha="right", va="center", linespacing=1.0,
            fontsize=6.3, color=xp.FLOW_KINDS["fit"][1], style="italic")
    xp.flow_legend(ax, ("data", "process", "fit", "product", "external"), loc=(0.5, 0.87))
    return xp.finalize(fig, out / "analysis_flowchart", formats=fmt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE.parent / "gallery"))
    ap.add_argument("--formats", nargs="+", default=["png"])
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    xp.use_style()
    plt.rcParams["savefig.dpi"] = 220
    figs = [fig_transmission, fig_light_curves, fig_phase_curve_map, fig_mass_radius, fig_rv,
            fig_tp_corner, fig_flowchart]
    for f in figs:
        if a.only and f.__name__ not in a.only:
            continue
        print(f.__name__, "->", ", ".join(str(p.name) for p in f(out, tuple(a.formats))))


if __name__ == "__main__":
    main()
