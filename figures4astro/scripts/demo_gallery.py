"""Render one example of every exoplot figure type from synthetic data.

    python demo_gallery.py OUT_DIR [--journal aas] [--family sans|serif]

The data are fake (seeded); the point is to check the style end to end.
Composition curves here are toy power laws, labelled as such -- real figures
must load published tables.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exoplot as xp  # noqa: E402

rng = np.random.default_rng(42)


def trapezoid(t, depth, t14=1.6, t12=0.25):
    """Toy transit (hours)."""
    x = np.abs(t)
    f = np.ones_like(t)
    full = x <= t14 / 2 - t12
    ing = (x > t14 / 2 - t12) & (x < t14 / 2)
    f[full] -= depth
    f[ing] -= depth * (t14 / 2 - x[ing]) / t12
    return f


def fig_white_light(out, J):
    fig, axes = xp.flux_residual_axes(ncols=2, size=xp.figsize(J, "double", 0.38))
    for j, (name, d) in enumerate([("Visit 1", 0.0062), ("Visit 2", 0.0061)]):
        t = np.linspace(-2.3, 2.3, 1800)
        model = trapezoid(t, d)
        flux = model + rng.normal(0, 4e-4, t.size)
        tm = np.linspace(-2.3, 2.3, 2000)
        xp.plot_lightcurve(axes[0, j], axes[1, j], t, flux, tm, trapezoid(tm, d), model,
                           err=np.full(t.size, 4e-4), bin_width=0.1, journal=J,
                           shade=[(-0.8, -0.55), (0.55, 0.8)])
        axes[0, j].text(0.03, 0.08, f"{name}\nNIRSpec G395H", transform=axes[0, j].transAxes,
                        fontsize=7)
        xp.panel_label(axes[0, j], "ab"[j], loc="lower right")
        if j:
            axes[0, j].set_ylabel(""); axes[1, j].set_ylabel("")
    axes[1, 0].set_ylim(-2200, 2200)
    return xp.finalize(fig, out / "01_white_light")


def fig_spec_stack(out, J):
    edges = xp.constant_r_edges(2.87, 5.17, 6)
    cen = np.sqrt(edges[:-1] * edges[1:])
    t = np.linspace(-2.3, 2.3, 900)
    depths = 0.0060 + 2.0e-4 * np.exp(-0.5 * ((cen - 4.3) / 0.1) ** 2)
    models = [trapezoid(t, d) for d in depths]
    fluxes = [m + rng.normal(0, 8e-4, t.size) for m in models]
    fig = plt.figure(figsize=xp.figsize(J, "double", 0.55))
    gs = fig.add_gridspec(1, 3, width_ratios=(1.2, 1.2, 0.5), wspace=0.16)
    a0, a1, a2 = (fig.add_subplot(gs[k]) for k in range(3))
    off = 0.0095
    cols = xp.stack_lightcurves(a0, t, fluxes, models, cen, offset=off, bin_width=0.08,
                                label_edges=list(zip(edges[:-1], edges[1:])))
    for k, (f, m) in enumerate(zip(fluxes, models)):
        xb, yb, _ = xp.bin_data(t, (f - m) * 1e6, width=0.08)
        a1.plot(xb, yb - k * 1500, "o", ms=2.2, mfc=cols[k], mec="k", mew=0.25)
        a1.axhline(-k * 1500, color="0.4", lw=0.6, ls="--")
        a2.hist((f - m) / 8e-4, bins=25, range=(-4, 4), histtype="step", color=cols[k], lw=1.0)
    g = np.linspace(-4, 4, 200)
    a2.plot(g, (f.size * 8 / 25) * np.exp(-0.5 * g ** 2) / np.sqrt(2 * np.pi), color="k", lw=0.8, ls="--")
    a2.set_yticks([]); a2.set_xlabel(r"Residual / $\sigma$")
    a1.set_yticks([]); a1.set_ylabel("Residuals + offset")
    for a in (a0, a1):
        a.set_xlabel(xp.label("Time from mid-transit", "h", J))
    for a, l in zip((a0, a1, a2), "abc"):
        xp.panel_label(a, l, inside=False)
    return xp.finalize(fig, out / "02_spectroscopic_stack")


def fig_transmission(out, J):
    wl = np.geomspace(0.6, 5.2, 3000)
    base = 15300 - 120 * np.log(wl)
    feat = (180 * np.exp(-0.5 * ((wl - 3.3) / 0.12) ** 2) + 140 * np.exp(-0.5 * ((wl - 4.3) / 0.07) ** 2)
            + 90 * np.exp(-0.5 * ((wl - 1.4) / 0.08) ** 2) + 80 * np.exp(-0.5 * ((wl - 2.3) / 0.15) ** 2))
    model = base + feat
    edges = xp.constant_r_edges(0.62, 5.15, 18)
    wb = np.sqrt(edges[:-1] * edges[1:])
    hw = np.vstack([wb - edges[:-1], edges[1:] - wb])
    err = 35 + 25 * (wb > 4)
    yb = np.interp(wb, wl, model) + rng.normal(0, err)
    fig, (ax, axc) = plt.subplots(2, 1, figsize=xp.figsize(J, "single", 1.25), sharex=True,
                                  gridspec_kw=dict(hspace=0.0, height_ratios=(1.4, 1)))
    xp.plot_spectrum(ax, wb, yb, err, wl_halfwidth=hw, journal=J, models={
        "Retrieved": dict(wl=wl, y=model, lo1=model - 40, hi1=model + 40, lo2=model - 80,
                          hi2=model + 80, color="#2c7bb6"),
        "Flat line": dict(wl=wl, y=np.full_like(wl, np.mean(yb)), color="0.5", ls=":")})
    ax.legend(loc="lower left", ncol=2)
    ax.set_xlabel("")
    # per-species contributions, filled in their fixed colours
    axc.fill_between(wl, base.min() - 50, base, color="0.75", lw=0, label="Haze / continuum")
    for sp, (c0, s, a) in {"H2O": (1.4, 0.08, 90), "CH4": (3.3, 0.12, 180), "CO2": (4.3, 0.07, 140)}.items():
        axc.fill_between(wl, base, base + a * np.exp(-0.5 * ((wl - c0) / s) ** 2),
                         color=xp.MOLECULE_COLORS[sp], alpha=0.6, lw=0, label=xp.mol(sp))
    axc.plot(wl, model, color="k", lw=0.8)
    axc.set_ylim(base.min() - 50, None)
    axc.legend(loc="upper right", ncol=2)
    xp.log_wavelength_axis(axc)
    axc.set_ylabel(xp.label("Transit depth", "ppm", J))
    axc.set_xlabel(xp.label("Wavelength", r"$\mu$m", J))
    xp.panel_label(ax, "a"); xp.panel_label(axc, "b")
    return xp.finalize(fig, out / "03_transmission_spectrum")


def fig_phase_curve(out, J):
    ph = np.linspace(0, 1, 2500)
    refl = 60 * np.clip(np.cos(2 * np.pi * (ph - 0.5)), 0, None) ** 1.2
    emis = 140 + 130 * np.cos(2 * np.pi * (ph - 0.5) + 0.25)
    tot = refl + emis
    ecl = np.abs(ph - 0.5) < 0.03
    tot_e = np.where(ecl, 0.0, tot)
    tr = (ph < 0.02) | (ph > 0.98)
    y = tot_e + rng.normal(0, 90, ph.size)
    fig, axes = xp.flux_residual_axes(size=xp.figsize(J, "single", 0.95))
    xp.plot_phase_curve(axes[0, 0], axes[1, 0], ph, y, np.full(ph.size, 90.0), ph, tot_e, tot_e,
                        components={"Emission": dict(y=emis, lo=emis - 15, hi=emis + 15, color=xp.DAYSIDE),
                                    "Reflection": dict(y=refl, lo=refl * 0.7, hi=refl * 1.3,
                                                       color=xp.REFLECTED)},
                        transit_mask=tr, journal=J)
    axes[0, 0].legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3)
    axes[0, 0].set_xlim(0, 1)
    axes[1, 0].set_ylim(-250, 250)
    return xp.finalize(fig, out / "04_phase_curve")


def fig_tp(out, J):
    P = np.logspace(-5, 2, 200)
    T = 1400 + 350 * np.tanh(np.log10(P) + 1.5)
    fig, (a, b) = plt.subplots(1, 2, figsize=xp.figsize(J, "single", 0.8), sharey=True,
                               gridspec_kw=dict(wspace=0.06))
    xp.plot_tp(a, T, P, T - 80, T + 80, T - 170, T + 170, color="k", label_="Retrieved", journal=J)
    xp.plot_tp(a, T + 120 * np.exp(-((np.log10(P) + 1) ** 2)), P, color=xp.DAYSIDE, ls="--",
               label_="GCM dayside", journal=J)
    for k, (sp, lo, hi) in enumerate([("H2O", 1e-3, 1e-1), ("CO2", 1e-4, 1e-2)]):
        xp.contribution_bar(a, 1050 + 40 * k, lo, hi, xp.MOLECULE_COLORS[sp])
    a.legend(loc="lower left")
    for sp, x0 in (("H2O", -3.3), ("CO", -3.6), ("CH4", -6.0)):
        prof = x0 - 0.0 * np.log10(P) - (1.5 if sp == "CH4" else 0) * (np.log10(P) < -2)
        b.plot(prof, P, color=xp.MOLECULE_COLORS[sp], lw=1.1, label=xp.mol(sp))
    b.set_xlabel(r"log$_{10}$ VMR"); b.set_xlim(-8, -2.5); b.legend(loc="lower left")
    xp.panel_label(a, "a", loc="upper right"); xp.panel_label(b, "b", loc="upper right")
    return xp.finalize(fig, out / "05_tp_profile")


def fig_posteriors(out, J):
    fig, axes = plt.subplots(2, 3, figsize=xp.figsize(J, "double", 0.42),
                             gridspec_kw=dict(hspace=0.55, wspace=0.08))
    specs = [("H2O", rng.normal(-3.2, 0.4, 20000), None), ("CO2", rng.normal(-4.6, 0.3, 20000), None),
             ("CH4", -6 - np.abs(rng.normal(0, 2, 20000)), "upper"),
             ("CO", rng.uniform(-12, -1, 20000) ** 1, "upper"),
             ("[M/H]", rng.normal(2.1, 0.25, 20000), None), ("C/O", rng.normal(0.45, 0.08, 20000), None)]
    for ax, (name, s, lim) in zip(axes.flat, specs):
        c = xp.MOLECULE_COLORS.get(name, "0.3")
        xp.plot_posterior(ax, s, color=c, limit=lim)
        ax.set_xlabel(rf"log$_{{10}}$ {xp.mol(name)}" if name in xp.MOLECULE_COLORS else name)
    return xp.finalize(fig, out / "06_posteriors")


def fig_rv(out, J):
    P = 11.19
    t = np.sort(rng.uniform(0, 600, 160))
    inst = np.where(rng.random(t.size) < 0.6, "NIRPS", "HARPS")
    err = np.where(inst == "NIRPS", 0.9, 1.5)
    K = 1.3
    rv = K * np.sin(2 * np.pi * t / P) + rng.normal(0, err)
    ph = ((t / P + 0.5) % 1) - 0.5
    mp = np.linspace(-0.5, 0.5, 300)
    mk = K * np.sin(2 * np.pi * mp)
    rvm = K * np.sin(2 * np.pi * t / P)
    draws = [(K + rng.normal(0, 0.15)) * np.sin(2 * np.pi * mp) for _ in range(40)]
    fig = plt.figure(figsize=xp.figsize(J, "double", 0.36))
    gs = fig.add_gridspec(2, 2, height_ratios=(3, 1), width_ratios=(1, 1.15), hspace=0, wspace=0.28)
    a, r = fig.add_subplot(gs[0, 0]), None
    r = fig.add_subplot(gs[1, 0], sharex=a)
    a.tick_params(labelbottom=False)
    xp.plot_rv_phased(a, r, ph, rv, err, inst, mp, mk, resid=rv - rvm, draws=draws, journal=J)
    a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=6)
    a.set_ylim(-6, 6); r.set_ylim(-5, 5)
    p = fig.add_subplot(gs[:, 1])
    period = np.geomspace(1.1, 300, 4000)
    power = 0.04 * rng.random(period.size) + 0.32 * np.exp(-0.5 * ((period - P) / 0.08) ** 2) \
        + 0.12 * np.exp(-0.5 * ((period - 79) / 3) ** 2)
    xp.plot_periodogram(p, period, power, fap={"10%": 0.12, "1%": 0.15, "0.1%": 0.18},
                        mark={r"$P_{\rm b}$": P, r"$P_{\rm rot}$": 79}, journal=J)
    p.legend(loc="upper left")
    xp.panel_label(a, "a", loc="upper right"); xp.panel_label(p, "b", loc="upper right")
    return xp.finalize(fig, out / "07_radial_velocity")


def fig_mass_radius(out, J):
    Mg = np.geomspace(0.5, 30, 50)
    curves = {  # TOY power laws for the demo only; load Zeng+2016/2019 tables in real work
        "100% Fe (toy)": dict(M=Mg, R=0.77 * Mg ** 0.27, color="#7f7f7f", ls="-"),
        "Earth-like (toy)": dict(M=Mg, R=Mg ** 0.27, color="#8c564b", ls="-"),
        "50% H$_2$O (toy)": dict(M=Mg, R=1.25 * Mg ** 0.27, color="#1f78b4", ls="--"),
        "1% H/He, 500 K (toy)": dict(M=Mg, R=1.8 * Mg ** 0.12, color="#fdae61", ls="-."),
    }
    bgM = 10 ** rng.uniform(-0.3, 1.5, 400); bgR = np.clip(bgM ** 0.4 * rng.lognormal(0, 0.18, 400), 0.6, 5)
    M = 10 ** rng.uniform(0.1, 1.2, 14); R = M ** 0.3 * rng.lognormal(0.2, 0.15, 14)
    Teq = rng.uniform(300, 1200, 14)
    fig, ax = plt.subplots(figsize=xp.figsize(J, "single", 0.85))
    sc = xp.plot_mass_radius(ax, M, R, 0.15 * M, 0.05 * R, c=Teq, curves=curves,
                             background=(bgM, bgR), highlight={"TOI-0000 b": (6.0, 2.4)}, journal=J)
    ax.set_xlim(0.5, 30); ax.set_ylim(0.6, 4.2)
    xp.colorbar(fig, sc, ax, label=xp.label(r"$T_{\rm eq}$", "K", J))
    ax.legend(loc="upper left", fontsize=6)
    return xp.finalize(fig, out / "08_mass_radius")


def fig_map(out, J):
    lon = np.linspace(-180, 180, 145); lat = np.linspace(-90, 90, 73)
    L, B = np.meshgrid(np.radians(lon), np.radians(lat))
    T = 900 + 900 * np.clip(np.cos(L - 0.35) * np.cos(B), 0, None) ** 0.6 + 60 * np.cos(2 * B)
    fig, ax = plt.subplots(figsize=xp.figsize(J, "single", 0.62))
    pm = xp.plot_map(ax, lon, lat, T, contours=[1200, 1500], journal=J)
    xp.colorbar(fig, pm, ax, label=xp.label("Temperature", "K", J))
    return xp.finalize(fig, out / "09_temperature_map")


def fig_stack(out, J):
    rows = []
    teqs = np.sort(rng.uniform(300, 1200, 7))
    for k, teq in enumerate(teqs):
        edges = np.concatenate([xp.constant_r_edges(2.87, 3.72, 25), xp.constant_r_edges(3.82, 5.17, 25)])
        wl = np.sqrt(edges[:-1] * edges[1:])
        wl = wl[(wl < 3.72) | (wl > 3.82)]
        y = rng.normal(0, 0.6, wl.size) + (1.2 * np.exp(-0.5 * ((wl - 4.3) / 0.08) ** 2) if k % 2 else 0)
        rows.append(dict(wl=wl, y=y, yerr=np.full(wl.size, 0.6), name=f"Planet {k + 1}", value=teq))
    fig, ax = plt.subplots(figsize=xp.figsize(J, "single", 1.3))
    xp.offset_stack(ax, rows, offset=4.0, label_x=5.22,
                    fmt=lambda r: f"{r['name']}\n{r['value']:.0f} K")
    ax.set_xlim(2.8, 5.95)
    ax.set_xlabel(xp.label("Wavelength", r"$\mu$m", J))
    ax.set_ylabel(r"Relative altitude in scale heights, $z(\lambda)/H_{\rm eq}$")
    return xp.finalize(fig, out / "10_population_stack")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--journal", default="aas", choices=list(xp.JOURNAL_WIDTHS))
    ap.add_argument("--family", default="sans", choices=["sans", "serif"])
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    xp.use_style(family=a.family)
    for f in (fig_white_light, fig_spec_stack, fig_transmission, fig_phase_curve, fig_tp,
              fig_posteriors, fig_rv, fig_mass_radius, fig_map, fig_stack):
        paths = f(out, a.journal)
        print(f.__name__, "->", ", ".join(p.name for p in paths))


if __name__ == "__main__":
    main()
