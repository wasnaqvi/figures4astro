# exoplot API

`scripts/exoplot.py` — numpy + matplotlib only (>= 3.6). Copy it into the
project, then:

```python
import matplotlib
matplotlib.use("Agg")                      # headless / scripts
import matplotlib.pyplot as plt
import exoplot as xp
xp.use_style()                             # once; family="serif" or usetex=True optional
```

## Style and geometry

| call | does |
|---|---|
| `use_style(fontsize=8, family="sans", usetex=False)` | apply rcParams: closed box, inward ticks on four sides, minor ticks, no caps, TrueType PDF, `axes.ymargin` 0.08 |
| `style_rc(...)` | the same dict, for `with mpl.rc_context(xp.style_rc()):` |
| `figsize(journal="aas", span="single", aspect=0.75, width=None)` | (w, h) inches at printed size. `span` is `"single"`, `"double"`, or a fraction of full width |
| `label(name, unit, journal)` | `"Wavelength (μm)"` or `"Wavelength [μm]"` per journal |
| `mol("CO2")` | `CO$_2$` (mathtext) |
| `panel_label(ax, "a", style="nature"/"paren"/"boxed", loc="upper left", inside=True)` | panel tag |
| `colorbar(fig, mappable, ax, label=None)` | colorbar with inward ticks and a thin outline |
| `finalize(fig, stem, formats=("pdf","png","svg"), dpi=300)` | align y labels, save each format, close, return paths |

## Constants

- `JOURNAL_WIDTHS`: `{journal: (single_in, full_in)}` for aas, aanda, mnras, nature.
- `MOLECULE_COLORS`, `MOLECULE_TEX`; `INSTRUMENT_COLORS`, `INSTRUMENT_MARKERS`.
- `DAYSIDE`, `NIGHTSIDE`, `REFLECTED`, `MODEL_GREY`, `RAW_GREY`.
- `LAYER_COLORS`: interior layers (H-He sand, H₂O blue, silicate brown, Fe grey).
- `CMAPS`: `temperature, teq, wavelength, phase, diverging, cloud, density, teq_stack`.
- `FLOW_KINDS`: face / edge colours and legend labels for flow-chart boxes.

## Data helpers

| call | returns |
|---|---|
| `bin_data(x, y, yerr=None, width=None, nbins=None, edges=None)` | `(xb, yb, eb)`. Inverse-variance mean when `yerr` is given, else mean and standard error. Empty bins dropped |
| `constant_r_edges(wmin, wmax, R)` | wavelength bin edges at constant λ/Δλ |

## Marks

| call | draws |
|---|---|
| `raw_and_binned(ax, x, y, yerr=None, bin_width=None, nbins=None, ...)` | faint rasterized raw dots (excluded from axis limits) plus open-circle bins. Returns the bins |
| `data_points(ax, x, y, yerr=None, xerr=None, color="k", open_=True)` | open-circle error bars, no caps |
| `model_band(ax, x, median, lo1, hi1, lo2, hi2, color, label)` | median line with 1σ / 2σ same-hue envelopes |
| `zero_line(ax, y=0)` | dashed grey reference line |
| `upper_limit(ax, x, y, color="k", frac=0.08)` | bar plus downward arrow. Call after the y-limits are final |
| `log_wavelength_axis(ax, ticks=(0.6,...,5), lim=None)` | log x with plain-number ticks |
| `shade_bands(ax, {"H2O": [(1.3, 1.5)], ...})` | pale vertical spans in species colours, labelled along the top |
| `contribution_bar(ax, x, p_lo, p_hi, color)` | vertical bar showing the pressures a band probes |

## Figure builders

| call | figure |
|---|---|
| `flux_residual_axes(fig=None, ncols=1, height_ratios=(3,1), size=None)` | `(fig, axes[2, ncols])`, flux over residual, zero gap, colliding tick labels pruned |
| `plot_lightcurve(ax_flux, ax_res, t, flux, model_t, model_f, model_at_data, err=None, bin_width=None, shade=None)` | white light: raw + bins + model, residuals in ppm with the rms printed, grey spans for `shade` |
| `stack_lightcurves(ax, t, fluxes, models, wavelengths, offset, bin_width, label_edges)` | spectroscopic channels offset downward, blue at top. Returns colours for matching residual panels |
| `plot_spectrum(ax, wl, depth, err, wl_halfwidth=None, models={...}, unit="ppm", log_x=True)` | transmission or emission spectrum. A `models` entry with `lo1`/`hi1`/... is drawn as a band, otherwise as a line (`ls` respected) |
| `plot_phase_curve(ax, ax_res, phase, fp, err, model_phase, model_fp, model_at_data, components=None, bin_n=40, transit_mask=None)` | Fp/F⋆ with component curves (emission, reflection) and bands, mid-eclipse line, residuals |
| `plot_tp(ax, T, P_bar, lo1, hi1, lo2, hi2, color, label_, ls)` | T-P profile, log P increasing downward, envelopes in T |
| `plot_posterior(ax, samples, color, label_, bins=40, filled=True, limit=None)` | filled step histogram, y ticks hidden. `limit="upper"/"lower"` adds a 95% arrow |
| `plot_rv_phased(ax, ax_res, phase, rv, err, instrument, model_phase, model_rv, resid=None, draws=None, nbins=10)` | per-instrument markers and colours, Keplerian, posterior draws, large open phase bins, O − C panel |
| `plot_periodogram(ax, period, power, fap={"1%": 0.15}, mark={"P_b": 11.19})` | GLS on log period with FAP lines and marked periods |
| `plot_mass_radius(ax, M, R, Merr, Rerr, c=Teq, curves={...}, highlight={...}, background=(M_all, R_all))` | grey archive population, coloured sample, composition curves, starred targets. Returns the scatter for a colorbar |
| `plot_map(ax, lon, lat, field, cmap, terminators=True, substellar=True, contours=None)` | equirectangular map with terminators at ±90° |
| `offset_stack(ax, rows, offset, cmap=CMAPS["teq_stack"], by_rank=True, label_x=None, fmt=None)` | Sing-style population stack, first row at the top, points only |
| `flow_canvas(fig=None, size=None)` | `(fig, ax)`: a frameless axis in unit coordinates for flow charts |
| `flow_box(ax, x, y, title, detail=None, kind="process", w=0.16, h=0.13, dashed=False)` | rounded box with a bold title and a smaller detail line. `kind` is data / process / fit / product / external. Returns anchors |
| `flow_arrow(ax, a, b, side_a="right", side_b="left", dashed=False, rad=0, label=None)` | arrow between box anchors; `rad` bends it, dashed for optional inputs |
| `flow_legend(ax, kinds, loc=(0.5, 0.02))` | swatch legend of box kinds |

## Contract notes

- Every builder labels its own axes using the journal's brackets. Pass
  `journal=` when it is not AAS.
- Builders never call `plt.show()` or save. `finalize` saves.
- `plot_mass_radius` takes composition curves as data. Load the published
  tables (Zeng et al. 2016/2019, Lopez & Fortney 2014, Aguichine et al. 2021,
  Turbet et al. 2020 for steam worlds); never type in a power-law
  approximation (the demo's toy curves say "toy" in their labels).
- Set limits **after** the builders when you call `upper_limit`: it sizes
  its arrow from the current y-range.
