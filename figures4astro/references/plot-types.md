# Recipes by figure type

Each recipe names its model figure (arXiv ID, figure number) so the original
can be pulled up. Builder calls refer to [api.md](api.md);
`scripts/demo_gallery.py` has a runnable version of every recipe marked (demo).

---

## 1. White-light transit / eclipse with residuals (demo 01)

*Models*: 2508.08416 Fig 1 (TRAPPIST-1 d, two visits side by side);
2608.17207 Fig 1 (four visits, with an Hα activity panel underneath);
2510.11397 Fig 1 (three eclipses, raw / detrended / residual rows).

- One column per visit (`flux_residual_axes(ncols=n)`), shared y within a row.
- Rows from top: raw flux with the systematics model (optional), detrended flux
  with the astrophysical model, residuals (ppm). Add a stellar-activity row
  (Hα, FWHM) when activity is part of the story.
- x is `Time from mid-transit (h)` or `(min)`. Absolute BJD only when the
  timing is the result.
- Grey spans for ingress/egress or eclipse duration; dotted line at
  mid-eclipse.
- Binned points at about 5-10 min; raw at alpha 0.2-0.3. Print the residual rms
  in the corner of the residual panel.
- Visit name and instrument as a two-line text in an empty corner. No titles.
- Integrations you masked: draw them in light grey, not deleted (2608.17207).

## 2. Spectroscopic light-curve stack (demo 02)

*Models*: 2512.10876 Fig 1 (a: offset light curves coloured by wavelength
with the band printed above each; b: residuals with per-row ppm scales; c:
residual histograms against a unit Gaussian); 2508.08416 Fig 1 bottom.

- Bluest channel at the top, colour from `CMAPS["wavelength"]`. Binned points
  with thin black edges, black model line.
- Offset ≥ depth + 4σ + room for the band label. Labels sit just above
  each baseline at the left.
- The residual column uses the same colours and offsets, with a dashed zero
  line per row. The histogram column overlays a dashed N(0, 1); a mismatch
  shows red noise or underestimated errors at a glance.
- Show 8-16 channels in the main text; all of them go in the appendix.

## 3. Transmission spectrum (demo 03)

*Models*: 2512.10876 Fig 2 (a: data + retrieved median + 1σ/2σ; b: the same
data over filled per-species contributions; c: comparison with an alternative
model); 2508.08416 Fig 4 (flat data against single-species and 1× solar
models); 2606.19487 Fig 2 (limb spectra with filled species regions).

- `plot_spectrum` on a log wavelength axis. Data: open circles with x error
  bars as the bin half-widths.
- Retrieved model as a band (median plus 1σ plus 2σ). Alternative or forward
  models as thin lines in species colours (a single-species model takes that
  species' colour).
- A **contribution panel** below (shared x, no gap): continuum/haze in grey,
  each species filled between the continuum and its model in its fixed
  colour. This shows the reader which bump is which molecule.
- If space allows only one panel, use `shade_bands` (pale spans, species
  names along the top) instead of the contribution panel.
- Several visits or reductions: offset vertically with a label per row
  (2608.17207 Fig 2a), or overlay with distinct markers. Never overlay two
  binnings of the same data.
- Flat-line and 1× solar clear-atmosphere references are what make a "flat"
  spectrum informative. Include them when the result is a non-detection.
- A near-flat spectrum on an absolute depth axis can look featureless. Also
  plot it in scale heights (`(depth − mean) / (2 H R_p / R⋆²)`) when comparing
  planets.

## 4. Emission spectrum / eclipse depths

*Models*: 2510.11397 Fig 2 (MIRI eclipse depth against bare-rock and
atmosphere models, with a side panel comparing the reductions);
2509.02120 Fig 3 (TRAPPIST-1 b: airless vs atmosphere models with data
crosses); 2605.04936 Fig 2 (day/night/east/west spectra with blackbody fits,
and deviation from blackbody below).

- Fp/F⋆ (ppm) vs wavelength. Bare-rock (A_B = 0) and blackbody references
  dashed black or grey; atmosphere models in colour.
- Photometric points: horizontal bar for the bandpass, open marker.
- Comparing reductions of the same eclipse: a narrow side panel with one
  column per pipeline (2510.11397) rather than overplotting.
- Phase-resolved emission: one colour per orbital phase from
  `CMAPS["phase"]` with a colorbar in degrees (2607.26182 Fig 1), or the
  `DAYSIDE`/`NIGHTSIDE` pair for just two.

## 5. Phase curve (demo 04)

*Models*: 2509.09760 Fig 5 (WASP-121 b: data, total model, separate
reflection and emission components with bands, residuals); 2608.21519
Figs 1-2 (TOI-561 b NRS1/NRS2 phase curves, brightness maps, Bond albedo and
ε histograms); 2605.04936 Fig 1 (phase curve above a wavelength × phase
flux-ratio image).

- Fp/F⋆ (ppm) vs orbital phase 0-1 with eclipse at 0.5, or −0.5 to 0.5 with
  transit at 0. Pick one per paper.
- Transits masked or clipped, not shown at full depth (it flattens the curve).
- Component curves: emission `DAYSIDE`, reflection `REFLECTED`, each with a
  1σ band. Total model black.
- Spectroscopic phase curves: an `imshow` of Fp/F⋆ over (phase, λ) with
  `inferno` under the white-light curve (2605.04936), or a grid of small
  panels per channel.
- Derived quantities (offset, amplitude, A_B, ε) vs wavelength as open
  circles, with the bands shaded where molecules absorb (2605.01589 Fig 6).

## 6. T-P and abundance profiles (demo 05)

*Models*: 2510.24708 Fig 4 (retrieved T-P with envelopes against GCM and
other retrievals; log H₂O vs P beside it); 2605.01589 Fig 4 (day/night T-P
for many GCM variants, with contribution-function bars beside the axis);
2401.06608 Figs 1-6 (Leconte: P on log y, several narrow panels sharing it).

- `plot_tp`: P (bar) log, increasing downward; median plus 1σ/2σ shaded in x.
- GCM or forward profiles as dashed lines (dayside orange, nightside blue).
- **Contribution bars** (`contribution_bar`) at the left edge, one per band,
  in band colours: they show where the spectrum is sensitive.
- Add abundance panels sharing y (`sharey=True`, `wspace` about 0.06),
  x = log₁₀ VMR, species colours.
- Condensation curves dotted, labelled at their end.

## 7. Retrieval posteriors (demo 06)

*Models*: 2512.10876 Fig 3 (step histograms per model variant, colour by
retrieval setup, median ± 1σ bars and arrows for limits above);
2508.08416 Fig 4 lower (posterior grid with the species name in its colour
in the corner); 2605.04936 Fig 3 (posteriors on a pressure axis).

- 1-D marginals: `plot_posterior` (filled step histogram, no y ticks). Upper
  limits get a 95% arrow; quote the number in the text.
- Comparing retrievals: same panel, one colour per setup, unfilled steps
  except the preferred one. Put the legend in the first panel only.
- Joint 2-D constraints (e.g. P_surf vs metallicity, 2508.08416 Fig 3):
  filled density plus 1/2/3σ contour lines labelled inline, regime boundaries
  dashed with boxed text.
- Corner plots: the `corner` package, `corner.corner(samples, labels=...,
  color="k", plot_datapoints=False, fill_contours=True, levels=(0.39, 0.86),
  smooth=1)` (1σ and 2σ in 2-D), at most about 8 parameters in the main text.
  Overlay a second posterior with a second colour on the same figure
  (2606.19487 Fig 5).

## 8. Radial velocities (demo 07)

*Models*: 2507.21751 Figs 1-6 (Proxima, NIRPS: stacked ΔRV / ΔFWHM / ΔCRX /
BERV time series, GP fit with periodograms beside each row, phase-folded b
and d with large open phase bins and O − C); 2604.07447 Figs 2-5 (multi-
instrument phase folds, periodogram sequence "raw → activity removed →
planet removed"); 2609.03006 Figs 2-6 (TRAPPIST-1 SPIRou + NIRPS, broken
time axes, K posteriors against predictions).

- **Time series**: `BJD − 2450000 (d)`; one marker shape and colour per
  instrument (`INSTRUMENT_MARKERS/COLORS`); GP mean as a line with a 1σ band;
  season gaps with broken axes (see common-patterns.md) or year ticks along
  the top.
- **Activity indicators** stacked below the RVs with shared x (FWHM, BIS,
  CRX, dLW, Hα), each with its own periodogram to the right.
- **Phase fold**: `plot_rv_phased`. Thin grey posterior draws, black
  Keplerian, large open navy phase bins, O − C panel. One panel per planet,
  with other planets and the GP subtracted (say so in the caption).
- **Periodograms**: log period, FAP 10/1/0.1% as coloured dashed/dotted
  lines, vertical lines at P_b, P_rot, 1 yr / 1 d aliases. Show the sequence
  of periodograms as signals are removed, top to bottom.

## 9. Mass-radius and population diagrams (demo 08)

*Models*: 2609.04318 Fig 1 (Helled group: log-log M-R, archive coloured by
instellation, inset zoom on the targets, H/He and Fe curves labelled along the
curve); 2604.07447 Fig 6 (NIRPS super-Earths: Earth-like plus water-fraction
curves at three temperatures, points coloured by T_eq); 2602.11271 Fig 3
(density vs mass with composition lines, targets as stars); 2509.02120
Fig 1 (cosmic shoreline: XUV vs escape velocity with a blurred band).

- Background population small and light grey (or coloured by a third
  variable with a colorbar). The paper's planets have black edges, error bars,
  and targets as gold stars with names.
- Composition curves from published tables, each in a distinct colour **and**
  dash. Label them in a legend or along the curve (rotated text), not both.
- Inset zoom (`ax.inset_axes` plus `ax.indicate_inset_zoom`) when the targets
  crowd in a small region (2609.04318).
- Regime boundaries (the radius valley, the cosmic shoreline) as shaded bands
  with low-alpha gradient edges, not hard lines.
- Trend plots (2608.15717 Fig 2): log axes, colour by host type, statistic
  (Kendall τ, p) in a small boxed corner label.

## 10. Interior structure

*Models*: 2609.04318 Figs 2-4 (wedge "pie-slice" interiors with layer
names, a table of r/R, P, T at each boundary, and a mass-fraction bar
beneath; alternative solutions in a grey side box; T-P adiabats with σ
bands).

- Draw wedges with `matplotlib.patches.Wedge` (centre at the bottom, opening
  about 40°), radii from the model, layer colours from `LAYER_COLORS`, layer
  names centred in the layer. Composition gradients as hatched blends
  between two layer colours.
- Boundary values as a small aligned table to the right, with leader lines.
- Mass fractions as a horizontal stacked bar under the wedge, percentages
  inside; a tiny layer's value goes in a callout box beneath it.

## 11. Maps: GCM fields and eclipse maps (demo 09)

*Models*: 2608.17207 Fig 6 right (cloud column maps, lon −180..180,
terminators dashed); 2510.24708 Figs 1, 6 (eclipse maps in Mollweide grids
per wavelength, `inferno`); 2606.19487 Fig 2 (polar view of the limb
temperature at t2, t3); 2401.06608 Fig 1 (cloud-resolving x-z slices: `RdBu_r`
vertical velocity, `inferno` humidity).

- Equirectangular `plot_map` for GCM output: ticks every 90° and 45°, dashed
  terminators at ±90°, substellar point marked, contours labelled inline.
- Eclipse / brightness maps: `projection="mollweide"` or `"hammer"`,
  graticule thin grey, only the visible hemisphere coloured, colorbar
  horizontal underneath.
- Many wavelengths: a grid of small maps sharing one colorbar.
- Signed fields centred on zero with `TwoSlopeNorm(vcenter=0)` and `RdBu_r`.

## 12. Population spectral stack, Sing style (demo 10)

*Models*: Sing et al. 2016 Nature Fig 1 (rules in SKILL.md, House rules).

- `offset_stack` with rows sorted coolest first (the first row is drawn at the
  top). Colour by T_eq rank on `CMAPS["teq_stack"]`.
- y = relative altitude in scale heights, constant offset (widen a gap only
  where rows would collide); name and T_eq in a right-hand gutter inside the
  frame in the row's colour.
- Points only, one resolution per figure (R≈25 per detector), detector gap
  kept, no detector offset, no group labels, no in-figure caption.
- Tall single-column (aspect 1.2-1.6) or full-page for more than about 12 rows.

## 13. Discovery-paper staples

*Models*: 2603.18233, 2603.10905 (TOI discoveries: TESS TPF with the
aperture and Gaia neighbours, multi-sector light curves, phase folds,
high-resolution imaging contrast curves, SED, RV fold).

- TPF: `imshow` of the pixel flux (`viridis`), the aperture outlined in red,
  Gaia sources as red dots scaled by Δmag, a legend for Δm, and a compass.
- Contrast curve: Δmag vs separation (arcsec), with the speckle image as an
  inset.
- SED: log λF_λ vs log λ, the model in black, photometry by survey in distinct
  marker shapes and colours, horizontal bars for filter widths.
- Folded photometry from several instruments: offset rows, one colour per
  instrument, models in black (2602.11271 Fig 1).

## 14. Flow charts (pipelines, analysis overviews)

*Model*: `examples/showcase.py` `fig_flowchart` (JWST time-series pipeline).

- `flow_canvas` + `flow_box` + `flow_arrow`, full page width, aspect about 0.5.
- Colour encodes the **kind of step**, not the order: data (grey), processing
  (blue), fitting / inference (orange), science products (green), external
  inputs (white, dashed). Show the legend once, above the boxes.
- Each box: a bold title of 1-3 words, then at most two short detail lines
  (about 22 characters each at 8 pt single-width boxes).
- Main flow solid, left to right; wrap to a second row right to left (a
  "snake") rather than shrinking boxes. Optional or external inputs dashed,
  entering from below. Feedback loops curved and dashed, in the colour of the
  step that triggers them, with a two-word italic label.
- Keep box sizes identical within a row; align centres on a grid.
- No drop shadows, gradients or icons. A flow chart in a paper is a figure,
  so it follows the same font sizes as the data figures.

