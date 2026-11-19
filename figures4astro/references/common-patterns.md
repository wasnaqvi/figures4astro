# Common patterns

Layout techniques that recur across the exoplanet figures studied, with code.

## Flux over residual, no gap

```python
fig, axes = xp.flux_residual_axes(ncols=2, size=xp.figsize("aas", "double", 0.38))
```
`hspace=0`, x shared, and the flux panel's lowest and the residual panel's
highest y tick label pruned, so "0.992" and "2000" do not print on top of
each other at the seam. Height ratio 3:1 (use 2:1 when the residuals are the
point).

## Raw under binned, without raw setting the limits

`raw_and_binned` adds the raw points with `ax.add_artist`, which leaves the
data limits alone. Autoscaling then follows the binned points and models. If
you plot raw points yourself with `ax.plot`, set the limits by hand: a few
5σ outliers will otherwise squash the signal into the middle 20% of the panel.

## Vertical offset stacks

Light curves by wavelength, spectra by visit, spectra by planet. Pick a
constant offset ≥ (signal amplitude + 4σ + label height) and keep it
constant; say "offset for clarity" in the caption. Label each row at its
baseline in the row's colour. Order:
- wavelength stacks: blue at the top
- population stacks: coolest at the top
- multi-visit spectra: chronological from the top

## Shared colorbar for a grid

```python
fig, axes = plt.subplots(2, 3, figsize=xp.figsize("aas", "double", 0.5), constrained_layout=True)
ims = [xp.plot_map(a, lon, lat, f, vmin=vmin, vmax=vmax) for a, f in zip(axes.flat, fields)]
xp.colorbar(fig, ims[0], list(axes.flat), label=xp.label("Temperature", "K"), location="bottom", aspect=40)
```
Same `vmin`/`vmax` on every panel. A grid with different colour scales per
panel invites a false comparison.

## Inset zoom

```python
ins = ax.inset_axes([0.55, 0.08, 0.42, 0.42])
xp.plot_mass_radius(ins, M, R, Merr, Rerr, c=Teq, log=False)
ins.set_xlim(4, 11); ins.set_ylim(1.9, 2.8); ins.set_xlabel(""); ins.set_ylabel("")
ax.indicate_inset_zoom(ins, edgecolor="k", lw=0.6)
```

## Broken time axis (RVs over several seasons)

```python
fig, (a, b) = plt.subplots(1, 2, sharey=True, figsize=xp.figsize("aas", "double", 0.3),
                           gridspec_kw=dict(wspace=0.04, width_ratios=(1, 1)))
a.set_xlim(1640, 1830); b.set_xlim(3240, 3275)
a.spines.right.set_visible(False); b.spines.left.set_visible(False)
a.tick_params(right=False, which="both"); b.tick_params(left=False, labelleft=False, which="both")
d = dict(marker=[(-1, -1.5), (1, 1.5)], ms=6, ls="none", color="k", mec="k", mew=0.8, clip_on=False)
a.plot([1, 1], [0, 1], transform=a.transAxes, **d); b.plot([0, 0], [0, 1], transform=b.transAxes, **d)
fig.supxlabel("Time (BJD − 2457000)")
```
Width ratios proportional to each segment's time span keep the time scale
equal.

## Molecule bands and labels

`xp.shade_bands(ax, {"H2O": [(1.3, 1.5), (1.8, 2.05), (2.6, 2.9)], "CH4": [(3.2, 3.45)],
"CO2": [(4.2, 4.45)], "CO": [(4.5, 4.9)]})`, alpha about 0.12-0.15 so error bars
stay readable on top. Alternatively, coloured bars along the bottom of the axes
with the species name in its colour (2510.24708 Fig 3).

## Legends that do not cover data

1. Empty corner inside the axes, `frameon=False`.
2. Above the axes: `ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4)`.
3. Grouped legends (Models / Instruments / Opacities) as separate
   `ax.legend` calls with titles, stacked outside the right edge
   (2605.01589). Keep the first one with `ax.add_artist(leg1)`.
4. For many lines, label them directly at their right ends in their own
   colour and drop the legend.

## Limits and detections

- Upper limits: `xp.upper_limit(ax, x, ylim95)` after the limits are set.
- Posterior histograms: 95% arrow from `plot_posterior(..., limit="upper")`.
- Non-detection bands: shade the excluded region pale grey with a short
  label ("excluded at 3σ") only if the caption cannot carry it.

## Year markers on BJD axes

```python
top = ax.secondary_xaxis("top")
top.set_ticks([bjd(2023, 1, 1) - 2450000, bjd(2024, 1, 1) - 2450000], labels=["2023", "2024"])
```
Or short vertical ticks with text inside the frame (2507.21751 Fig 1).

## Panel tags outside the frame

`xp.panel_label(ax, "a", inside=False)` puts the tag above the top-left
corner. Use it when data fill the panel (dense stacks, maps).
