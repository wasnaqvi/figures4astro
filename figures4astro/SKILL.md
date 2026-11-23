---
name: figures4astro
description: >-
  Publication-ready matplotlib figures for exoplanet science papers, talks and
  proposals: transit / eclipse / phase-curve light curves with residuals,
  stacked spectroscopic light curves, transmission and emission spectra with
  retrieved models and molecule bands, T-P profiles, retrieval posteriors,
  radial velocities and periodograms, mass-radius and population diagrams,
  interior structure, GCM / eclipse maps, Sing-style population spectral
  stacks, and pipeline / analysis flow charts. House style distilled from recent figures by the Cowan, Dang,
  Artigau, Helled, Queloz, Bouchy, Komacek and Leconte groups, with a tested
  helper module (scripts/exoplot.py). Use whenever making or polishing such a
  figure in matplotlib. Not for interactive dashboards (Plotly, Bokeh),
  throwaway exploratory plots, or Illustrator-first schematics.
---

# Exoplanet figure making

Open `references/` only as needed. Start with the table, then the recipe for
the figure type in [references/plot-types.md](references/plot-types.md).

## Workflow

1. **Pin the target**: journal (AAS / A&A / MNRAS / Nature), single or double
   column, and what the figure must show in one sentence. If any of this would
   change the layout and is not stated, ask.
2. **Copy `scripts/exoplot.py` into the project** (next to the figure script,
   e.g. `figures/exoplot.py`). Do not import it from `~/.claude`: a paper's
   figures must rebuild from the repository alone.
3. `xp.use_style()` once, then build at **printed size** with
   `xp.figsize(journal, "single"|"double", aspect)`. Never draw big and shrink.
   **Text must be big**: `xp.use_style(fontsize=10.5-11)`, so axis labels and
   tick labels are 10-11 pt and legends, object names and colour-bar labels
   9.5-10.5 pt on the page; nothing below 9 pt (7.5 pt only for a number set
   inside a marker). The exoplot default of 8 pt reads as too small in
   Overleaf.
4. Follow the recipe for the figure type; reuse the fixed colour maps
   (`MOLECULE_COLORS`, `INSTRUMENT_COLORS`, `DAYSIDE`/`NIGHTSIDE`, `CMAPS`).
5. `xp.finalize(fig, "figures/<stem>")` writes **PDF + PNG + SVG**.
6. **Look at the PNG before calling it done** (Read the image). Check: every
   label and legend entry readable at printed size, nothing
   clipped, no tick labels colliding where panels meet, legends off the data,
   raw scatter not setting the limits, every axis has a unit.

`python scripts/demo_gallery.py OUT_DIR` renders one synthetic example of
every figure type — run it to see the style or to check a change to exoplot.

## The house style in eight rules

1. **Closed box, ticks inward on all four sides, minor ticks on.** This is the
   one thing every group studied shares (Cowan/Piaulet-Ghorayeb, NIRPS,
   Helled, Leconte, Queloz). Not the top/right-spineless "data-science" look.
2. **Data are open circles** (white face, black edge, black error bars, no
   caps). Unbinned data go underneath as faint rasterized dots and never set
   the axis limits.
3. **Every fit gets a residual panel** sharing x, no gap, in ppm (light curves)
   or m s⁻¹ (RVs, labelled O − C), with a dashed zero line.
4. **Models are lines; uncertainty is shading of the same hue**, 1σ darker and
   2σ lighter. Never draw a model as points or data as a line.
5. **Colours carry meaning and keep it across the paper**: H₂O blue, CO₂ red,
   CH₄ purple, CO green, SO₂ orange; dayside orange / nightside blue; NRS1 blue
   / NRS2 red; NIRPS dark red / HARPS orange. Spectroscopic channels blue →
   red with wavelength.
6. **Wavelength axes are log-scaled with plain-number ticks** (0.6, 0.8, 1,
   1.5, 2, 3, 4, 5) when the range exceeds about a factor of three.
7. **Panel tags** are bold lowercase `a b c` (Nature, AAS) or `(a)` (A&A,
   MNRAS); one convention per paper.
8. **Units on every axis** in the journal's brackets: `Wavelength (μm)` for
   AAS/MNRAS, `Wavelength [μm]` for A&A (`xp.label()` does this).

## House rules

- **Labels and legends must be big and visible.** Axis labels, tick labels,
  legend text, colour-bar labels and object names are all at least 9.5 pt at
  printed size (10-11 pt for axes), and markers scale up with them (legend
  handles at the size of the plotted markers). If a figure goes in at less
  than full `\textwidth`, scale the fonts up by the same factor. When the
  bigger text crowds, make the figure taller, shorten or wrap labels
  (`"NIRSpec\nBOTS"`) or retune label offsets; never shrink the text back.
- Deliver **PNG + PDF + SVG**, always.
- **Points are never joined by lines**; a line means a model.
- **Never overlay two binnings** (e.g. R=100 and R=25) in one panel: make
  one figure per resolution.
- **No captions, notes or annotations inside the figure** beyond panel tags,
  legends and object names. Explanations belong in the LaTeX caption.
- **Population spectra stacks** (Sing et al. 2016 style): coolest at the TOP,
  coloured by T_eq **rank** on a blue → orange ramp (`CMAPS["teq_stack"]`),
  constant offsets, y label `Relative altitude in scale heights, z(λ)/H_eq`,
  name + T_eq only beside each row, thin separators between temperature
  groups but no group labels, R≈25 per detector with the NRS1/NRS2 gap kept.
- Features should be discernible by eye. If a spectrum is genuinely flat, add a
  reference model (e.g. clear solar H₂/He) rather than display tricks.
- Show the rendered figure before declaring it finished.

## Related files

| File | Open when |
|------|-----------|
| [references/plot-types.md](references/plot-types.md) | Recipe for a specific figure (light curve, spectrum, phase curve, T-P, posteriors, RV, M-R, interior, map, stack, flow chart) |
| [references/api.md](references/api.md) | `exoplot` function signatures and constants |
| [references/design-theory.md](references/design-theory.md) | Why the style is what it is: geometry, fonts, colour, export |
| [references/common-patterns.md](references/common-patterns.md) | Layout tricks: offset stacks, insets, broken axes, limits, shared colorbars |
| [references/demos.md](references/demos.md) | The 24 papers studied, figure by figure, with arXiv IDs |

Adapted from ChenLiu-1996/figures4papers (CC BY-NC 4.0): the generic
bar/heatmap house style was replaced with exoplanet-specific conventions,
helpers and recipes.
