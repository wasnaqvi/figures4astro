# Design theory

What the exoplanet groups studied actually do, and the reasoning that turns
it into rules. Evidence: about 120 figures from 24 papers, 2024-2026 (see
[demos.md](demos.md)).

## 1. Geometry: design at printed size

The best figures studied (Piaulet-Ghorayeb et al. on TRAPPIST-1 d and the
temperate sub-Neptunes; Suárez Mascareño et al. on Proxima with NIRPS) are
drawn at their final column width, so 8 pt in the script is 8 pt on the page.
The worst have 20 pt axis labels on a single-column figure, or illegible 4 pt
ticks on a full-page grid, because they were drawn large and scaled down.

| journal | single column | full width | units in labels |
|---|---|---|---|
| AAS (ApJ, AJ, PSJ) | 3.50 in | 7.10 in | `(μm)` |
| A&A | 88 mm = 3.46 in | 180 mm = 7.09 in (120 mm also allowed) | `[μm]` |
| MNRAS | 84 mm = 3.31 in | 177 mm = 6.97 in | `(μm)` |
| Nature / Nat. Astron. | 89 mm | 183 mm | `(μm)`, bold `a b c` tags |

Base font 10.5-11 pt for a paper figure; nothing on the page below 9 pt
(axis and tick labels 10-11 pt, legends, object names and colour-bar labels
9.5-10.5 pt). The 7-9 pt range the studied figures use matches the caption
text on paper but reads as too small in Overleaf and on screen, which is how
the figures are actually read and refereed; the exoplot default (8 pt) is
only a floor for dense multi-panel grids. Scale markers and legend handles
up with the text. Aspect ratios that recur:
light curve + residuals 0.9 (single column); two visits side by side
0.35-0.4 (full width); a spectrum 0.6-0.75; a T-P profile taller than wide.

## 2. Frame and ticks

Every group studied uses a **closed box with inward ticks on all four sides
and minor ticks**. It is the convention of observational astronomy and helps
in practice: right and top ticks let the reader read a value off a point far
from the left axis, which matters for 1000-ppm-scale spectra. Spines are thin
(0.6-0.8 pt) so the box frames the data rather than competing with it.

Grids: almost never. Where they appear (the Komacek-group Spitzer catalogue,
the Helled giant-planet mixing figure) they compete with error bars. Use
reference lines instead: zero residual, mid-eclipse, terminators, FAP levels.

## 3. Marks

| element | mark |
|---|---|
| final / binned data | open circle `o`, white face, black edge 0.7 pt, ms 3-4.5, error bars black, **no caps** |
| unbinned data | `.`, ms 1-2, grey or channel colour, alpha 0.15-0.3, **rasterized** |
| data from another work | different marker shape (s, D, ^) and colour, or grey crosses |
| best-fit / median model | solid line, 1-1.2 pt |
| model families | thin lines (0.8-1 pt), differing in colour **and** dash style |
| 1σ / 2σ envelopes | same-hue `fill_between`, alpha about 0.45 / 0.22, no edge |
| limits | bar plus arrow (`xp.upper_limit`); say 2σ or 3σ in the caption |
| reference lines | dashed or dotted, 0.6-0.8 pt, grey |

Rasterize dense scatters (`rasterized=True`) so the PDF stays small and
opens fast, while axes and text stay vector. A 50 000-point light curve as
vector paths is a 20 MB PDF.

## 4. Colour

Colour encodes a variable or a category whose meaning is fixed for the whole
paper.

- **Species** (`MOLECULE_COLORS`): H₂O `#1f77b4` blue, CO₂ red, CH₄ purple,
  CO green, SO₂ orange, NH₃ pink, H₂S olive. These follow the tab10 mapping
  most of the groups converged on (Piaulet-Ghorayeb, Komacek, Leconte).
- **Day / night**: orange `#e66101` / deep blue `#2c3e91`. Reflected light in
  a mid blue, distinct from nightside.
- **Instruments**: NRS1 blue / NRS2 red (Dang-group TOI-561 b); NIRPS dark red
  / HARPS orange (NIRPS GTO papers); NIRISS navy, MIRI brown, TESS blue,
  Spitzer red, HST orange (Komacek WASP-121 b figures).
- **Continuous variables**:
  - temperature maps and brightness: `inferno`
  - planets coloured by T_eq or instellation: `RdYlBu_r`, with black marker
    edges (Helled M-R, NIRPS M-R)
  - spectroscopic channels: `exo_wavelength` (turbo without its near-black
    ends, blue → red with wavelength). A rainbow is acceptable only here, where
    hue order is wavelength order. `Spectral_r` was tried and rejected: its pale
    yellow middle makes the central channels vanish on white.
  - orbital phase: cyclic `twilight_shifted`
  - signed fields (winds, w, residual maps): `RdBu_r`, centred on zero
  - cloud and condensate columns: `Purples` (Piaulet-Ghorayeb cloudy GCM)
  - population stacks: the custom `exo_teq_stack` blue → purple → orange
    ramp, coloured by rank. A diverging map with a pale middle made one planet
    illegible in the demo.
- Avoid `jet`, pure `rainbow` for non-wavelength data, and pastel-on-white
  for data points.

## 5. Typography

Two families are in use:
- **Sans** (DejaVu Sans / Helvetica): Cowan, Dang, NIRPS, Helled, Queloz
  groups. This is the default (`family="sans"`, mathtext `dejavusans`).
- **Serif** matching LaTeX body text (STIX or Computer Modern): Leconte,
  parts of Komacek's work, the JWST Weather Report. Use `family="serif"` or
  `usetex=True` when the paper's other figures are serif, or for A&A if you
  want the figure to blend with the text.

One family per paper. Chemistry uses mathtext subscripts (`xp.mol("CO2")`
gives CO$_2$); symbols are italic (`$R_{\rm p}/R_\star$`,
`$T_{\rm eq}$`), descriptive subscripts upright (`\rm`).

## 6. Labels and units

- Wavelength in μm; transit and eclipse depth in ppm (or % above about
  1%); time as `Time from mid-transit (h)` or `BJD_TDB − 2460000`, with the
  offset in the label and never in an offset text; RV in m s⁻¹; pressure in
  bar on a log axis increasing downward; masses and radii in M⊕/R⊕ below about
  30 M⊕, otherwise M_J/R_J.
- `axes.formatter.useoffset = False`: a "+2.46e6" offset label hides the
  time stamp. Subtract a round number yourself and say so in the label.
- Name objects as the archive does (`TOI-270 d`, with a space before the
  letter).

## 7. Export

`xp.finalize` writes PDF (vector, TrueType fonts embedded via
`pdf.fonttype 42`, so the journal can edit text), PNG at 300 dpi (slides,
GitHub, quick look) and SVG with live text (`svg.fonttype none`, for
Illustrator or Inkscape touch-ups). `bbox_inches="tight"` with a 0.02 in pad.
`fig.align_ylabels()` is applied so stacked panels' y labels line up.

## 8. What not to copy

Seen in the sample and deliberately not adopted:
- Black blobs of 10⁵ unbinned points with binned points lost on top
  (Spitzer catalogue). Faint raw points plus open-circle bins instead.
- Six-instrument legends inside the data area covering points. Put them
  above the axes (`bbox_to_anchor=(0.5, 1.0)`, `ncol=N`) or in an empty corner.
- Corner plots of more than about 8 parameters in a main-text figure. Show the
  4-6 parameters that carry the claim and send the full corner to the appendix.
- Titles above single-panel figures. That job belongs to the caption. Text
  inside a panel names the object, visit or instrument only.
