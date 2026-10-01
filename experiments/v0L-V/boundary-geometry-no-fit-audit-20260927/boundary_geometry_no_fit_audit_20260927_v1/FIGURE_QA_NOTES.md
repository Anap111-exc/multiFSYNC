# Figure QA notes

## Export bundle

- Figure 1 and Figure 2 were exported as editable PDF and SVG, plus 600 dpi PNG and LZW-compressed TIFF.
- Final width is 183 mm. Figure heights are 145 mm and 105 mm, respectively.
- PNG and TIFF dimensions are 4322 x 3425 pixels for Figure 1 and 4322 x 2480 pixels for Figure 2; embedded raster resolution is 600 dpi.
- The SVG files retain text objects (`<text>` elements), and `pdftotext` can recover the PDF labels.

## Automated checks

- The rendering script parsed successfully with `Rscript parse()` and completed without an export error.
- The static figure preflight returned 18 PASS, 2 WARN and 0 FAIL.
- One warning is the validator's generic request to run an R parse check; that parse check passed.
- The other warning says final width was not statically detected. The renderer supplies the width through the `width_mm` function argument; both calls use 183 mm.
- The PDF content-stream font-size audit reports 1 pt. This is a known false positive for the R Cairo device in this figure: Cairo writes a 1 pt text state under a scaling transformation. The renderer's smallest configured text is 7 pt, the SVG retains editable text, and visual inspection at final size found no undersized or unreadable labels.

## Visual inspection

- No clipping, panel overlap, truncated labels or legend ambiguity was found.
- The representative functions remain fully visible, including all endpoint regions.
- Grey shading identifies the two audited 10% boundary intervals.
- The extreme non-periodic curves approach zero smoothly at both ends; no seam spike is visible.
- The periodic curves close smoothly at the seam; high boundary-energy examples show a broad boundary level or peak rather than a discontinuity.

## Reproduction boundary

`render_boundary_geometry_figures_20260927_v1.R` reads the completed CSV tables only. It does not regenerate functions, fit a model, call CAVI or ELBO, or read sealed scientific outcomes.
