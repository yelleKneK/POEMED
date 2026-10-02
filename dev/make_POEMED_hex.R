# Canonical generator for the POEMED hex logo (man/figures/logo.svg).
# Azurite stone ground, brass trim. The picture is the method's story: an
# exposure X (a bold capital in square brackets) reaching an outcome Y (italic)
# through many candidate mediators (a tall fan with a vertical ellipsis for
# "many more"), a few of them active (brass). The X->M side fans OUT to the
# boxes; the M->Y side mirrors it in reverse, converging to a single point at Y
# with one arrowhead (same size as the X->M heads). Tagline names the method;
# wordmark POEMED = POwer-Enhanced MEDiation.
# Run from the package root:
#   Rscript dev/make_POEMED_hex.R          # writes man/figures/logo.svg
#   Rscript dev/make_POEMED_hex.R --png    # also renders logo.png (1200 x 1260)
#                                          # and logo-hires.png (4000 x 4200)
# The PNG step uses the family renderer, ../KenKelleyHexFamily/_render.R
# (headless Chrome at the true aspect; qlmanage forces a square thumbnail).

GROUND <- "#1c5a83"   # azurite
METAL  <- "#c69c3f"   # brass (copper is DMAR's trim; POEMED uses brass)
INK    <- "#e3e9ee"; WHITE <- "#ffffff"

stealth <- function(x, y, ang, col, len = 2.4, w = 1.4) {
  P <- list(c(0, 0), c(-len, w), c(-len * 0.55, 0), c(-len, -w))
  rp <- vapply(P, function(p) c(x + p[1] * cos(ang) - p[2] * sin(ang), y + p[1] * sin(ang) + p[2] * cos(ang)), numeric(2))
  sprintf('<polygon points="%s" fill="%s"/>', paste(sprintf("%.2f,%.2f", rp[1, ], rp[2, ]), collapse = " "), col)
}
arrow <- function(sx, sy, bx, by, col, wln, len = 2.4) {
  ang <- atan2(by - sy, bx - sx); ex <- bx - len * cos(ang); ey <- by - len * sin(ang)
  paste0(sprintf('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="%g" stroke-linecap="round" opacity="0.5"/>', sx, sy, ex, ey, col, wln), stealth(bx, by, ang, col, len = len))
}
boxc <- function(cx, cy, w, fill, stroke, sw) sprintf('<rect x="%.2f" y="%.2f" width="%g" height="%g" fill="%s" stroke="%s" stroke-width="%g"/>', cx - w / 2, cy - w / 2, w, w, fill, stroke, sw)
# cy is %.1f, not %g: the shipped art carries cy="96.0" and cy="112.0", and %g
# drops a trailing .0. Byte-for-byte reproduction of the shipped logo needs it.
ellipsis <- function(cx, ys, col) paste0(vapply(ys, function(y) sprintf('<circle cx="%g" cy="%.1f" r="1.05" fill="%s"/>', cx, y, col), ""), collapse = "")

# Exposure as a bold capital X inside hand-drawn square brackets, centered
# vertically on the X (the bracket glyphs from a font sit unevenly, so the
# brackets are drawn as paths to share the X's vertical center exactly).
xnode <- function(cx, cy) {
  bx <- 6.2; bh <- 4.6; tk <- 1.7   # half-width, half-height, serif length
  br <- function(x, dir) sprintf('<path d="M%.1f,%.1f L%.1f,%.1f L%.1f,%.1f L%.1f,%.1f" fill="none" stroke="%s" stroke-width="0.9" stroke-linecap="round"/>',
                                 x + dir * tk, cy - bh, x, cy - bh, x, cy + bh, x + dir * tk, cy + bh, INK)
  paste0(br(cx - bx, 1), br(cx + bx, -1),
         sprintf('<text x="%g" y="%g" text-anchor="middle" fill="%s" font-size="8" font-weight="800">X</text>', cx, cy + 2.7, INK))
}
ynode <- function(cx, cy) sprintf('<text x="%g" y="%g" text-anchor="middle" fill="%s" font-size="8.6" font-weight="400" font-style="italic">Y</text>', cx, cy + 2.9, INK)

fan_body <- function() {
  # Reconciled 2026-08-05 with the shipped logo.svg, which this script had
  # drifted behind. Two changes: the ellipsis is SIX dots in two groups, one
  # above and one below the direct effect (the mediators continue in both
  # directions), not three in a single run; and a direct X -> Y effect is drawn
  # straight across at yc in the same thin line style, which the older version
  # omitted entirely.
  top <- c(64, 71, 78, 85, 92); bot <- c(116, 123, 130, 137, 144); yc <- 104
  ell <- c(96.0, 98.8, 101.6, 106.4, 109.2, 112.0)
  ys <- c(top, bot)
  ina <- paste0(vapply(ys, function(y) arrow(54, yc, 97, y, INK, 0.4), ""), collapse = "")
  # The M->Y side is the reverse of the X->M fan: every candidate mediator sends
  # a faint line that CONVERGES to one point just left of Y, where a single
  # arrowhead -- identical in size to the X->M heads -- enters Y. A converging
  # fan cannot give each mediator its own separated head the way the diverging
  # X->M fan does (equal heads would pile up at one spot), so the shared head
  # keeps the picture clean and balances the single [X] node on the left.
  V <- 145.8                                   # convergence point; head tip at V + 2.4
  out <- paste0(
    paste0(vapply(ys, function(y) sprintf(
      '<line x1="103" y1="%g" x2="%.2f" y2="%g" stroke="%s" stroke-width="0.4" opacity="0.5" stroke-linecap="round"/>',
      y, V, yc, INK), ""), collapse = ""),
    stealth(V + 2.4, yc, 0, INK),
    # the DIRECT effect: X straight to Y at the same height, in the same thin
    # line style, ending at the convergence point so it shares that arrowhead
    sprintf('<line x1="54" y1="%g" x2="%.2f" y2="%g" stroke="%s" stroke-width="0.4" opacity="0.5" stroke-linecap="round"/>',
            yc, V, yc, INK))
  bt <- paste0(vapply(seq_along(top), function(i) boxc(100, top[i], 4.4, if (i %in% c(2, 4)) METAL else "none", INK, 0.7), ""), collapse = "")
  bb <- paste0(vapply(seq_along(bot), function(i) boxc(100, bot[i], 4.4, if (i == 3) METAL else "none", INK, 0.7), ""), collapse = "")
  paste0(ina, out, bt, ellipsis(100, ell, INK), bb, xnode(46, yc), ynode(152, yc))
}
svg <- paste0(
  '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 200 210" font-family="\'Helvetica Neue\',Arial,sans-serif">',
  '<polygon points="100,5 186.6,55 186.6,155 100,205 13.4,155 13.4,55" fill="', GROUND, '" stroke="', METAL, '" stroke-width="9" stroke-linejoin="round"/>',
  # Title Case at the family tagline size and calibrated baselines (TAG_Y = 34.9,
  # second line +7.4). This script had drifted behind the shipped logo.svg and
  # would have reverted the tagline to ALL CAPS at 5.8; reconciled 2026-08-05.
  # font-family is stated explicitly: the wrapper's declaration is stripped when
  # the inner SVG is embedded elsewhere (the family sheet), where it would fall
  # back to serif.
  sprintf('<text x="100" y="34.9" text-anchor="middle" fill="%s" font-family="\'Helvetica Neue\',Arial,sans-serif" font-size="6.6" font-weight="700" letter-spacing="0.1">Power Enhancement</text>', WHITE),
  sprintf('<text x="100" y="42.3" text-anchor="middle" fill="%s" font-family="\'Helvetica Neue\',Arial,sans-serif" font-size="6.6" font-weight="700" letter-spacing="0.1">High-Dimensional Mediation</text>', WHITE),
  fan_body(),
  # Wordmark. POEMED takes the family's six-letter size (24, as SIEVER does)
  # with tight tracking. It is placed so that the gap between the E and the
  # M sits on the hexagon's center line: text-anchor="middle" centers the
  # word's advance width, but M and D are wider than P and O, so the anchor
  # moves right by 2.33 (x = 102.33). The baseline (y = 170.4) is then set
  # so both end letters just meet the border. MEASURED 2026-09-29 on the
  # 1200 x 1260 render, white ink against the inner edge of the border: E-M
  # gap center x = 100.00; the P's lower left clears by +0.55 and the D's
  # lower right by +0.77 (the family target is 0.5 to 1.0, "just kissing");
  # ink x = 49.83 to 155.00; cap top 152.75, 6.2 units below the lowest
  # mediator box. font-family is stated on the element (family rule: a
  # wrapper's declaration is stripped when the inner SVG is embedded
  # elsewhere). If the size, tracking, or position changes, re-render and
  # re-measure.
  sprintf('<text x="102.33" y="170.4" text-anchor="middle" fill="%s" font-family="\'Helvetica Neue\',Arial,sans-serif" font-size="24" font-weight="800" letter-spacing="0.5">POEMED</text>', WHITE),
  '</svg>')

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
writeLines(svg, "man/figures/logo.svg")
cat("wrote man/figures/logo.svg (azurite + brass, POEMED, apex M->Y)\n")

# PNGs only on request, so that an ordinary run rewrites the SVG and nothing
# else (tools/release_checks.R runs this script in a scratch directory and
# compares the SVG byte for byte). Chrome sizes its window to the SVG's
# width and height, so rasterize from a dimensioned copy; the shipped
# logo.svg stays viewBox-only.
if ("--png" %in% commandArgs(trailingOnly = TRUE)) {
  renderer <- file.path("..", "KenKelleyHexFamily", "_render.R")
  if (!file.exists(renderer)) stop("the family renderer is not at ", renderer, call. = FALSE)
  source(renderer)
  tmp <- file.path("man", "figures", ".logo_render.svg")
  writeLines(sub("<svg ", '<svg width="200" height="210" ', svg, fixed = TRUE), tmp)
  render_svg(tmp, file.path("man", "figures", "logo.png"), 200, 210, 6)          # 1200 x 1260
  render_svg(tmp, file.path("man", "figures", "logo-hires.png"), 200, 210, 20)   # 4000 x 4200
  unlink(tmp)
}
