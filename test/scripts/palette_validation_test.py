"""Palette validation for the consensus report's categorical hues.

Guards Finding 8 (RECOMMENDATIONS.md): the light and dark palettes must be
distinct (a dark-specific derivation, not the light hexes reused), keep enough
colour-blind separation between adjacent slots to be told apart, and — in dark
mode, where the old reused palette failed — clear the 3:1 surface-contrast floor.

The colour maths (WCAG relative luminance/contrast, CIE-Lab CIE76 deltaE, and a
Machado-2009 CVD simulation) is implemented here in pure Python so the check runs
in CI with no extra dependency. Thresholds are set below the current palettes'
measured values so the test fails only on a real regression, not on the exact
number matching any particular tool.
"""

import math
import unittest

from viralunity.scripts.python.generate_consensus_report import (
    MAX_CATEGORICAL,
    PALETTE_DARK,
    PALETTE_LIGHT,
)

# Card surfaces the hues are drawn against (from the template's CSS variables).
LIGHT_SURFACE = "#fcfcfb"
DARK_SURFACE = "#1a1a19"

# Machado 2009 severity-1.0 matrices (deuteranopia / protanopia approximations).
DEUTERANOPIA = [
    [0.367322, 0.860646, -0.227968],
    [0.280085, 0.672501, 0.047413],
    [-0.011820, 0.042940, 0.968881],
]
PROTANOPIA = [
    [0.152286, 1.052583, -0.204868],
    [0.114503, 0.786281, 0.099216],
    [-0.003882, -0.048116, 1.051998],
]

# Floors: the current palettes measure above these (worst ADJACENT CVD deltaE
# ~27.6 light / ~13.6 dark; worst ALL-PAIRS CVD deltaE ~11.5 light / ~4.6 dark;
# min contrast ~2.1 light-on-light / ~3.5 dark-on-dark). Each floor sits below the
# measured value so the test trips only on a real regression, not on an exact
# number. Adjacency is the primary guarantee (the palette is drawn in fixed
# order); the all-pairs and light-on-light floors are weaker regression guards —
# tightening them would mean re-designing the palette (a dataviz decision).
MIN_ADJACENT_CVD_DELTA_E = 10.0
MIN_ALLPAIRS_CVD_DELTA_E = 4.0
MIN_DARK_CONTRAST = 3.0
MIN_LIGHT_CONTRAST = 2.0


def _hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i : i + 2], 16) / 255 for i in (0, 2, 4))


def _linear(c):
    return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4


def _luminance(rgb):
    r, g, b = (_linear(c) for c in rgb)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _contrast(a, b):
    la, lb = _luminance(a), _luminance(b)
    return (max(la, lb) + 0.05) / (min(la, lb) + 0.05)


def _lab(rgb):
    r, g, b = (_linear(c) for c in rgb)
    x = r * 0.4124 + g * 0.3576 + b * 0.1805
    y = r * 0.2126 + g * 0.7152 + b * 0.0722
    z = r * 0.0193 + g * 0.1192 + b * 0.9505
    xn, yn, zn = 0.95047, 1.0, 1.08883

    def f(t):
        return t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116

    fx, fy, fz = f(x / xn), f(y / yn), f(z / zn)
    return (116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz))


def _delta_e76(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(_lab(a), _lab(b))))


def _simulate_cvd(rgb, matrix):
    out = [sum(matrix[i][j] * rgb[j] for j in range(3)) for i in range(3)]
    return tuple(min(1.0, max(0.0, c)) for c in out)


def _min_adjacent_cvd_delta_e(palette, matrix):
    sim = [_simulate_cvd(_hex_to_rgb(h), matrix) for h in palette]
    return min(_delta_e76(sim[i], sim[i + 1]) for i in range(len(sim) - 1))


def _min_allpairs_cvd_delta_e(palette, matrix):
    sim = [_simulate_cvd(_hex_to_rgb(h), matrix) for h in palette]
    return min(_delta_e76(sim[i], sim[j]) for i in range(len(sim)) for j in range(i + 1, len(sim)))


class TestPaletteValidation(unittest.TestCase):
    def test_palettes_are_well_formed_and_distinct(self):
        for pal in (PALETTE_LIGHT, PALETTE_DARK):
            self.assertEqual(len(pal), MAX_CATEGORICAL)
            for hexval in pal:
                self.assertRegex(hexval, r"^#[0-9a-fA-F]{6}$")
        # Dark mode uses a dark-specific derivation, not the light hexes reused.
        self.assertNotEqual(PALETTE_LIGHT, PALETTE_DARK)

    def test_adjacent_hues_stay_separable_under_colour_blindness(self):
        for name, pal in (("light", PALETTE_LIGHT), ("dark", PALETTE_DARK)):
            for cvd, matrix in (("deuteranopia", DEUTERANOPIA), ("protanopia", PROTANOPIA)):
                worst = _min_adjacent_cvd_delta_e(pal, matrix)
                self.assertGreaterEqual(
                    worst,
                    MIN_ADJACENT_CVD_DELTA_E,
                    f"{name} palette: worst adjacent {cvd} deltaE {worst:.1f} "
                    f"< floor {MIN_ADJACENT_CVD_DELTA_E}",
                )

    def test_non_adjacent_hues_do_not_collide_under_colour_blindness(self):
        # Adjacency is the main guarantee, but any two hues in the palette can end
        # up compared (a run with several samples), so guard the worst all-pairs
        # separation from regressing below its (weaker) current floor.
        for name, pal in (("light", PALETTE_LIGHT), ("dark", PALETTE_DARK)):
            for cvd, matrix in (("deuteranopia", DEUTERANOPIA), ("protanopia", PROTANOPIA)):
                worst = _min_allpairs_cvd_delta_e(pal, matrix)
                self.assertGreaterEqual(
                    worst,
                    MIN_ALLPAIRS_CVD_DELTA_E,
                    f"{name} palette: worst all-pairs {cvd} deltaE {worst:.1f} "
                    f"< floor {MIN_ALLPAIRS_CVD_DELTA_E}",
                )

    def test_dark_palette_clears_the_contrast_floor(self):
        surface = _hex_to_rgb(DARK_SURFACE)
        for hexval in PALETTE_DARK:
            ratio = _contrast(_hex_to_rgb(hexval), surface)
            self.assertGreaterEqual(
                ratio,
                MIN_DARK_CONTRAST,
                f"dark hue {hexval} contrast {ratio:.2f} < floor {MIN_DARK_CONTRAST}",
            )

    def test_light_palette_clears_the_light_surface_contrast_floor(self):
        # The light hues are drawn on the light card too; guard their (lower)
        # contrast against that surface from regressing.
        surface = _hex_to_rgb(LIGHT_SURFACE)
        for hexval in PALETTE_LIGHT:
            ratio = _contrast(_hex_to_rgb(hexval), surface)
            self.assertGreaterEqual(
                ratio,
                MIN_LIGHT_CONTRAST,
                f"light hue {hexval} contrast {ratio:.2f} < floor {MIN_LIGHT_CONTRAST}",
            )


if __name__ == "__main__":
    unittest.main()
