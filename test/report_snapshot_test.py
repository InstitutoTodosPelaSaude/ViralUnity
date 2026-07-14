"""HTML-content snapshot + JS-source guards for the consensus report.

Regenerates the report from the committed fixtures (an Illumina paired-end run
with a known x2 mapping case, a Nanopore single-end run with a known
zero-coverage run at positions 1-5, and a segmented Influenza-like run) and
asserts the classes of defect the UX review fixed cannot silently return:
grouped counts, a mapping-rate column, honest zeros (no clamp), linear-default
depth, both progressive-disclosure accordions, and — critically — that the scale
toggle rebuilds via Plotly.react instead of a yaxis.type relayout (the ~6 s freeze
a static snapshot alone would miss). Metadata is taken from each fixture's config
YAML, exercising that path end to end.
"""

import json
import os
import re
import tempfile
import unittest

from viralunity.scripts.python.generate_consensus_report import (
    build_report_metadata,
    load_run_config,
    render_report,
    write_report,
)

FIXTURE_ROOT = os.path.join(os.path.dirname(__file__), "fixtures", "report")


def _render(fixture):
    d = os.path.join(FIXTURE_ROOT, fixture)
    cfg = os.path.join(d, "config.yml")
    metadata = build_report_metadata(load_run_config(cfg) if os.path.isfile(cfg) else None, d)
    return render_report(d, metadata)


def _coverage_json(html):
    m = re.search(r'id="coverageData"[^>]*>(.*?)</script>', html, re.S)
    return json.loads(m.group(1))


def _app_js(html):
    """The report's own <script> tail (the last one — the huge first block is the
    vendored plotly.js, which contains unrelated tokens like ``Math.max(1,``)."""
    return re.findall(r"<script>(.*?)</script>", html, re.S)[-1]


class TestReportSnapshot(unittest.TestCase):
    def _assert_common_invariants(self, html):
        app = _app_js(html)
        # honest zeros: the log clamp is gone; the line breaks at true zeros.
        self.assertNotIn("Math.max(1,", app)
        self.assertIn("connectgaps: false", app)
        # both progressive-disclosure accordions are built.
        self.assertIn("By sample", app)
        self.assertIn("By segment", app)
        # the scale toggle rebuilds via react, NOT a yaxis.type relayout (guards
        # the ~6 s freeze). yaxis.type appears only unquoted in a comment; a real
        # relayout would quote it as an object key.
        self.assertIn("Plotly.react", app)
        self.assertIn("Log10", app)
        self.assertNotIn("'yaxis.type'", app)
        self.assertNotIn('"yaxis.type"', app)
        # the mapped column is a rate, not a raw count (server-rendered header).
        self.assertIn(">Mapped %<", html)
        # depth is linear by default: the old baked "(log)" axis title is gone.
        self.assertNotIn("Depth (log)", html)

    def test_illumina_paired_snapshot(self):
        html = _render("unsegmented")
        self._assert_common_invariants(html)
        # declared metadata surfaces platform + layout + primer scheme.
        self.assertIn("paired-end", html)
        self.assertIn("schemes/sarscov2.primers.bed", html)
        # 6-7 digit counts are grouped.
        self.assertIn("10,000", html)
        # paired x2 reconciliation: sample-A mapped 9000 / (2 x 9500) = 47.4%.
        self.assertIn("47.4%", html)

    def test_nanopore_single_end_snapshot(self):
        html = _render("nanopore")
        self._assert_common_invariants(html)
        # single-end run: no doubling, and the report says QC was not run.
        self.assertIn("single-end", html)
        self.assertIn("no read-QC step", html)
        # 256410 / 462000 = 55.5% (undoubled denominator).
        self.assertIn("55.5%", html)
        self.assertIn("462,000", html)
        # the known zero-coverage run is preserved verbatim in the embedded data
        # (not clamped away): barcode05 begins with five depth-0 positions.
        cov = _coverage_json(html)
        first_five = cov["barcode05"][0]["y"][:5]
        self.assertEqual(first_five, [0.0, 0.0, 0.0, 0.0, 0.0])

    def test_segmented_snapshot(self):
        html = _render("segmented")
        self._assert_common_invariants(html)
        # a segment selector is offered for the by-segment view.
        self.assertIn("segmentSelect", html)

    def test_writes_single_self_contained_file(self):
        with tempfile.TemporaryDirectory() as d:
            dest = os.path.join(d, "report.html")
            fixture = os.path.join(FIXTURE_ROOT, "nanopore")
            metadata = build_report_metadata(
                load_run_config(os.path.join(fixture, "config.yml")), fixture
            )
            write_report(fixture, dest, metadata)
            self.assertTrue(os.path.isfile(dest))


if __name__ == "__main__":
    unittest.main()
