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
    config = load_run_config(cfg) if os.path.isfile(cfg) else None
    metadata = build_report_metadata(config, d)
    return render_report(d, metadata, config)


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
        # both progressive-disclosure accordions are server-rendered (structure
        # present in the HTML, not reconstructed at load by JS).
        self.assertIn('id="by-sample-details"', html)
        self.assertIn('id="by-segment-details"', html)
        self.assertIn("<summary>", html)
        self.assertIn("By sample", html)
        self.assertIn("By segment", html)
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
        # KPI summary tiles lead the report (server-rendered, five of them).
        self.assertIn('id="kpi-grid"', html)
        self.assertIn("Samples analyzed", html)
        for key in ("samples", "pass_count", "below_warn", "median_coverage", "mean_depth"):
            self.assertIn(f'data-kpi="{key}"', html)
        # threshold labels are derived from the params, not hardcoded 90/70.
        self.assertIn("&ge;90% coverage", html)
        self.assertIn("Below 70%", html)
        # these fixtures all ship a config, so the run-parameters drawer is present.
        self.assertIn('id="params-btn"', html)
        self.assertIn('id="params-drawer"', html)
        # the assembly-stats table carries its search box + low-coverage filter,
        # a worst-first pre-sort, and visual status encoding.
        self.assertIn('id="statsSearch"', html)
        self.assertIn('id="lowCovToggle"', html)
        self.assertIn('id="statsRowCount"', html)
        self.assertIn("cov-dot", html)
        self.assertIn("cov-bar", html)
        # throughput is the 3-series stacked bar with an Absolute/Percent toggle.
        self.assertIn('id="throughput-scale"', html)
        self.assertIn("Removed by QC", html)
        self.assertIn("QC-passed, unmapped", html)
        # the coverage heatmap replaces the old aggregated line overlay.
        self.assertIn('id="vu-heatmap"', html)
        self.assertIn('id="heatmapData"', html)
        self.assertIn("Coverage heatmap", html)
        self.assertNotIn('id="aggregated-card"', html)

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
        # KPI tiles carry the actual server-rendered values (not just the keys):
        # 2 samples, 1 at >=90% (sample-A 0.995; sample-B 0.70 warns but is not
        # below the 0.70 warn cutoff, so below_warn=0), median coverage
        # median(0.995,0.70)=0.8475 -> 84.7%, mean depth mean(512.4,45.7)=279.05
        # -> 279.1x (with the x suffix). pass_sub = 1/2 -> "50% of run".
        self.assertIn('data-kpi="samples">2<', html)
        self.assertIn('data-kpi="pass_count">1<', html)
        self.assertIn('data-kpi="below_warn">0<', html)
        self.assertIn('data-kpi="pass_sub">50% of run<', html)
        self.assertIn('data-kpi="median_coverage">84.7%<', html)
        self.assertIn('data-kpi="mean_depth">279.1×<', html)

    def test_nanopore_single_end_snapshot(self):
        html = _render("nanopore")
        self._assert_common_invariants(html)
        # single-end run: no doubling, and the report says QC was not run.
        self.assertIn("single-end", html)
        self.assertIn("no read-QC step", html)
        # 256410 / 462000 = 55.5% (undoubled denominator).
        self.assertIn("55.5%", html)
        self.assertIn("462,000", html)
        # KPI tiles: single sample at 0.982 breadth / 845.3x depth.
        self.assertIn('data-kpi="median_coverage">98.2%<', html)
        self.assertIn('data-kpi="mean_depth">845.3×<', html)
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
        # segmented runs get the Global | Per-segment KPI switch.
        self.assertIn('id="kpi-scope"', html)
        self.assertIn("Per segment", html)
        # the parameters drawer shows the config verbatim, incl. a Resources group.
        self.assertIn("minimum_depth", html)
        self.assertIn(">42<", html)  # the fixture's distinctive minimum_depth value
        self.assertIn(">Resources</h3>", html)
        # length-weighted whole-genome KPIs land in the right tiles: one sample,
        # weighted breadth (0.97,0.92 over equal-length tracks) = 0.945 -> 94.5%,
        # weighted depth (220.5,88.2) = 154.35 -> 154.3x.
        self.assertIn('data-kpi="samples">1<', html)
        self.assertIn('data-kpi="median_coverage">94.5%<', html)
        self.assertIn('data-kpi="mean_depth">154.3×<', html)
        # the "N segments each" subtitle appears under the Samples tile.
        self.assertIn("2 segments each", html)

    def test_unsegmented_has_no_per_segment_kpi_switch(self):
        # a single-genome run has nothing to split by, so no Global/Per-segment toggle.
        self.assertNotIn('id="kpi-scope"', _render("unsegmented"))

    def _annotation_json(self, html):
        m = re.search(r'id="annotationData"[^>]*>(.*?)</script>', html, re.S)
        self.assertIsNotNone(m, "annotationData block missing")
        return json.loads(m.group(1))

    def test_annotation_tracks_present_when_staged(self):
        html = _render("annotated")
        ann = self._annotation_json(html)
        self.assertTrue(ann, "expected a non-empty annotation model")
        labels = [ln["label"] for seg in ann.values() for ln in seg["lanes"]]
        kinds = [ln["kind"] for seg in ann.values() for ln in seg["lanes"]]
        self.assertIn("Genes", labels)
        self.assertTrue(any(lbl.startswith("Pool") for lbl in labels))
        self.assertIn("gene", kinds)
        self.assertIn("primer", kinds)
        # the client draws the tracks on a second y-axis via a shared helper.
        app = _app_js(html)
        self.assertIn("annotationTraces", app)
        self.assertIn("yaxis2", app)
        # Genes/Primers on-off chips are wired (built from the present kinds).
        self.assertIn("buildTrackToggles", app)
        self.assertIn("annotKinds", app)
        # track hues are theme-driven from server-set body attributes.
        self.assertIn("data-track-genes", html)
        self.assertIn("data-track-primers", html)

    def test_no_annotation_fixtures_have_empty_annotation_block(self):
        # Regression fence: runs without staged annotation embed an empty model,
        # draw no tracks, and stay byte-compatible with the pre-feature report.
        for fixture in ("unsegmented", "nanopore", "segmented"):
            html = _render(fixture)
            self.assertEqual(self._annotation_json(html), {}, fixture)

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
