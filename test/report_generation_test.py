"""End-to-end test for the interactive consensus HTML report.

Unlike the inline-fixture unit tests under test/scripts/, this exercises the full
render against realistic nested output trees under test/fixtures/report/ (an
``assembly/coverage_stats/...`` layout and its segmented ``assembly/{segment}/...``
variant), matching resolve_basewise_path's exact structure.
"""

import os
import re
import tempfile
import unittest

from plotly.offline import get_plotlyjs

from viralunity.scripts.python.generate_consensus_report import write_report

FIXTURE_ROOT = os.path.join(os.path.dirname(__file__), "fixtures", "report")
EXTERNAL_RESOURCE_TAG = re.compile(r"<(?:script[^>]+src|link[^>]+href)\s*=", re.IGNORECASE)


class _ReportCase:
    """Shared assertions; subclasses set fixture_dir + expected_samples."""

    fixture_dir = None
    expected_samples = ()

    def _render(self):
        with tempfile.TemporaryDirectory() as d:
            dest = os.path.join(d, "report.html")
            write_report(self.fixture_dir, dest)
            self.assertTrue(os.path.isfile(dest))
            # exactly one output file is produced.
            self.assertEqual(os.listdir(d), ["report.html"])
            with open(dest) as fh:
                return fh.read()

    def test_self_contained_offline(self):
        html = self._render()
        # No external resource is loaded: no <script src=>/<link href=> tags.
        self.assertIsNone(EXTERNAL_RESOURCE_TAG.search(html))
        # Plotly is inlined exactly once; removing that single vendored block,
        # OUR content must carry no http(s) URL. (The spec's literal "no http
        # anywhere" is unattainable because plotly.js embeds xmlns/map URLs.)
        remainder = html.replace(get_plotlyjs(), "", 1)
        self.assertNotIn("http://", remainder)
        self.assertNotIn("https://", remainder)

    def test_contains_sample_names_verbatim(self):
        html = self._render()
        for sample in self.expected_samples:
            self.assertIn(sample, html)

    def test_contains_coverage_guides(self):
        html = self._render()
        self.assertIn("20x", html)
        self.assertIn("100x", html)

    def test_has_prerendered_plotly_figures(self):
        html = self._render()
        # reads histogram + depth bar + >=1 aggregated panel are pre-rendered.
        self.assertGreaterEqual(html.count("Plotly.newPlot"), 3)


class TestUnsegmentedReport(_ReportCase, unittest.TestCase):
    fixture_dir = os.path.join(FIXTURE_ROOT, "unsegmented")
    expected_samples = ("sample-A", "sample-B")


class TestSegmentedReport(_ReportCase, unittest.TestCase):
    fixture_dir = os.path.join(FIXTURE_ROOT, "segmented")
    expected_samples = ("sample-A",)

    def test_segment_panels_present(self):
        html = self._render()
        # one aggregated panel per segment, toggled by the segment <select>.
        self.assertIn('data-key="S1"', html)
        self.assertIn('data-key="S2"', html)


if __name__ == "__main__":
    unittest.main()
