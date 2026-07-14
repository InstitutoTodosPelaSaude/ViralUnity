#!/usr/bin/env python3
"""Empirical before->after analysis of the consensus report UX pass.

Renders the committed report fixtures (``test/fixtures/report/``) with BOTH the
pre-UX generator (extracted from git at the ``PRE_UX_REF`` commit) and the current
generator, then extracts concrete, checkable signals from each rendered report so
the improvement is quantitative rather than only visual. Regenerate the results
table in ``EMPIRICAL_ANALYSIS.md`` with:

    python report-ux-improvements/empirical_analysis.py

Writes the rendered before/after HTML reports to ``--out`` (default: a temp dir)
so they can be opened and compared in a browser.
"""

import argparse
import importlib.util
import json
import os
import re
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIX = os.path.join(REPO, "test", "fixtures", "report")
PRE_UX_REF = "6277c03"  # last commit before the UX/data-viz pass
MODULE = "viralunity/scripts/python/generate_consensus_report.py"
TEMPLATE = "viralunity/scripts/python/templates/report_template.html.j2"

sys.path.insert(0, REPO)
from viralunity.scripts.python.generate_consensus_report import (  # noqa: E402
    build_report_metadata,
    load_run_config,
    render_report,
)


def load_pre_ux_generator(workdir):
    """Materialize the pre-UX generator + template from git and import it."""
    pkg = os.path.join(workdir, "old")
    os.makedirs(os.path.join(pkg, "templates"), exist_ok=True)
    for src, dst in [(MODULE, "gen_old.py"), (TEMPLATE, "templates/report_template.html.j2")]:
        blob = subprocess.check_output(["git", "-C", REPO, "show", f"{PRE_UX_REF}:{src}"])
        with open(os.path.join(pkg, dst), "wb") as fh:
            fh.write(blob)
    spec = importlib.util.spec_from_file_location("gen_old", os.path.join(pkg, "gen_old.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _app_js(html):
    return re.findall(r"<script>(.*?)</script>", html, re.S)[-1]


def _coverage(html):
    m = re.search(r'id="coverageData"[^>]*>(.*?)</script>', html, re.S)
    return json.loads(m.group(1)) if m else {}


def signals(html):
    app = _app_js(html)
    first = next(iter(_coverage(html).values()), [])
    y0 = first[0]["y"] if first else []
    body = re.sub(r"<style>.*?</style>", "", html, flags=re.S).replace(" ", "")
    return {
        "grouped thousands sep": bool(re.search(r">\d{1,3}(,\d{3})+<", html)),
        "mapped shown as %": ">Mapped %<" in html,
        "baked (log) axis title": "Depth (log)" in html,
        "average-depth chart": "Average depth" in html,
        "'Sequencing throughput'": "Sequencing throughput" in html,
        "zero clamped to 1 (JS)": "Math.max(1," in app,
        "honest-zero gap (JS)": "connectgaps: false" in app,
        "toggle via react (JS)": "Plotly.react" in app,
        "yaxis.type relayout flip": ("'yaxis.type'" in app) or ('"yaxis.type"' in app),
        "By sample accordion": "By sample" in app,
        "By segment accordion": "By segment" in app,
        "fixed 900px width": "width:900px" in body,
        "min depth in coverage": (min(y0) if y0 else None),
    }


def mapped_cells(html):
    body = re.search(r"<tbody>(.*?)</tbody>", html, re.S).group(1)
    return re.findall(r'title="[^"]*">([^<]+)<', body)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None, help="Directory for rendered before/after HTML.")
    args = ap.parse_args()
    out = args.out or tempfile.mkdtemp(prefix="report_ux_")
    os.makedirs(out, exist_ok=True)

    with tempfile.TemporaryDirectory() as work:
        old = load_pre_ux_generator(work)
        cases = [
            (
                "unsegmented",
                {
                    "platform": "illumina",
                    "library_layout": "paired",
                    "primer_scheme": "schemes/sarscov2.primers.bed",
                    "qc_performed": True,
                },
            ),
            (
                "segmented",
                {
                    "platform": "illumina",
                    "library_layout": "paired",
                    "primer_scheme": None,
                    "qc_performed": True,
                },
            ),
        ]
        cols, table = [], {}
        for fixture, meta in cases:
            d = os.path.join(FIX, fixture)
            before, after = old.render_report(d), render_report(d, meta)
            open(os.path.join(out, f"before_{fixture}.html"), "w").write(before)
            open(os.path.join(out, f"after_{fixture}.html"), "w").write(after)
            for tag, html in [(f"before {fixture}", before), (f"after {fixture}", after)]:
                cols.append(tag)
                table[tag] = signals(html)
        nd = os.path.join(FIX, "nanopore")
        nano = render_report(
            nd, build_report_metadata(load_run_config(os.path.join(nd, "config.yml")), nd)
        )
        open(os.path.join(out, "after_nanopore.html"), "w").write(nano)
        cols.append("after nanopore")
        table["after nanopore"] = signals(nano)

    width = max(len(k) for k in table[cols[0]])
    header = f"{'signal':{width}} | " + " | ".join(f"{c:^16}" for c in cols)
    print(header)
    print("-" * len(header))
    for key in table[cols[0]]:
        print(f"{key:{width}} | " + " | ".join(f"{str(table[c][key]):^16}" for c in cols))
    print("\nMapping-rate cells (after):")
    print("  unsegmented:", mapped_cells(open(os.path.join(out, "after_unsegmented.html")).read()))
    print("  nanopore   :", mapped_cells(nano))
    print(f"\nRendered reports written to: {out}")


if __name__ == "__main__":
    main()
