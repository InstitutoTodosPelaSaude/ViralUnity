"""Guard the report generator's minimal-environment constraint.

``generate_consensus_report.py`` runs under Snakemake's ``envs/report.yaml`` conda
env, which installs only pandas/plotly/jinja2 — NOT the ``viralunity`` package and
NOT PyYAML. The module therefore must (a) never ``import viralunity`` at module
scope and (b) import ``yaml`` lazily, only on the CLI ``--config-file`` path.

The rest of the suite runs in the full dev env (both available), so a regression
that hoisted either import would pass every other test yet crash the real
``generate_html_report`` rule at runtime. This test reproduces the constraint by
importing the module as a standalone file with ``viralunity`` and ``yaml`` blocked
from the import system, then rendering a fixture (no config -> no yaml needed).
"""

import os
import subprocess
import sys
import unittest

from viralunity.scripts.python import generate_consensus_report as _gcr

MODULE_PATH = _gcr.__file__
FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "report", "unsegmented")

# Child program: block the two modules the report.yaml env lacks, load the
# generator by path, and render a report. Any hoisted import trips exec_module.
_CHILD = r"""
import importlib.util
import sys

class _Blocker:
    def __init__(self, names):
        self.names = set(names)
    def find_spec(self, fullname, path, target=None):
        if fullname.split(".")[0] in self.names:
            raise ModuleNotFoundError(fullname + " is blocked (report.yaml env)")
        return None

BLOCKED = {"viralunity", "yaml"}
sys.meta_path.insert(0, _Blocker(BLOCKED))
for _m in [m for m in sys.modules if m.split(".")[0] in BLOCKED]:
    del sys.modules[_m]

mod_path, fixture = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location("gcr_standalone", mod_path)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)                       # must not import viralunity/yaml

meta = mod.build_report_metadata(None, fixture)    # config=None -> yaml never touched
html = mod.render_report(fixture, meta, None)
assert 'id="kpi-grid"' in html, "report did not render"
print("OK")
"""


class TestReportEnvConstraint(unittest.TestCase):
    def test_imports_and_renders_without_viralunity_or_yaml(self):
        proc = subprocess.run(
            [sys.executable, "-c", _CHILD, MODULE_PATH, FIXTURE],
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            proc.returncode,
            0,
            msg=f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}",
        )
        self.assertIn("OK", proc.stdout)

    def test_yaml_is_only_imported_lazily_on_the_config_path(self):
        # The single ``import yaml`` must live inside load_run_config, not at
        # module scope (checked structurally as a fast, env-independent guard).
        with open(MODULE_PATH) as fh:
            lines = fh.read().splitlines()
        module_scope_yaml = [ln for ln in lines if ln.startswith(("import yaml", "from yaml"))]
        self.assertEqual(module_scope_yaml, [], "yaml must not be imported at module scope")
        self.assertTrue(
            any("import yaml" in ln for ln in lines),
            "expected a lazy 'import yaml' inside the CLI config loader",
        )


if __name__ == "__main__":
    unittest.main()
