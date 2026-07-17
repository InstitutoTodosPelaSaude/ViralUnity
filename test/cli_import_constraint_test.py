"""Guard the base CLI against eagerly importing report-only dependencies.

The interactive HTML report needs ``plotly`` and ``jinja2``, but the pipeline
builds it in a separate ``envs/report.yaml`` conda env — the base ``viralunity``
env does not need them. They must therefore be imported lazily (inside the
``report`` command callback), so that importing ``viralunity.viralunity_cli`` and
running the other subcommands (``create-samplesheet``, ``meta``, ``consensus``,
``setup`` …) works even when plotly/jinja2 are absent.

A regression that hoisted the report import back to module scope would pass every
other test (the dev env has plotly) yet crash every ``viralunity`` invocation in a
pre-1.4.0 environment. This test reproduces that environment by blocking plotly and
jinja2 from the import system, then importing the CLI and exercising a non-report
command.
"""

import subprocess
import sys
import unittest

# Child program: block the report-only modules, import the CLI, and drive a
# non-report subcommand through click's test runner. Any hoisted plotly/jinja2
# import trips the import of viralunity_cli.
_CHILD = r"""
import sys

class _Blocker:
    def __init__(self, names):
        self.names = set(names)
    def find_spec(self, fullname, path, target=None):
        if fullname.split(".")[0] in self.names:
            # Set `name` like a real missing-module error so code that inspects
            # `e.name` (e.g. the report command's dependency hint) behaves as it
            # would in a genuinely plotly-less environment.
            raise ModuleNotFoundError(
                fullname + " is blocked (report-only dep)", name=fullname
            )
        return None

BLOCKED = {"plotly", "jinja2"}
sys.meta_path.insert(0, _Blocker(BLOCKED))
for _m in [m for m in sys.modules if m.split(".")[0] in BLOCKED]:
    del sys.modules[_m]

from viralunity.viralunity_cli import cli  # must not import plotly/jinja2
assert "plotly" not in sys.modules, "importing the CLI eagerly imported plotly"
assert "jinja2" not in sys.modules, "importing the CLI eagerly imported jinja2"

from click.testing import CliRunner
for sub in ("meta", "create-samplesheet", "consensus"):
    result = CliRunner().invoke(cli, [sub, "--help"])
    assert result.exit_code == 0, sub + " --help failed:\n" + result.output

# report --help must also work without the deps (the callback never runs)...
assert CliRunner().invoke(cli, ["report", "--help"]).exit_code == 0, "report --help failed"
# ...but actually running report without the deps must fail with the friendly hint,
# not a raw ModuleNotFoundError traceback.
result = CliRunner().invoke(cli, ["report", "--input", "/nonexistent/vu/dir"])
assert result.exit_code != 0, "report should fail when deps are missing"
assert "dependencies are missing" in result.output, (
    "expected the friendly missing-deps message, got:\n" + result.output
)

print("OK")
"""


class TestCliImportConstraint(unittest.TestCase):
    def test_cli_imports_and_runs_without_plotly_or_jinja2(self):
        proc = subprocess.run(
            [sys.executable, "-c", _CHILD],
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            proc.returncode,
            0,
            msg=f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}",
        )
        self.assertIn("OK", proc.stdout)

    def test_report_deps_are_not_imported_at_module_scope(self):
        # The report generator import must live inside the report() callback, not
        # at viralunity_report module scope (fast, env-independent structural guard).
        from viralunity import viralunity_report as _vr

        with open(_vr.__file__) as fh:
            lines = fh.read().splitlines()
        offending = [
            ln
            for ln in lines
            if ln.startswith(("import ", "from ")) and "generate_consensus_report" in ln
        ]
        self.assertEqual(
            offending,
            [],
            "generate_consensus_report must be imported lazily inside report(), "
            "not at module scope",
        )


if __name__ == "__main__":
    unittest.main()
