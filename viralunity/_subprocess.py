"""Shared subprocess helper for ViralUnity CLI commands.

Internal module — not part of the public API. Used by the `get-databases`
and `build-deacon-index` subcommands to invoke external tools (`wget`,
`tar`, `diamond`, `deacon`, `datasets`, ...) with consistent echoing and
error handling.
"""

from __future__ import annotations

import subprocess
from typing import Optional, Sequence

import click

#: Default wall-clock ceiling for a single external command (seconds). Long
#: enough for large downloads/DB builds, but bounded so a hung tool cannot block
#: a run (or a service worker) indefinitely. Override per-call via ``timeout``.
DEFAULT_TIMEOUT = 6 * 60 * 60  # 6 hours


def run_command(
    cmd: Sequence[object],
    cwd: Optional[str] = None,
    timeout: Optional[float] = DEFAULT_TIMEOUT,
) -> None:
    """Run a subprocess command, streaming output and raising on failure.

    Args:
        cmd: Argument list (no shell). Each element is coerced to ``str``
            when echoed; pass paths and flags as separate elements.
        cwd: Optional working directory in which to run the command.
        timeout: Wall-clock ceiling in seconds. Defaults to
            ``DEFAULT_TIMEOUT``; pass ``None`` to disable.

    Raises:
        click.ClickException: If the subprocess exits non-zero or times out.
            The exit code/timeout and command are included in the message so
            the failure surfaces cleanly in the CLI output.
    """
    click.echo(f"$ {' '.join(str(c) for c in cmd)}")
    try:
        result = subprocess.run(list(cmd), cwd=cwd, timeout=timeout)
    except subprocess.TimeoutExpired as e:
        raise click.ClickException(
            f"Command timed out after {e.timeout:.0f}s: " f"{' '.join(str(c) for c in cmd)}"
        ) from e
    if result.returncode != 0:
        raise click.ClickException(
            f"Command failed with exit code {result.returncode}: "
            f"{' '.join(str(c) for c in cmd)}"
        )
