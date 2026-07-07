# Releasing ViralUnity

This document describes how to cut a versioned release. Only project
maintainers need to read it.

## Versioning

ViralUnity follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
The version string is single-sourced from
[`viralunity/__init__.py`](viralunity/__init__.py): the value of
`__version__` is read dynamically by `pyproject.toml` through
`[tool.setuptools.dynamic]`, so updating it in one place propagates to
the Python package metadata.

The `Dockerfile`'s `LABEL version="X.Y.Z"` is not derived
automatically — it must be bumped in lockstep.

## Cutting a release (X.Y.Z)

1. **Bump `__version__`** in `viralunity/__init__.py` to `X.Y.Z`.

2. **Bump the Dockerfile LABEL** at `Dockerfile:2` to
   `LABEL version="X.Y.Z"`. Keep the `description` line aligned with
   `viralunity._description` (this should already match).

3. **Update `CHANGELOG.md`**: rename the most recent unreleased section
   (or add a new one) to `## [X.Y.Z] - YYYY-MM-DD`, grouping entries
   under `### Added` / `### Changed` / `### Fixed` / `### Refactored`.

4. **Commit** the three edits as a single commit:

   ```bash
   git add viralunity/__init__.py Dockerfile CHANGELOG.md
   git commit -m "release: X.Y.Z"
   ```

5. **Tag**:

   ```bash
   git tag -a vX.Y.Z -m "Release vX.Y.Z"
   ```

6. **Push** the branch and the tag:

   ```bash
   git push origin main
   git push origin vX.Y.Z
   ```

   Pushing the `vX.Y.Z` tag triggers the `release.yaml` workflow, which builds
   the sdist + wheel and **publishes to PyPI** via trusted publishing. Stage on
   TestPyPI first — see [Publishing to PyPI](#publishing-to-pypi) below.

7. **Build and publish the Docker image** (if applicable):

   ```bash
   docker build -t institutotodospelasaude/viralunity:X.Y.Z .
   docker tag institutotodospelasaude/viralunity:X.Y.Z \
              institutotodospelasaude/viralunity:latest
   docker push institutotodospelasaude/viralunity:X.Y.Z
   docker push institutotodospelasaude/viralunity:latest
   ```

## Publishing to PyPI

Publishing is automated by [`.github/workflows/release.yaml`](.github/workflows/release.yaml)
using **PyPI trusted publishing** (OIDC) — no API tokens are stored in the repo.

### One-time setup (per maintainer account)

Before the first release, register this repository as a *trusted publisher* on
both indexes and create the matching GitHub deployment environments:

1. On <https://test.pypi.org> and <https://pypi.org>, add a *pending publisher*
   under the account/project publishing settings:
   - Owner: `InstitutoTodosPelaSaude`, repo: `ViralUnity`
   - Workflow: `release.yaml`
   - Environment: `testpypi` (on test.pypi.org) / `pypi` (on pypi.org)
2. In GitHub → Settings → Environments, create environments named `testpypi`
   and `pypi` (add required reviewers on `pypi` if you want a manual approval
   gate before the real upload).

### Each release

1. **Stage on TestPyPI first.** Trigger `release.yaml` manually
   (Actions → Release → *Run workflow*, or
   `gh workflow run release.yaml`). The `workflow_dispatch` run uploads to
   TestPyPI. Verify a clean install:

   ```bash
   pip install --index-url https://test.pypi.org/simple/ \
     --extra-index-url https://pypi.org/simple/ viralunity==X.Y.Z
   viralunity --version   # -> viralunity, version X.Y.Z
   ```

2. **Promote to PyPI.** Push the `vX.Y.Z` tag (step 6 above). The tag push runs
   the `publish-pypi` job, which uploads the same build to PyPI. The first real
   upload of a version cannot be redone, which is why TestPyPI staging comes
   first.

## Verifying

After the tag is pushed:

```bash
# PyPI install (once the publish-pypi job completes)
pip install viralunity==X.Y.Z
# Python package version
viralunity --version
# -> viralunity, version X.Y.Z

# Docker image labels
docker inspect institutotodospelasaude/viralunity:X.Y.Z \
  | jq '.[0].Config.Labels'
# -> { "version": "X.Y.Z", "description": "A pipeline for viral metagenomics analysis." }
```

`viralunity.__version__`, the Dockerfile `LABEL version`, the git tag,
and the topmost `CHANGELOG.md` section should all agree.

## Hotfix releases

For a patch release on top of a tag (e.g. `v1.1.0` → `v1.1.1`):

1. Branch from the tag: `git checkout -b release/1.1.1 v1.1.0`.
2. Apply the fix, then follow steps 1–7 above with `X.Y.Z = 1.1.1`.
3. Open a PR back into `main` so the fix lands there too.

## Backfilling a missed tag

If a release commit landed on `main` but no tag was created, the tag
can be added retroactively at that commit:

```bash
git tag -a vX.Y.Z <commit-sha> -m "Release vX.Y.Z"
git push origin vX.Y.Z
```
