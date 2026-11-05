# Contributing to pygbrowse

## Quality gate

Both of these must pass before you commit. Run them in the `pygbrowse` mamba environment
(`mamba env create -f environment.yaml`):

```bash
mamba activate pygbrowse
ruff check pygbrowse tests
python -m pytest tests/
```

`ruff` is the single linter and formatter; its config lives in `pyproject.toml`. Do not
introduce a second linter. `pre-commit install` wires ruff into `git commit`.

## Changelog

Update [CHANGELOG.md](CHANGELOG.md) in the same commit as the change. New entries go under
`## Unreleased`, in the Keep a Changelog category that fits (`Added`, `Changed`,
`Deprecated`, `Removed`, `Fixed`, `Security`). If your change breaks a caller, it also
needs an entry under `### Breaking changes` with a *Migration:* note.

## Versioning policy

This project follows [Semantic Versioning 2.0.0](https://semver.org/spec/v2.0.0.html)
and is at **major version 0**, which under semver §4 means the public API is not yet
declared stable.

- **`0.x` releases may contain breaking changes in a MINOR bump.**
- **Every breaking change is listed under a `### Breaking changes` heading** in that
  release's section, with a *Migration:* note per item. Treat its absence as a claim that
  the release is backward compatible.
- **PATCH bumps are reserved for changes that cannot break a caller**: bug fixes that do
  not alter documented behavior, and internal or developer-tooling changes.
- **Public** means anything importable from `pygbrowse` or its modules without a leading
  underscore.
- **Deprecation is best-effort while on `0.x`.** Pin an exact tag if you need stability.

### Tag format

Tags are `vMAJOR.MINOR.PATCH`, annotated, e.g. `v0.4.0`.

## Cutting a release

1. Confirm the working tree is clean and the quality gate passes.
2. In `CHANGELOG.md`, rename `## Unreleased` to `## [MAJOR.MINOR.PATCH] - YYYY-MM-DD`
   and open a fresh empty `## Unreleased` above it.
3. Check the bump size against the policy above.
4. Commit, then tag: `git tag -a vX.Y.Z -m "vX.Y.Z - <summary>"`.
5. `pip install -e .` and confirm `pygbrowse.__version__` reports the new version.
6. Push the branch and the tag. A pushed tag should not be moved or deleted.
