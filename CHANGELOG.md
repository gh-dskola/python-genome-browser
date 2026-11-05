# Changelog

All notable changes to pygbrowse are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions come from git tags via [setuptools-scm](https://github.com/pypa/setuptools-scm).

See [CONTRIBUTING.md](CONTRIBUTING.md#versioning-policy) for the versioning policy:
this project is on `0.x`, so a MINOR bump may contain breaking changes. Every break is
listed under a `### Breaking changes` heading with a *Migration:* note.

---

## Unreleased

### Added

- `BigWigDataSource` reads signal tracks from BigWig files (`pyBigWig` dependency).
- `pyproject.toml` with setuptools-scm versioning, `environment.yaml`, ruff lint config,
  `.pre-commit-config.yaml`.
- Docs: `CLAUDE.md`, `CONTRIBUTING.md`, `PAPERCUTS.md`.

### Removed

- `setup.py` and `build_and_upload.sh` (no PyPI publishing).

## [0.3.3] - 2018-10-12

### Added

- `GeneModelPlot(genes_include=...)` restricts the plot to the named genes.

### Fixed

- Multiple intervals with the same name are all drawn.

## [0.3.2] - 2018-09-18

### Fixed

- The `alpha` parameter of `WigPlot` takes effect.

## [0.3.1] - 2018-08-09

### Fixed

- Initial exons and five-prime UTRs were missed on + strand genes.

## [0.3.0] - 2018-08-05

### Added

- Automatic y-axis scaling and improved y-tick labels for `BedPlot`.

### Changed

- `WigPlot` is solid-filled by default.

## [0.2.0] - 2018-08-02

### Changed

- Gene model backend reads GFF3 through tabix, for speed and memory.
- Transcripts can be filtered by support level (default: all).

## [0.1.0] - 2018-07-24

- Initial release.
