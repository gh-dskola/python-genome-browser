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
- `IntervalData.query(chrom, start, end)` returns the overlapping intervals.
- `TagDirectory(read_handling="starts")` counts read starts only.
- `GenomeBrowser.visualize(style=...)` overrides the track style for one figure.
- Data sources accept `pathlib.Path` as well as `str` paths.
- pytest suite on synthetic data (`tests/`).

### Breaking changes

- **seaborn is no longer used.** `GenomeBrowser.visualize` applies its style inside a
  `matplotlib.rc_context` and leaves global rcParams alone.
  *Migration:* replace `seaborn_style=seaborn.axes_style(...)` with `style={<rcParams>}`.
- **`WigPlot(smoothing_bandwidth=...)` removed**; it duplicated the data-source option.
  *Migration:* pass `smoothing_bandwidth` to the data source (`BigWigDataSource`, `TagDirectory`).
- **`Gff3Annotations(outgoing_chromosome_name_converter=...)` removed**; it was never used.
  *Migration:* drop the argument.
- **`HicDataDir` has no default file template, and `binsize` is now `bin_size`.**
  *Migration:* pass the template explicitly; rename the keyword.
- **`GeneModelPlot(genes_include=...)` matches whole gene names**, not substrings.
  *Migration:* list every gene name to show.

### Fixed

- `convert_chromosome_name` no longer reads `X`, `V`, `XI` etc. as Roman numerals outside
  the `yeast` dialect (`X` → `chr10` before); `chrM` ↔ `MT` replaces a misspelled Drosophila
  name.
- Tracks no longer fall off the bottom of the figure at six or more rows.
- `GeneModelPlot` draws only each gene's own UTRs and CDSs (it drew every overlapping
  gene's components on every row), lays out only the genes in `genes_include`, and sizes
  + strand 3' UTRs in data units.
- `BedPlot` draws intervals that span the whole window, and `patch_kwargs` no longer
  change the class default.
- `GenomeBrowser.visualize` returns a new figure on every call.
- `HicPlot` no longer writes into the source matrix (it failed under pandas Copy-on-Write)
  and masks diagonals on both sides.
- `HicDataDir` queries work (they raised `NameError`).
- `TagDirectory`: no infinite loop at end of file; works with `normalize_to=None`; queries
  past the last tag return zeros; read coverage is half-open (it covered one extra base).
- `utilities.binary_search_tag_file` finds the first tag at or after the target (it returned
  `None` for most targets).
- `SparseVectors` queries starting before the first point keep the leading points.
- `WigPlot` scaling and centering no longer modify the data source's vectors.
- `IntervalData`: named DataFrames without a `strand` column keep their columns; BED files
  with 3 to 6+ columns load; an unknown `format` raises `ValueError`.
- `Gff3Annotations` parses trailing `;`, `=` inside values, %-escapes, comma-separated
  `Parent`, and genes without `Name` (the ID is used); orphan features are logged.
- `bgzip_gff` / `generate_tabix_index` use pysam instead of a shell pipeline and raise on
  failure instead of printing.

### Removed

- `romannumerals.py`; `yeast` conversion uses a 16-entry table.
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
