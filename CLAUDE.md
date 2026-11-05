# CLAUDE.md

pygbrowse draws matplotlib figures of genomic data as stacked browser tracks. Data
sources (`pygbrowse/datasources.py`) read files and answer region queries. Plot classes
(`pygbrowse/plots.py`) draw one track each. `GenomeBrowser` stacks the tracks into one
figure.

## Environment

```bash
mamba env create -f environment.yaml   # creates env `pygbrowse`, installs the package editable
mamba run -n pygbrowse python -m pytest tests/
mamba run -n pygbrowse ruff check pygbrowse tests
```

Both commands must pass before a commit. The gate is 0 failures and 0 skips. Tests build
all their data (GFF3, BigWig, HOMER tag directory, BED) in `tmp_path`, so no test needs a
file from outside the repo.

## Coordinate convention

- Every `query(chrom, start, end)` takes and returns **0-based, half-open** coordinates.
- Conversions from 1-based formats live in one place each, with a comment:
  - GFF3 (1-based, closed) → 0-based half-open in `Gff3Annotations`.
  - HOMER tag files (1-based starts) → 0-based in `TagDirectory` and
    `utilities.binary_search_tag_file`.
- BED and BigWig are already 0-based, half-open.
- Chromosome names in queries are UCSC style (`chr1`, `chrX`, `chrM`). Data sources
  convert to the file's dialect with `utilities.convert_chromosome_name`.

## Layout

- `datasources.py`: vector sources (`BigWigDataSource`, `TagDirectory`, `SparseVectors`),
  `IntervalData`, `Gff3Annotations`, Hi-C matrices (`HicDataDir`, `HicDataDict`).
- `plots.py`: `GenomeBrowser`, `GeneModelPlot`, `WigPlot`, `BedPlot`, `InteractionPlot`,
  `HicPlot`.
- `_arcs.py`: polygon arcs for `InteractionPlot`. matplotlib's `Arc` patch is not used.
- `utilities.py`: chromosome-name conversion, tag-file binary search, GFF3 bgzip/tabix.
- `pygbrowse_demonstration.ipynb`: demo notebook. Its data files are not in the repo, so
  it cannot be run as a check.

## Conventions

- Logging via `logging.getLogger(__name__)`; no `print`.
- ruff with `C901` at `max-complexity = 10`.
- Update `CHANGELOG.md` (`## Unreleased`) with every change. Breaking changes need a
  `### Breaking changes` entry with a *Migration:* note (see `CONTRIBUTING.md`).

## Gotchas

- `python -m pytest` puts the current directory first on `sys.path`, so run from the repo
  root it always imports the working-tree package. To test an old checkout, run pytest
  from another directory with `PYTHONPATH` pointing at that checkout.
- `pygbrowse.__version__` comes from git tags via setuptools-scm; without an installed
  package (plain `PYTHONPATH=.`) it reads `0.0.0.dev0`.
