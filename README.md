![alt pygbrowse](pygbrowse_logo_1_med_flat.png)

# python-genome-browser
Tools for making plots of genomic datasets in a genome-browser-like format.

## What does it do?

The python genome browser, AKA pygbrowse, makes it easy to generate aligned plots of genomic data. This is similar to software such as the [Integrated Genome Viewer (IGV) from the Broad Institute](http://software.broadinstitute.org/software/igv/) or webserver-based public genome browsers such as the [UCSC genome browser](http://genome.ucsc.edu) or the [WashU Epigenome Browser](http://epigenomegateway.wustl.edu/).

## Typical use cases:

1. Visualizing data inside a Jupyter notebook.
2. Preparing figures for presentation or publication.

## Installation

```bash
git clone https://github.com/guardant/python-genome-browser.git
cd python-genome-browser
mamba env create -f environment.yaml   # creates the `pygbrowse` env and installs the package
mamba activate pygbrowse
```

Requires Python 3.10 or later. `pip install -e .` inside any environment with the
dependencies in `pyproject.toml` also works.

## Quick start

```python
from pygbrowse import datasources, plots, utilities

# Gene models need a sorted, bgzipped, tabix-indexed GFF3 (Ensembl layout).
utilities.bgzip_gff("Mus_musculus.GRCm38.93.chr.gff3.gz", "genes.gff3.bgz")
utilities.generate_tabix_index("genes.gff3.bgz")

genes = plots.GeneModelPlot(data=datasources.Gff3Annotations("genes.gff3.bgz"))
signal = plots.WigPlot(data=datasources.BigWigDataSource("signal.bw"), label="ChIP-seq")
peaks = plots.BedPlot(data=datasources.IntervalData("peaks.bed"), label="Peaks")

browser = plots.GenomeBrowser([[signal], [peaks], [genes]])
fig = browser.visualize("chr4", 62494071, 62546083)
fig.savefig("region.png", dpi=300, bbox_inches="tight")
```

Each inner list is one row; tracks in the same row share an Axes. `visualize` returns a
new `matplotlib.figure.Figure` and leaves matplotlib's global settings unchanged. Pass
`style={...}` (rcParams) to change the look of one figure. To see progress messages,
call `logging.basicConfig(level=logging.INFO)`.

## Coordinates and chromosome names

- Every `query(chrom, start, end)` and every `visualize(chrom, start, end)` window is
  **0-based, half-open** (`[start, end)`), as in BED.
- Query chromosomes are UCSC style (`chr1`, `chrX`, `chrM`). Sources convert to the file's
  naming: `Gff3Annotations` expects Ensembl contigs (`1`, `X`, `MT`) by default, and
  `BigWigDataSource(chromosome_dialect=...)` accepts `ucsc`, `ensembl` or `yeast`.
- 1-based formats are converted on read: GFF3 (1-based, closed) in `Gff3Annotations`, HOMER
  tag files (1-based starts) in `TagDirectory`.

## Supported formats

| Data | Source class | Track class |
|---|---|---|
| Gene models: GFF3, bgzipped + tabix-indexed | `Gff3Annotations` | `GeneModelPlot` |
| Signal: BigWig | `BigWigDataSource` | `WigPlot` |
| Signal: HOMER tag directory | `TagDirectory` | `WigPlot` |
| Signal: in-memory `pandas.Series` per chromosome | `SparseVectors` | `WigPlot` |
| Intervals: BED (3+ columns), HOMER peaks, HOMER annotated peaks, DataFrame | `IntervalData` | `BedPlot` |
| Pairwise interactions: DataFrame with `chr1, start1, end1, start2, end2` | (the DataFrame) | `InteractionPlot` |
| Hi-C matrices: per-chromosome TSV files or in-memory DataFrames | `HicDataDir`, `HicDataDict` | `HicPlot` |

## Testing

```bash
mamba activate pygbrowse
ruff check pygbrowse tests
python -m pytest tests/
```

The tests build all their data (GFF3, BigWig, HOMER tag directory, BED) in a temporary
directory. See [CONTRIBUTING.md](CONTRIBUTING.md) for the changelog and versioning rules.

## Why not just use those tools you just mentioned?

For a couple of reasons:

1. Speed and ease of use: The workflow for converting data to the appropriate formats, uploading the data to a remote webserver (or providing URLs for hosted data), then selecting the region to visualize and waiting for it to be rendered is slow and cumbersome with many manual interventions needed. Although APIs exist for some of these tools, creating a flexible end-to-end automated visualizer that takes arbitrary data and displays it on your screen would require quite a bit of custom scripting.
2. Flexibility: Different genome browsers natively accept different subsets of the data types commonly used in genomics, and interconverting them is tedious. In some cases certain data types are not supported at all (e.g. long-range interactions such as HiC are not suported by the UCSC genome browser). Pygbrowse (will) natively support most of the most common data formats, removing one or more "data munging" operations from your workflow.
3. Transparency: Current genome browsers have very specific requirements for the format of their input data, not all of which are as well-documented as we might like. In addition, very little feedback is provided regarding such errors, leading the user into a painful trial-and-error process in order to load their data. Pygbrowse strives to be more flexible in the format of its input data and to provide helpful feedback when problems do occur. 
4. Beauty: They say that beauty is in the eye of the beholder but to _this_ beholder the default outputs of the available genome browsers are aesthetically-lacking. They often require extensive manipulation in Adobe Illustrator or Ink to prepare them for publication or even to be legible in a presentation slide. Pygbrowse is designed from the ground up for generating static figures with proportions that can be easily scaled for common use cases.

## Wait, you're calling this a genome browser but I can't really browse around like I can with those other tools, can I?

Well, you _can_ browse, in a sense, by calling the visualizer with different sets of coordinates. But it's not really designed for the kind of interactive, dynamic browsing experience provided by the other tools -- as mentioned earlier it is optimzed for rendering static "browser snapshots". That being said, however, we may add controls to provide interactive browser functionality when used inside a Jupyter notebook in a future version . . .

Logo makes use of clipart by Paula Washington at [AnimalsClipart.com](http://animalsclipart.com/pig-design/)
