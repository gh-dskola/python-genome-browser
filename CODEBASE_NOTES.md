# pygbrowse Codebase Notes

## Overview

**python-genome-browser (pygbrowse)** is a Python package for creating publication-quality genomic data visualizations in a genome-browser-like format, similar to IGV or UCSC Genome Browser but optimized for Jupyter notebooks and static figures.

- **Version**: 0.3.3
- **Author**: Dylan Skola
- **License**: MIT
- **Repository**: https://github.com/phageghost/python-genome-browser

## Architecture

The package follows a composition-based architecture with clear separation between data loading and visualization.

### File Structure

```
pygbrowse/
├── __init__.py           # Package initialization
├── plots.py              # Visualization classes (824 lines)
├── datasources.py        # Data source classes (406 lines)
├── utilities.py          # Helper functions (309 lines)
├── intervaloverlaps.py   # Interval overlap computation (74 lines)
└── romannumerals.py      # Roman numeral conversion (113 lines)
```

## Data Sources (`datasources.py`)

### Gene Models
- **`Gff3Annotations`** - Loads GFF3 files (Ensembl format) with tabix indexing for efficient random access
  - Filters by source, gene type, transcript type, component type
  - Supports chromosome name conversion (Ensembl ↔ UCSC ↔ yeast)
  - Tracks hierarchy: genes → transcripts → components (exons/UTRs/CDS)

### Vector/Signal Data
- **`TagDirectory`** - Queries HOMER tag directories with normalization
- **`SparseVectors`** - Memory-efficient storage for sparse genomic vectors

### Interval Data
- **`IntervalData`** - Flexible loader supporting:
  - BED files
  - HOMER peak files
  - HOMER annotated peak files
  - Pandas DataFrames
  - Auto-detects and normalizes column names

### Matrix Data
- **`HicDataDir`** - Loads Hi-C data from per-chromosome files
- **`HicDataDict`** - In-memory Hi-C data from dictionary of DataFrames

## Plot Classes (`plots.py`)

### Core Classes
- **`GenomeBrowser`** - Master orchestrator that renders multi-track visualizations
- **`_BrowserSubPlot`** - Abstract base class for all subplot types

### Visualization Types
- **`GeneModelPlot`** - Transcripts with exons, UTRs, CDS, and strand-specific chevrons
  - Parameters: feature heights, chevron styling, truncation visualization
  - `genes_include` parameter to filter specific genes
  - Automatically arranges overlapping genes at different vertical levels

- **`WigPlot`** - Signal/coverage tracks
  - Features: solid fill vs line plots, alpha transparency, smoothing, auto-scaling

- **`BedPlot`** - Interval/peak data
  - Features: color mapping, value display, custom patch styling

- **`InteractionPlot`** - Long-range genomic interactions using arcs

- **`HicPlot`** - Hi-C matrix visualization with rotated display
  - Features: transform functions, colormap support, diagonal masking

### Utility Functions
- Arc drawing functions (`draw_arc()`, `draw_visible_arc()`, `draw_arc_interaction()`)
- `compute_ax_row_positions()` - Calculates subplot layout geometry
- `match_ylims()` - Synchronizes y-axis limits across subplots

## Utilities (`utilities.py`)

- `roundto()` - Rounds numbers to nearest increment
- `convert_chromosome_name()` - Converts between naming conventions (ucsc, ensembl, yeast)
- `binary_search_tag_file()` - Efficient binary search in tag files
- `bgzip_gff()` & `generate_tabix_index()` - GFF3 file preprocessing
- `gaussian_kernel()` - Generates smoothing kernels
- `add_label()` & `adjust_limits()` - Matplotlib axes management

## Typical Workflow

```python
import pygbrowse

# 1. Load data sources
genes = pygbrowse.datasources.Gff3Annotations(gff3_file)
signal = pygbrowse.datasources.TagDirectory(tag_dir)
peaks = pygbrowse.datasources.IntervalData(bed_file)

# 2. Create plot objects
gene_plot = pygbrowse.plots.GeneModelPlot(genes)
wig_plot = pygbrowse.plots.WigPlot(signal)
bed_plot = pygbrowse.plots.BedPlot(peaks)

# 3. Create browser with nested plot structure
# Each inner list becomes a row; multiple plots in a row share the axis
browser = pygbrowse.plots.GenomeBrowser([
    [gene_plot],
    [wig_plot],
    [bed_plot]
])

# 4. Visualize a genomic region
fig = browser.visualize(
    chrom='chr1',
    start=1000000,
    end=2000000,
    fig_width=12,
    row_heights=[3, 2, 2]
)
```

## Key Design Patterns

1. **Tabix-based indexing** - GFF3 files indexed for efficient random access without loading entire genome
2. **Streaming data access** - Tag directories return only requested genomic regions
3. **Flexible input formats** - Auto-detection and normalization of various file formats
4. **Composition over inheritance** - Data sources and plot types are independent, mixable components
5. **Automatic geometry management** - Intelligently positions overlapping features

## Dependencies

- numpy - Numerical arrays
- scipy - Signal processing, statistics
- pandas - DataFrames for tabular data
- matplotlib - Base plotting library
- seaborn - Plot styling
- pysam - Interface to tabix-indexed files
- intervaltree - Efficient interval overlap queries

## Important Files

- `pygbrowse_demonstration.ipynb` - Demonstration notebook showing the package in action
- `TODO.md` - Planned enhancements (GTF support, BAM files, interactive browsing, etc.)
