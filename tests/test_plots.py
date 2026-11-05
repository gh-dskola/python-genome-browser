"""Tests for pygbrowse.plots: rendering through GenomeBrowser.visualize on synthetic data."""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy
import pandas
import pytest
from conftest import EXPECTED_COMPONENT_RECTANGLES, GFF3_ROWS, write_gff3

from pygbrowse import datasources, plots

FIG_WIDTH = 12


@pytest.fixture(autouse=True)
def close_figures() -> None:
    """Close every figure after each test."""
    yield
    plt.close("all")


def rectangles(ax: matplotlib.axes.Axes) -> list[matplotlib.patches.Rectangle]:
    """Return the Rectangle patches drawn on an Axes."""
    return [patch for patch in ax.patches if isinstance(patch, matplotlib.patches.Rectangle)]


def render(track: plots._BrowserSubPlot, chrom: str, start: int, end: int) -> plt.Axes:
    """Render a single track and return its Axes."""
    figure = plots.GenomeBrowser([[track]]).visualize(chrom, start, end, fig_width=FIG_WIDTH)
    return figure.axes[0]


def interval_frame() -> pandas.DataFrame:
    """Named-column interval DataFrame with one interval spanning chr1 [0, 5000)."""
    return pandas.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "chromStart": [100, 0],
            "chromEnd": [200, 5000],
            "strand": ["+", "+"],
            "score": [1.0, 2.0],
        },
        index=pandas.Index(["peak", "wide"], name="name"),
    )


# --- GeneModelPlot ----------------------------------------------------------------------


def test_gene_components_drawn_once_per_gene(gff3_fpath: Path) -> None:
    """Each gene draws only its own components, so overlapping genes add no extra boxes."""
    ax = render(plots.GeneModelPlot(datasources.Gff3Annotations(gff3_fpath)), "chr1", 0, 1000)
    assert (
        len(rectangles(ax))
        == EXPECTED_COMPONENT_RECTANGLES["GENEA"] + (EXPECTED_COMPONENT_RECTANGLES["GENEB"])
    )


def test_gene_component_rectangles_stay_inside_their_component(tmp_path: Path) -> None:
    """No UTR or CDS box of a lone + strand gene extends past the component it draws."""
    genea_rows = [row for row in GFF3_ROWS if row[4].endswith(("GA", "TA", "GENEA"))]
    gff3_fpath = write_gff3(tmp_path / "genea.gff3", genea_rows)
    ax = render(plots.GeneModelPlot(datasources.Gff3Annotations(gff3_fpath)), "chr1", 0, 1000)
    # GFF 1-based closed [start, end] -> 0-based half-open [start - 1, end).
    component_bounds = [
        (start - 1, end)
        for feature_type, start, end, _, _ in genea_rows
        if feature_type in ("CDS", "five_prime_UTR", "three_prime_UTR")
    ]
    assert len(rectangles(ax)) == EXPECTED_COMPONENT_RECTANGLES["GENEA"]
    for rectangle in rectangles(ax):
        left, right = rectangle.get_x(), rectangle.get_x() + rectangle.get_width()
        assert any(
            bound_start - 1e-6 <= left and right <= bound_end + 1e-6
            for bound_start, bound_end in component_bounds
        ), f"rectangle [{left}, {right}) is outside every component"


def test_genes_include_filters_before_layout(gff3_fpath: Path) -> None:
    """Only the included gene is drawn and it takes a single row."""
    track = plots.GeneModelPlot(datasources.Gff3Annotations(gff3_fpath), genes_include=["GENEB"])
    ax = render(track, "chr1", 0, 1000)
    assert [text.get_text() for text in ax.texts] == ["GENEB"]
    assert ax.get_ylim() == (-0.5, 0.5)


def test_genes_include_is_exact_match(gff3_fpath: Path) -> None:
    """genes_include matches whole names, not substrings."""
    track = plots.GeneModelPlot(datasources.Gff3Annotations(gff3_fpath), genes_include=["GENEAB"])
    ax = render(track, "chr1", 0, 1000)
    assert list(ax.texts) == []


# --- BedPlot ----------------------------------------------------------------------------


def test_bed_plot_draws_interval_spanning_window() -> None:
    """An interval that covers the whole window is drawn."""
    ax = render(plots.BedPlot(datasources.IntervalData(interval_frame())), "chr1", 1000, 2000)
    assert len(rectangles(ax)) == 1


def test_bed_plot_patch_kwargs_do_not_leak() -> None:
    """Custom patch_kwargs on one BedPlot do not change the class default or other plots."""
    intervals = datasources.IntervalData(interval_frame())
    plots.BedPlot(intervals, patch_kwargs={"linewidth": 3})
    assert plots.BedPlot.DEFAULT_PATCH_KWARGS["linewidth"] == 1
    assert plots.BedPlot(intervals).patch_kwargs["linewidth"] == 1


# --- WigPlot ----------------------------------------------------------------------------


def test_wig_plot_does_not_mutate_source_vector() -> None:
    """Scaling and centering work on a copy of the queried vector."""
    series = pandas.Series(numpy.arange(100, dtype=float), index=numpy.arange(100))
    original = series.copy()
    track = plots.WigPlot(
        datasources.SparseVectors({"chr1": series}), scale_vector_to_plot=True, center_vector=True
    )
    render(track, "chr1", 0, 100)
    pandas.testing.assert_series_equal(series, original)


# --- HicPlot ----------------------------------------------------------------------------


def hic_matrix() -> pandas.DataFrame:
    """A 40 x 40 Hi-C matrix of ones over bins 0, 10, ..., 390."""
    bins = list(range(0, 400, 10))
    return pandas.DataFrame(numpy.ones((40, 40)), index=bins, columns=bins)


def test_hic_plot_leaves_source_matrix_unchanged() -> None:
    """Masking the diagonal does not write into the data source's matrix."""
    matrix = hic_matrix()
    track = plots.HicPlot(datasources.HicDataDict({"chr1": matrix}, 10))
    render(track, "chr1", 100, 300)
    assert (matrix.to_numpy() == 1).all()


def test_hic_plot_under_copy_on_write() -> None:
    """HicPlot renders with pandas Copy-on-Write enabled."""
    with pandas.option_context("mode.copy_on_write", True):
        track = plots.HicPlot(datasources.HicDataDict({"chr1": hic_matrix()}, 10))
        render(track, "chr1", 100, 300)


# --- GenomeBrowser ----------------------------------------------------------------------


def test_visualize_returns_a_new_figure_each_call() -> None:
    """Two calls return two figures."""
    browser = plots.GenomeBrowser([[plots.BedPlot(datasources.IntervalData(interval_frame()))]])
    first = browser.visualize("chr1", 0, 1000)
    second = browser.visualize("chr1", 0, 1000)
    assert first is not second


def test_visualize_leaves_global_rcparams_unchanged() -> None:
    """The track style applies to the figure only, not to matplotlib's global settings."""
    before = dict(matplotlib.rcParams)
    browser = plots.GenomeBrowser([[plots.BedPlot(datasources.IntervalData(interval_frame()))]])
    browser.visualize("chr1", 0, 1000)
    assert dict(matplotlib.rcParams) == before
