"""Tests for pygbrowse.utilities and the figure row layout."""

import pytest

from pygbrowse import plots, utilities


@pytest.mark.parametrize(
    ("chrom", "dialect", "expected"),
    [
        ("X", "ucsc", "chrX"),
        ("chrX", "ensembl", "X"),
        ("V", "ucsc", "chrV"),
        ("XI", "ucsc", "chrXI"),
        ("1", "ucsc", "chr1"),
        ("chr1", "ucsc", "chr1"),
        ("chr10", "ensembl", "10"),
        ("chrM", "ensembl", "MT"),
        ("MT", "ucsc", "chrM"),
        ("chrXI", "yeast", "XI"),
        ("chr11", "yeast", "XI"),
        ("11", "yeast", "XI"),
        ("XVI", "yeast", "XVI"),
        ("chr16", "yeast", "XVI"),
        ("chrM", "yeast", "M"),
    ],
)
def test_convert_chromosome_name(chrom: str, dialect: str, expected: str) -> None:
    """Names convert between dialects without Roman-numeral misreads outside yeast."""
    assert utilities.convert_chromosome_name(chrom, dialect=dialect) == expected


def test_convert_chromosome_name_rejects_unknown_dialect() -> None:
    """An unknown dialect raises ValueError."""
    with pytest.raises(ValueError, match="dialect"):
        utilities.convert_chromosome_name("chr1", dialect="refseq")


@pytest.mark.parametrize("num_rows", range(1, 21))
@pytest.mark.parametrize("ax_spacing", [0.0, 0.05, 0.1])
def test_row_positions_fit_in_figure(num_rows: int, ax_spacing: float) -> None:
    """Every row lies inside [0, 1] in figure coordinates and rows do not overlap."""
    bottoms, heights = plots.compute_ax_row_positions(
        row_heights=[1] * num_rows, ax_spacing=ax_spacing
    )
    tops = [bottom + height for bottom, height in zip(bottoms, heights)]
    assert min(bottoms) >= -1e-9
    assert max(tops) <= 1 + 1e-9
    for upper_bottom, lower_top in zip(bottoms, tops[1:]):
        assert lower_top <= upper_bottom + 1e-9


@pytest.mark.parametrize(("num", "nearest", "expected"), [(14, 10, 10), (15, 10, 20), (0, 5, 0)])
def test_roundto(num: float, nearest: float, expected: int) -> None:
    """Rounds to the nearest multiple; halves round up."""
    assert utilities.roundto(num, nearest) == expected


def test_gaussian_kernel_is_symmetric_and_peaks_in_the_middle() -> None:
    """The kernel has 2 * 3 * sd + 1 weights, peaks at the center, and normalizes to 1."""
    kernel = utilities.gaussian_kernel(sd=2, normalize=True)
    assert len(kernel) == 13
    assert kernel.argmax() == 6
    assert kernel.max() == 1.0
    assert list(kernel) == list(kernel[::-1])
