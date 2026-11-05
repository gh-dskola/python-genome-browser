"""Helpers: chromosome-name conversion, HOMER tag-file search, GFF3 bgzip/tabix, axes tweaks."""

import gzip
import logging
import os
import tempfile
from pathlib import Path

import matplotlib.axes
import numpy
import pysam
from scipy.stats import norm

logger = logging.getLogger(__name__)


def roundto(num: float, nearest: float) -> int:
    """Round ``num`` to the nearest multiple of ``nearest`` (halves round up).

    Parameters
    ----------
    num : float
        Value to round.
    nearest : float
        Increment to round to.

    Returns
    -------
    int
        The rounded value.
    """
    return int((num + (nearest / 2)) // nearest * nearest)


YEAST_ROMAN_NUMERALS = {
    str(number): numeral
    for number, numeral in enumerate(
        [
            "I",
            "II",
            "III",
            "IV",
            "V",
            "VI",
            "VII",
            "VIII",
            "IX",
            "X",
            "XI",
            "XII",
            "XIII",
            "XIV",
            "XV",
            "XVI",
        ],
        start=1,
    )
}
UCSC_MITOCHONDRIAL_ALIASES = {"M", "MT", "dmel_mitochondrion_genome"}


def convert_chromosome_name(chrom_string: str, dialect: str = "ucsc") -> str:
    """Convert a chromosome name to the given naming dialect.

    Parameters
    ----------
    chrom_string : str
        Chromosome name in any dialect, e.g. ``chr1``, ``1``, ``chrX``, ``MT``, ``XI``.
    dialect : str
        ``ucsc`` (``chr1``, ``chrX``, ``chrM``), ``ensembl`` (``1``, ``X``, ``MT``) or
        ``yeast`` (Roman numerals ``I`` to ``XVI``). Roman numerals are only produced for
        ``yeast``; other dialects keep ``X``, ``V`` etc. as letters.

    Returns
    -------
    str
        The chromosome name in the requested dialect.

    Raises
    ------
    ValueError
        If ``dialect`` is not one of the three above.
    """
    bare_name = chrom_string[3:] if chrom_string[:3].lower() == "chr" else chrom_string
    if dialect == "ucsc":
        return "chrM" if bare_name in UCSC_MITOCHONDRIAL_ALIASES else f"chr{bare_name}"
    if dialect == "ensembl":
        return "MT" if bare_name == "M" else bare_name
    if dialect == "yeast":
        return YEAST_ROMAN_NUMERALS.get(bare_name, bare_name)
    raise ValueError(f"Unknown dialect {dialect!r}; expected 'ucsc', 'ensembl' or 'yeast'")


def binary_search_tag_file(tag_fpath: str | Path, search_target: int) -> int | None:
    """Find the first tag in a sorted HOMER tag file whose start is >= ``search_target``.

    Parameters
    ----------
    tag_fpath : str or Path
        HOMER ``<chrom>.tags.tsv`` file, sorted by start. After stripping, field 2 holds
        the 1-based start; blank lines are skipped.
    search_target : int
        1-based position to search for.

    Returns
    -------
    int or None
        Byte offset of the start of that tag's line, or None if every tag starts before
        ``search_target``.
    """
    with open(tag_fpath, "rb") as tag_file:

        def first_line_from(offset: int) -> tuple[int, bytes]:
            """Return the offset and content of the first non-blank line at or after offset."""
            tag_file.seek(max(offset - 1, 0))
            if offset > 0:
                tag_file.readline()
            while True:
                line_offset, line = tag_file.tell(), tag_file.readline()
                if not line or line.strip():
                    return line_offset, line

        low, high = 0, os.path.getsize(tag_fpath)
        while low < high:
            middle = (low + high) // 2
            _, line = first_line_from(middle)
            if not line or int(line.split()[1]) >= search_target:
                high = middle
            else:
                low = middle + 1
        line_offset, line = first_line_from(low)
        return line_offset if line else None


def _gff3_sort_key(line: str) -> tuple[str, int]:
    """Sort GFF3 rows like ``sort -k1,1 -k4,4n``: contig, then numeric start."""
    fields = line.split("\t", 4)
    return fields[0], int(fields[3])


def bgzip_gff(gff3_fpath: str | Path, bgzipped_fpath: str | Path) -> Path:
    """Sort a GFF3 file by contig and start, and write it block-gzipped.

    Parameters
    ----------
    gff3_fpath : str or Path
        Input GFF3; read as gzip if the name ends in ``.gz``, otherwise as plain text.
    bgzipped_fpath : str or Path
        Output path for the sorted, block-gzipped GFF3. Overwritten if it exists.

    Returns
    -------
    Path
        ``bgzipped_fpath``.

    Raises
    ------
    ValueError
        If input and output are the same file.
    FileNotFoundError
        If ``gff3_fpath`` does not exist.
    """
    gff3_fpath, bgzipped_fpath = Path(gff3_fpath), Path(bgzipped_fpath)
    if gff3_fpath.resolve() == bgzipped_fpath.resolve():
        raise ValueError(f"Source and destination are the same file: {gff3_fpath}")
    if not gff3_fpath.is_file():
        raise FileNotFoundError(f"GFF3 file not found: {gff3_fpath}")

    opener = gzip.open if gff3_fpath.suffix == ".gz" else open
    # ponytail: holds the whole GFF3 in memory; stream through `sort` if that ever runs out.
    with opener(gff3_fpath, "rt") as gff3_file:
        lines = gff3_file.readlines()
    header_lines = [line for line in lines if line.startswith("#")]
    feature_lines = sorted(
        (line for line in lines if line.strip() and not line.startswith("#")), key=_gff3_sort_key
    )
    logger.info("Read %d feature lines from %s", len(feature_lines), gff3_fpath)

    with tempfile.TemporaryDirectory(dir=bgzipped_fpath.parent) as scratch_path:
        sorted_fpath = Path(scratch_path) / "sorted.gff3"
        sorted_fpath.write_text("".join(header_lines + feature_lines))
        pysam.tabix_compress(str(sorted_fpath), str(bgzipped_fpath), force=True)
    logger.info("Wrote sorted, block-gzipped GFF3 to %s", bgzipped_fpath)
    return bgzipped_fpath


def generate_tabix_index(target_fpath: str | Path) -> Path:
    """Build a tabix index (``.tbi``) next to a block-gzipped, sorted GFF3 file.

    Parameters
    ----------
    target_fpath : str or Path
        Block-gzipped GFF3, e.g. the output of ``bgzip_gff``. Its directory must be writeable.

    Returns
    -------
    Path
        The indexed file.
    """
    indexed_fpath = Path(pysam.tabix_index(str(target_fpath), preset="gff", force=True))
    logger.info("Indexed %s", indexed_fpath)
    return indexed_fpath


def gaussian_kernel(sd: float, sd_cutoff: float = 3, normalize: bool = False) -> numpy.ndarray:
    """Return a Gaussian kernel for smoothing per-base signal.

    Parameters
    ----------
    sd : float
        Standard deviation in bp.
    sd_cutoff : float
        Half-width of the kernel in standard deviations.
    normalize : bool
        Scale so the peak is 1; otherwise the values are the normal PDF.

    Returns
    -------
    numpy.ndarray
        ``2 * int(sd_cutoff * sd) + 1`` kernel weights, centered on the middle element.
    """
    half_width = int(sd_cutoff * sd)
    kernel = norm(scale=sd).pdf(numpy.arange(-half_width, half_width + 1))
    return kernel / kernel.max() if normalize else kernel


def add_label(ax: matplotlib.axes.Axes, tick: float, tick_label: str, axis: str = "x") -> None:
    """Add a labeled tick to an axis, or relabel the tick if it exists.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to change.
    tick : float
        Tick position in data units.
    tick_label : str
        Label for the tick.
    axis : str
        ``x`` or ``y``.
    """
    if axis == "y":
        tick_getter, label_getter = ax.get_yticks, ax.get_yticklabels
        tick_setter, label_setter = ax.set_yticks, ax.set_yticklabels
    else:
        tick_getter, label_getter = ax.get_xticks, ax.get_xticklabels
        tick_setter, label_setter = ax.set_xticks, ax.set_xticklabels

    labels = dict(zip(tick_getter(), label_getter()))
    labels[tick] = tick_label
    new_ticks, new_labels = zip(*sorted(labels.items()))
    tick_setter(new_ticks)
    label_setter(new_labels)


def adjust_limits(
    ax: matplotlib.axes.Axes,
    new_position: float,
    axis: str = "y",
    padding_fraction: float = 0.1,
) -> None:
    """Widen an axis' limits, if needed, so that ``new_position`` is inside the content range.

    The current limits are assumed to be the content range padded on each side by
    ``padding_fraction`` of its span. If ``new_position`` lies outside that content range,
    the range is extended to it and padded again by the same fraction.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to change.
    new_position : float
        Position, in data units, that must fit.
    axis : str
        ``x`` or ``y``.
    padding_fraction : float
        Padding on each side as a fraction of the content span; must be below 0.5.

    Raises
    ------
    ValueError
        If ``padding_fraction`` is 0.5 or more.
    """
    if padding_fraction >= 0.5:
        raise ValueError(f"padding_fraction must be below 0.5, got {padding_fraction}")
    limit_getter, limit_setter = (
        (ax.get_ylim, ax.set_ylim) if axis == "y" else (ax.get_xlim, ax.set_xlim)
    )

    current_plot_min, current_plot_max = limit_getter()
    current_data_span = (current_plot_max - current_plot_min) / (1 + 2 * padding_fraction)
    current_pad = current_data_span * padding_fraction
    current_data_min = current_plot_min + current_pad
    current_data_max = current_plot_max - current_pad

    if current_data_min <= new_position <= current_data_max:
        return
    new_data_min = min(current_data_min, new_position)
    new_data_max = max(current_data_max, new_position)
    new_pad = (new_data_max - new_data_min) * padding_fraction
    limit_setter((new_data_min - new_pad, new_data_max + new_pad))
