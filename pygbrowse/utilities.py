import datetime
import gzip
import logging
import os
import tempfile
from pathlib import Path

import numpy
import pysam
from scipy.stats import norm

logger = logging.getLogger(__name__)


# ToDo: Bring back scale bar
# ToDo: Add option for solid fill of vectors


def roundto(num, nearest):
    """
    Rounds :param:`num` to the nearest increment of :param:`nearest`
    """
    return int((num + (nearest / 2)) // nearest * nearest)


YEAST_ROMAN_NUMERALS = {
    str(number): numeral
    for number, numeral in enumerate(
        ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV",
         "XV", "XVI"],
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
        HOMER ``<chrom>.tags.tsv`` file, sorted by start. Column 2 holds the 1-based start.
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
            """Return the offset and content of the first line starting at or after offset."""
            tag_file.seek(max(offset - 1, 0))
            if offset > 0:
                tag_file.readline()
            return tag_file.tell(), tag_file.readline()

        low, high = 0, os.path.getsize(tag_fpath)
        while low < high:
            middle = (low + high) // 2
            _, line = first_line_from(middle)
            if not line or int(line.split(b"\t")[1]) >= search_target:
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


def pretty_now():
    """
    Returns the current date/time in a nicely formatted string (without decimal seconds)
    """
    return datetime.datetime.strftime(datetime.datetime.now(), "%Y-%b-%d %H:%M:%S")


def log_print(message, tabs=1):
    """
    Print a chunk of text preceded by a timestamp and an optional number of tabs (default 1).

    :param message:
    :param tabs:
    :return:
    """
    print("{}{}{}".format(pretty_now(), "\t" * tabs, message))


def gaussian_kernel(sd, sd_cutoff=3, normalize=False):
    """
    Generate and return a numpy.Array whose elements are proportional to the PDF of a normal distribution
    having standard deviation :param:`sd`.

    :param sd:
    :param sd_cutoff:
    :param normalize:
    :return:
    """
    bw = sd_cutoff * sd * 2 + 1
    midpoint = sd_cutoff * sd
    kern = numpy.zeros(bw)
    frozen_rv = norm(scale=sd)
    for i in range(bw):
        kern[i] = frozen_rv.pdf(i - midpoint)
    if normalize:
        kern = kern / kern.max()
    return kern


def add_label(ax, tick, tick_label, axis="x"):
    """
    Updates the set of ticks and tick labels for the specified matplotlib.Axes object
    and axis.

    If the tick already exists, it's label will be updated. If not, it will be created and labeled
    appropriately.

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


def adjust_limits(ax, new_position, axis="y", padding_fraction=0.1):
    """
    If necessary adjusts the limits for the specified :param axis: on
    :param ax: to accomodate :param new_position: according to the
    following scheme:

    1. Assumes that the current limits are the
        smallest and largest content item minus / plus a padding equal to
        :param padding_fraction: * the span between the smallest
        and largest content item.
    2. If :param new_position: is beyond the inferred content limits,
        adjust the padding to :param padding_fraction: * the new content
        span, then adjust the plot limits to the new content limits
        minus / plus the new padding.
    """
    assert padding_fraction < 0.5, "padding_fraction must be below 0.5!"

    if axis == "y":
        limit_getter = ax.get_ylim
        limit_setter = ax.set_ylim
    else:
        limit_getter = ax.get_xlim
        limit_setter = ax.set_xlim

    current_plot_min, current_plot_max = limit_getter()
    current_plot_span = current_plot_max - current_plot_min
    current_data_span = current_plot_span / (1 + 2 * padding_fraction)
    current_pad = current_data_span * padding_fraction
    current_data_min = current_plot_min + current_pad
    current_data_max = current_plot_max - current_pad

    #     print(current_plot_min, current_plot_max, current_plot_span)
    #     print(current_data_min, current_data_max, current_data_span, current_pad)

    if new_position > current_data_max:
        new_data_min = current_data_min
        new_data_max = new_position

    elif new_position < current_data_min:
        new_data_min = new_position
        new_data_max = current_data_max
    else:
        # no changes needed
        return

    new_data_span = new_data_max - new_data_min
    new_pad = new_data_span * padding_fraction
    new_plot_min = new_data_min - new_pad
    new_plot_max = new_data_max + new_pad

    #     print(new_data_min, new_data_max, new_data_span, new_pad)
    #     print(new_plot_min, new_plot_max)

    limit_setter((new_plot_min, new_plot_max))


def diag_indices(n, k=0):
    """
    Return the indices corresponding to the kth diagonal of an n X n array
    in the form of a tuple of (x coords, y coords).

    Created since numpy does not provide this functionality.
    """
    if k <= 0:
        x_coords = numpy.arange(-k, n)
        y_coords = numpy.arange(0, n + k)
    else:
        x_coords = numpy.arange(0, n - k)
        y_coords = numpy.arange(k, n)

    return (x_coords, y_coords)
