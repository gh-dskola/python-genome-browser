"""Shared fixtures: synthetic GFF3, BigWig, HOMER tag directory and BED files built in tmp_path.

Coordinates in the fixture docstrings are given in each file format's native system;
the 0-based half-open equivalents the data sources must return are stated alongside.
"""

import matplotlib

matplotlib.use("Agg")

from pathlib import Path  # noqa: E402

import numpy  # noqa: E402
import pyBigWig  # noqa: E402
import pysam  # noqa: E402
import pytest  # noqa: E402

GFF3_HEADER = "##gff-version 3\n"

# Contig "1" (Ensembl dialect; queried as "chr1").
# GENEA + strand, GFF 101..600 -> 0-based [100, 600).
# GENEB - strand, GFF 401..900 -> 0-based [400, 900); overlaps GENEA.
# GENEC + strand, GFF 2001..2500 -> 0-based [2000, 2500); isolated.
GFF3_ROWS = [
    ("gene", 101, 600, "+", "ID=gene:GA;Name=GENEA"),
    ("mRNA", 101, 600, "+", "ID=transcript:TA;Parent=gene:GA"),
    ("five_prime_UTR", 101, 150, "+", "Parent=transcript:TA"),
    ("CDS", 151, 250, "+", "Parent=transcript:TA"),
    ("CDS", 351, 500, "+", "Parent=transcript:TA"),
    ("gene", 401, 900, "-", "ID=gene:GB;Name=GENEB"),
    ("mRNA", 401, 900, "-", "ID=transcript:TB;Parent=gene:GB"),
    ("three_prime_UTR", 401, 450, "-", "Parent=transcript:TB"),
    ("CDS", 451, 800, "-", "Parent=transcript:TB"),
    ("five_prime_UTR", 801, 900, "-", "Parent=transcript:TB"),
    ("three_prime_UTR", 501, 600, "+", "Parent=transcript:TA"),
    ("gene", 2001, 2500, "+", "ID=gene:GC;Name=GENEC"),
    ("mRNA", 2001, 2500, "+", "ID=transcript:TC;Parent=gene:GC"),
    ("CDS", 2001, 2500, "+", "Parent=transcript:TC"),
]

# Number of matplotlib Rectangles each gene's components draw: CDS -> 1, UTR -> 2.
EXPECTED_COMPONENT_RECTANGLES = {"GENEA": 6, "GENEB": 5, "GENEC": 1}


def write_gff3(gff3_fpath: Path, rows: list[tuple]) -> Path:
    """Write GFF3 rows on contig "1", then bgzip and tabix-index the file.

    Parameters
    ----------
    gff3_fpath : Path
        Destination for the uncompressed GFF3; the bgzipped copy is written next to it.
    rows : list of tuple
        (feature_type, start, end, strand, attributes) with GFF 1-based closed coordinates.

    Returns
    -------
    Path
        The bgzipped, tabix-indexed GFF3 file.
    """
    sorted_rows = sorted(rows, key=lambda row: row[1])
    lines = [
        f"1\tensembl\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n"
        for feature_type, start, end, strand, attributes in sorted_rows
    ]
    gff3_fpath.write_text(GFF3_HEADER + "".join(lines))
    return Path(pysam.tabix_index(str(gff3_fpath), preset="gff", force=True))


@pytest.fixture
def gff3_fpath(tmp_path: Path) -> Path:
    """Bgzipped, tabix-indexed GFF3 with the three genes described in GFF3_ROWS."""
    return write_gff3(tmp_path / "genes.gff3", GFF3_ROWS)


# chr1 is 10,000 bp; signal covers 0-based [0, 1000) with value position / 10.
# chr2 is 5,000 bp with no signal.
BIGWIG_CHROM_SIZES = [("chr1", 10000), ("chr2", 5000)]
BIGWIG_SIGNAL_END = 1000


@pytest.fixture
def bigwig_fpath(tmp_path: Path) -> Path:
    """BigWig with signal value position / 10 on chr1 [0, 1000)."""
    bigwig_fpath = tmp_path / "signal.bw"
    bigwig_file = pyBigWig.open(str(bigwig_fpath), "w")
    bigwig_file.addHeader(BIGWIG_CHROM_SIZES)
    bigwig_file.addEntries(
        "chr1",
        0,
        values=list(numpy.arange(BIGWIG_SIGNAL_END, dtype=float) / 10),
        span=1,
        step=1,
    )
    bigwig_file.close()
    return bigwig_fpath


# HOMER tag rows as (1-based start, strand, depth, read length) on chr1.
# Row 1 covers 0-based [100, 110), row 2 [105, 125), row 3 [500, 510).
TAG_ROWS = [(101, 0, 1.0, 10), (106, 1, 2.0, 20), (501, 0, 1.0, 10)]
TAG_TOTAL = 4


@pytest.fixture
def tag_directory_path(tmp_path: Path) -> Path:
    """HOMER tag directory with chr1.tags.tsv and a tagInfo.txt reporting TAG_TOTAL tags."""
    tag_directory_path = tmp_path / "tags"
    tag_directory_path.mkdir()
    (tag_directory_path / "tagInfo.txt").write_text(
        f"name\tUnique Positions\tTotal Tags\ngenome=test\t{len(TAG_ROWS)}\t{TAG_TOTAL}\n"
    )
    (tag_directory_path / "chr1.tags.tsv").write_text(
        "".join(
            f"chr1\t{start}\t{strand}\t{depth}\t{length}\n"
            for start, strand, depth, length in TAG_ROWS
        )
    )
    return tag_directory_path


# BED rows (0-based half-open): chrom, start, end, name, score, strand.
BED_ROWS = [
    ("chr1", 100, 200, "peak1", 5.0, "+"),
    ("chr1", 300, 400, "peak2", -3.0, "-"),
    ("chr1", 0, 5000, "wide", 1.0, "+"),
    ("chr2", 100, 200, "other", 2.0, "+"),
]


@pytest.fixture
def bed_fpath(tmp_path: Path) -> Path:
    """Six-column BED file with the rows in BED_ROWS."""
    bed_fpath = tmp_path / "peaks.bed"
    bed_fpath.write_text("".join("\t".join(map(str, row)) + "\n" for row in BED_ROWS))
    return bed_fpath
