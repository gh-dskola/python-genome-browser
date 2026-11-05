"""Tests for pygbrowse.datasources. All coordinates asserted here are 0-based half-open."""

import signal
from collections.abc import Iterator
from pathlib import Path

import numpy
import pandas
import pytest
from conftest import BIGWIG_SIGNAL_END, BED_ROWS, write_gff3

from pygbrowse import datasources, utilities


@pytest.fixture
def alarm() -> Iterator[None]:
    """Fail a test that runs longer than 5 seconds instead of hanging the suite."""

    def on_timeout(signum: int, frame: object) -> None:
        raise TimeoutError("test exceeded 5 s; probable infinite loop")

    previous_handler = signal.signal(signal.SIGALRM, on_timeout)
    signal.alarm(5)
    yield
    signal.alarm(0)
    signal.signal(signal.SIGALRM, previous_handler)


# --- TagDirectory -----------------------------------------------------------------------


def test_tag_directory_per_base_coverage(tag_directory_path: Path, alarm: None) -> None:
    """A read of length L at 1-based start S covers exactly 0-based [S-1, S-1+L)."""
    tags = datasources.TagDirectory(tag_directory_path, 4)
    coverage = tags.query("chr1", 90, 130)
    expected = pandas.Series(numpy.zeros(40), index=numpy.arange(90, 130))
    expected.loc[100:109] += 1.0
    expected.loc[105:124] += 2.0
    pandas.testing.assert_series_equal(coverage, expected)


def test_tag_directory_without_normalization(tag_directory_path: Path, alarm: None) -> None:
    """normalize_to=None returns raw depths."""
    tags = datasources.TagDirectory(tag_directory_path, None)
    assert tags.query("chr1", 100, 101).iloc[0] == 1.0


def test_tag_directory_normalizes_to_target(tag_directory_path: Path, alarm: None) -> None:
    """Depths scale by normalize_to / total tags."""
    tags = datasources.TagDirectory(tag_directory_path, 8)
    assert tags.query("chr1", 100, 101).iloc[0] == 2.0


def test_tag_directory_query_reaching_last_tag(tag_directory_path: Path, alarm: None) -> None:
    """A query that covers the file's last tag stops at end of file."""
    coverage = datasources.TagDirectory(tag_directory_path, 4).query("chr1", 400, 600)
    assert coverage.sum() == 10.0


def test_tag_directory_query_past_last_tag(tag_directory_path: Path, alarm: None) -> None:
    """A query downstream of every tag returns zeros."""
    coverage = datasources.TagDirectory(tag_directory_path, 4).query("chr1", 1000, 1100)
    assert len(coverage) == 100
    assert coverage.sum() == 0.0


def test_tag_directory_read_starts(tag_directory_path: Path, alarm: None) -> None:
    """read_handling='starts' counts each read once, at its 0-based start."""
    tags = datasources.TagDirectory(tag_directory_path, 4, read_handling="starts")
    coverage = tags.query("chr1", 100, 130)
    assert coverage.loc[100] == 1.0
    assert coverage.loc[105] == 2.0
    assert coverage.sum() == 3.0


@pytest.mark.parametrize(
    ("target", "expected_line"), [(50, 0), (101, 0), (102, 1), (106, 1), (107, 2), (501, 2)]
)
def test_binary_search_tag_file(
    tag_directory_path: Path, target: int, expected_line: int
) -> None:
    """Returns the byte offset of the first tag whose 1-based start is >= target."""
    tag_fpath = tag_directory_path / "chr1.tags.tsv"
    line_offsets = numpy.cumsum([0] + [len(line) for line in tag_fpath.read_bytes().splitlines(True)])
    assert utilities.binary_search_tag_file(tag_fpath, target) == line_offsets[expected_line]


def test_binary_search_tag_file_past_end(tag_directory_path: Path) -> None:
    """Returns None when every tag starts before the target."""
    assert utilities.binary_search_tag_file(tag_directory_path / "chr1.tags.tsv", 502) is None


# --- SparseVectors ----------------------------------------------------------------------


def test_sparse_vectors_query_from_before_first_point() -> None:
    """A query starting before the first point keeps the leading points."""
    series = pandas.Series([1.0, 2.0, 3.0], index=[10, 20, 30])
    result = datasources.SparseVectors({"chr1": series}).query("chr1", 0, 25)
    assert list(result.index[:2]) == [10, 20]


# --- IntervalData -----------------------------------------------------------------------


def test_interval_data_reads_bed(bed_fpath: Path) -> None:
    """A six-column BED file loads with BED column names and names as the index."""
    intervals = datasources.IntervalData(bed_fpath)
    assert len(intervals.data) == len(BED_ROWS)
    assert intervals.data.loc["peak1", "chromStart"] == 100
    assert intervals.data.loc["peak1", "chromEnd"] == 200
    assert intervals.data.loc["peak2", "score"] == -3.0
    assert intervals.data.loc["peak2", "strand"] == "-"


def test_interval_data_named_frame_without_strand_keeps_columns() -> None:
    """A DataFrame with BED column names but no strand column is not relabelled."""
    frame = pandas.DataFrame(
        {"chrom": ["chr1"], "chromStart": [10], "chromEnd": [20], "name": ["a"], "score": [7.0]}
    )
    intervals = datasources.IntervalData(frame)
    assert intervals.data["score"].tolist() == [7.0]
    assert intervals.data["name"].tolist() == ["a"]


def test_interval_data_unknown_format_raises(bed_fpath: Path) -> None:
    """An unknown file format fails at construction."""
    with pytest.raises(ValueError, match="format"):
        datasources.IntervalData(bed_fpath, "gff")


# --- Gff3Annotations --------------------------------------------------------------------


def test_gff3_hierarchy(gff3_fpath: Path) -> None:
    """Genes, transcripts and components link up, with 0-based half-open coordinates."""
    genes, transcripts, components, names_to_ids = datasources.Gff3Annotations(
        gff3_fpath
    ).query("chr1", 0, 1000)
    assert set(genes) == {"gene:GA", "gene:GB"}
    assert (genes["gene:GA"]["start"], genes["gene:GA"]["end"]) == (100, 600)
    assert genes["gene:GA"]["transcripts"] == ["transcript:TA"]
    assert len(transcripts["transcript:TA"]["components"]) == 4
    assert names_to_ids["GENEB"] == ["gene:GB"]
    component_starts = sorted(component["start"] for component in components.values())
    assert component_starts[0] == 100


def test_gff3_attribute_parsing(tmp_path: Path) -> None:
    """Trailing ';', '=' inside values, %-escapes, comma-separated Parent and missing Name."""
    gff3_fpath = write_gff3(
        tmp_path / "odd.gff3",
        [
            ("gene", 101, 600, "+", "ID=gene:GD;Name=GENE%3BD;note=a=b;"),
            ("gene", 201, 500, "+", "ID=gene:GE"),
            ("mRNA", 101, 600, "+", "ID=transcript:TD1;Parent=gene:GD"),
            ("mRNA", 101, 600, "+", "ID=transcript:TD2;Parent=gene:GD"),
            ("CDS", 151, 250, "+", "Parent=transcript:TD1,transcript:TD2;"),
        ],
    )
    genes, transcripts, _, names_to_ids = datasources.Gff3Annotations(gff3_fpath).query(
        "chr1", 0, 1000
    )
    assert genes["gene:GD"]["Name"] == "GENE;D"
    assert genes["gene:GD"]["note"] == "a=b"
    assert genes["gene:GE"]["Name"] == "gene:GE"
    assert names_to_ids["GENE;D"] == ["gene:GD"]
    assert len(transcripts["transcript:TD1"]["components"]) == 1
    assert len(transcripts["transcript:TD2"]["components"]) == 1


# --- GFF3 preparation helpers -----------------------------------------------------------


def test_bgzip_and_index_unsorted_gff3(tmp_path: Path) -> None:
    """An unsorted plain GFF3 is sorted, bgzipped and indexed so that it can be queried."""
    plain_fpath = tmp_path / "unsorted.gff3"
    plain_fpath.write_text(
        "##gff-version 3\n"
        "1\tensembl\tgene\t2001\t2500\t.\t+\t.\tID=gene:GC;Name=GENEC\n"
        "1\tensembl\tgene\t101\t600\t.\t+\t.\tID=gene:GA;Name=GENEA\n"
    )
    bgzipped_fpath = tmp_path / "sorted.gff3.gz"
    utilities.bgzip_gff(plain_fpath, bgzipped_fpath)
    utilities.generate_tabix_index(bgzipped_fpath)
    genes, _, _, _ = datasources.Gff3Annotations(bgzipped_fpath).query("chr1", 0, 3000)
    assert set(genes) == {"gene:GA", "gene:GC"}


def test_bgzip_gff_missing_input_raises(tmp_path: Path) -> None:
    """A missing input file raises instead of logging and returning."""
    with pytest.raises(FileNotFoundError):
        utilities.bgzip_gff(tmp_path / "absent.gff3", tmp_path / "out.gff3.gz")


def test_bgzip_gff_same_source_and_destination_raises(tmp_path: Path) -> None:
    """Writing the output over the input is refused."""
    plain_fpath = tmp_path / "a.gff3"
    plain_fpath.write_text("##gff-version 3\n")
    with pytest.raises(ValueError, match="same"):
        utilities.bgzip_gff(plain_fpath, plain_fpath)


# --- Hi-C -------------------------------------------------------------------------------


def test_hic_data_dir_query(tmp_path: Path) -> None:
    """HicDataDir reads one per-chromosome matrix and slices it to the rounded window."""
    bins = list(range(0, 100, 10))
    matrix = pandas.DataFrame(numpy.ones((10, 10)), index=bins, columns=bins)
    matrix.to_csv(tmp_path / "hic_chr1.tsv", sep="\t")
    hic = datasources.HicDataDir(str(tmp_path / "hic_{}.tsv"), 10)
    assert hic.bin_size == 10
    assert hic.query("chr1", 20, 50).shape == (4, 4)


# --- BigWigDataSource -------------------------------------------------------------------


def test_bigwig_query_values_and_padding(bigwig_fpath: Path) -> None:
    """Values are 0-based; positions without data and past the chromosome end are filled."""
    bigwig = datasources.BigWigDataSource(bigwig_fpath)
    result = bigwig.query("chr1", 995, 1005)
    assert list(result.index) == list(range(995, 1005))
    assert result.loc[995] == pytest.approx(99.5)
    assert result.loc[BIGWIG_SIGNAL_END] == 0.0
    padded = bigwig.query("chr2", 4990, 5010)
    assert len(padded) == 20
    assert padded.sum() == 0.0
