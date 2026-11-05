"""Data sources that answer region queries for the browser tracks.

Every data source has ``query(chrom, start, end)``. ``chrom`` is a UCSC-style name
(``chr1``); ``start`` and ``end`` are **0-based, half-open**. Sources convert to the
file's own chromosome naming and coordinate system internally:

- GFF3 (1-based, closed) is converted in ``Gff3Annotations``.
- HOMER tag files (1-based starts) are converted in ``TagDirectory``.
- BED and BigWig are already 0-based, half-open.

Exceptions: Hi-C matrix queries (``HicDataDir``, ``HicDataDict``) round both ends to bins
and include both end bins. ``SparseVectors`` also returns the nearest point outside the
window on each side.
"""

import abc
import contextlib
import logging
import os
from collections.abc import Callable, Iterator
from pathlib import Path
from urllib.parse import unquote

import numpy
import pandas
import pyBigWig
import pysam
from scipy.signal import convolve

from . import utilities

logger = logging.getLogger(__name__)

DEFAULT_TAG_COUNT_NORMALIZATION_TARGET = 10000000
DEFAULT_FEATURE_SOURCES = ("ensembl", "havana", "ensembl_havana")
DEFAULT_GENE_TYPES = (
    "gene",
    "RNA",
    "mt_gene",
    "lincRNA_gene",
    "miRNA_gene",
    "ncRNA_gene",
    "rRNA_gene",
    "snRNA_gene",
    "snoRNA_gene",
    "processed_transcript",
)
DEFAULT_TRANSCRIPT_TYPES = (
    "mRNA",
    "transcript",
    "lincRNA",
    "lnc_RNA",
    "miRNA",
    "ncRNA",
    "snRNA",
    "snoRNA",
)
DEFAULT_COMPONENT_TYPES = ("CDS", "three_prime_UTR", "five_prime_UTR")
CHROMOSOME_DIALECTS = ("ucsc", "ensembl", "yeast")

VectorTransform = Callable[[pandas.Series], pandas.Series]


class _ChromWrapper:
    """``source[chrom].loc[start:end]`` indexing for vector data sources."""

    def __init__(self, chrom: str, parent_data_source: "_VectorDataSource"):
        self.chrom = chrom
        self.parent_data_source = parent_data_source

    def __getitem__(self, key: slice | int) -> pandas.Series:
        """Query a slice ``[start:end]`` (0-based half-open) or a single position."""
        if isinstance(key, slice):
            start, end = key.start, key.stop
        else:
            start, end = key, key + 1
        return self.parent_data_source.query(self.chrom, start, end)


class _DataVector:
    def __init__(self, chrom: str, parent_data_source: "_VectorDataSource"):
        self.loc = _ChromWrapper(chrom=chrom, parent_data_source=parent_data_source)


class _VectorDataSource(abc.ABC):
    def __init__(self, transform: VectorTransform | None = None, smoothing_bandwidth: float = 0):
        """Base class for per-base signal sources.

        Parameters
        ----------
        transform : callable, optional
            Applied to every query result, after smoothing.
        smoothing_bandwidth : float
            Standard deviation (bp) of a Gaussian smoothing kernel; 0 disables smoothing.
        """
        self.transform = transform
        if smoothing_bandwidth:
            self.convolution_kernel = utilities.gaussian_kernel(smoothing_bandwidth)
        else:
            self.convolution_kernel = None

    @abc.abstractmethod
    def _query(self, chrom: str, start: int, end: int) -> pandas.Series:
        """Return the raw signal for a 0-based half-open region, indexed by position."""

    def query(self, chrom: str, start: int, end: int) -> pandas.Series:
        """Return the smoothed, transformed signal for a region.

        Parameters
        ----------
        chrom : str
            Chromosome, UCSC style.
        start, end : int
            Region, 0-based half-open.

        Returns
        -------
        pandas.Series
            Signal values indexed by 0-based position.
        """
        query_result = self._query(chrom, start, end)
        if self.convolution_kernel is not None:
            query_result = pandas.Series(
                convolve(query_result, self.convolution_kernel, mode="same"),
                index=query_result.index,
            )
        if self.transform:
            query_result = self.transform(query_result)
        return query_result

    def __getitem__(self, chrom: str) -> _DataVector:
        """Return an object whose ``.loc[start:end]`` queries ``chrom``."""
        return _DataVector(chrom=chrom, parent_data_source=self)


class SparseVectors(_VectorDataSource):
    def __init__(
        self,
        series_dict: dict[str, pandas.Series],
        transform: VectorTransform | None = None,
        smoothing_bandwidth: float = 0,
    ):
        """Signal held in memory as one sparse Series per chromosome.

        Parameters
        ----------
        series_dict : dict of str to pandas.Series
            Per-chromosome values indexed by sorted 0-based position.
        transform : callable, optional
            Applied to every query result.
        smoothing_bandwidth : float
            Gaussian smoothing standard deviation in bp; 0 disables smoothing.
        """
        super().__init__(transform=transform, smoothing_bandwidth=smoothing_bandwidth)
        self.data = series_dict

    def _query(self, chrom, start, end):
        this_chrom_vector = self.data[chrom]
        # One point either side of the window is kept so that line plots reach the edges.
        start_ipos = max(numpy.searchsorted(this_chrom_vector.index, start) - 1, 0)
        end_ipos = numpy.searchsorted(this_chrom_vector.index, end) + 1
        return this_chrom_vector.iloc[start_ipos:end_ipos]


class TagDirectory(_VectorDataSource):
    READ_HANDLING_MODES = ("reads", "starts")

    def __init__(
        self,
        tag_directory_path: str | Path,
        normalize_to: float | None = DEFAULT_TAG_COUNT_NORMALIZATION_TARGET,
        transform: VectorTransform | None = None,
        smoothing_bandwidth: float = 0,
        read_handling: str = "reads",
    ):
        """Read coverage from a HOMER tag directory (``<chrom>.tags.tsv`` + ``tagInfo.txt``).

        Tag lines are stripped and split on whitespace into chrom, 1-based start, strand,
        depth, read length; blank lines are skipped. Strand is ignored.

        Parameters
        ----------
        tag_directory_path : str or Path
            HOMER tag directory.
        normalize_to : float, optional
            Scale depths so that the directory's total tag count equals this value.
            ``None`` or 0 returns raw depths.
        transform : callable, optional
            Applied to every query result.
        smoothing_bandwidth : float
            Gaussian smoothing standard deviation in bp; 0 disables smoothing.
        read_handling : str
            ``reads`` spreads each tag over its read length; ``starts`` counts it once,
            at its start.

        Raises
        ------
        ValueError
            If ``read_handling`` is not one of ``READ_HANDLING_MODES``.
        """
        super().__init__(transform=transform, smoothing_bandwidth=smoothing_bandwidth)
        if read_handling not in self.READ_HANDLING_MODES:
            raise ValueError(
                f"read_handling must be one of {self.READ_HANDLING_MODES}, got {read_handling!r}"
            )
        self.tag_directory_path = Path(tag_directory_path)
        self.read_handling = read_handling
        self.normalization_factor = 1.0

        if normalize_to:
            tag_info_fpath = self.tag_directory_path / "tagInfo.txt"
            with open(tag_info_fpath) as tag_info_file:
                size_fields = tag_info_file.readlines()[1].strip().split("\t")
            num_tags = int(float(size_fields[2]))
            self.normalization_factor = normalize_to / num_tags
            logger.info(
                "Read %d total tags from %s; normalization factor %g",
                num_tags,
                tag_info_fpath,
                self.normalization_factor,
            )

    def _query(self, chrom, start, end):
        coverage = numpy.zeros(end - start)
        tag_fpath = self.tag_directory_path / f"{chrom}.tags.tsv"
        # HOMER tag starts are 1-based: 0-based position `start` is 1-based position start + 1.
        start_offset = utilities.binary_search_tag_file(tag_fpath, start + 1)

        if start_offset is not None:
            with open(tag_fpath, "rb") as tag_file:
                tag_file.seek(start_offset)
                # ponytail: reads starting left of the window are not seen; seeking back by the
                # maximum read length would catch reads that overhang into the window.
                for line in tag_file:
                    fields = line.split()
                    if not fields:
                        continue
                    # HOMER tag starts are 1-based; -1 gives the 0-based start.
                    read_start = int(fields[1]) - 1
                    if read_start >= end:
                        break
                    read_length = 1 if self.read_handling == "starts" else int(fields[4])
                    covered = slice(
                        max(read_start, start) - start, min(read_start + read_length, end) - start
                    )
                    coverage[covered] += float(fields[3])

        return pandas.Series(coverage * self.normalization_factor, index=numpy.arange(start, end))


class BigWigDataSource(_VectorDataSource):
    def __init__(
        self,
        bigwig_fpath: str | Path,
        chromosome_dialect: str = "ucsc",
        fill_value: float = 0.0,
        transform: VectorTransform | None = None,
        smoothing_bandwidth: float = 0,
    ):
        """Read per-base signal from a BigWig file, lazily, one query at a time.

        Parameters
        ----------
        bigwig_fpath : str or Path
            BigWig file.
        chromosome_dialect : str
            Chromosome naming in the file: ``ucsc``, ``ensembl`` or ``yeast``. Queries
            always use UCSC names and are converted.
        fill_value : float
            Value for positions without data and past the chromosome end.
        transform : callable, optional
            Applied to every query result.
        smoothing_bandwidth : float
            Gaussian smoothing standard deviation in bp; 0 disables smoothing.

        Raises
        ------
        FileNotFoundError
            If ``bigwig_fpath`` does not exist.
        ValueError
            If the file is not a BigWig or the dialect is unknown.

        Examples
        --------
        >>> signal = BigWigDataSource("signal.bw").query("chr1", 1000000, 1001000)
        """
        super().__init__(transform=transform, smoothing_bandwidth=smoothing_bandwidth)
        if not os.path.isfile(bigwig_fpath):
            raise FileNotFoundError(f"BigWig file not found: {bigwig_fpath}")
        if chromosome_dialect not in CHROMOSOME_DIALECTS:
            raise ValueError(
                f"chromosome_dialect must be one of {CHROMOSOME_DIALECTS}, "
                f"got {chromosome_dialect!r}"
            )
        self.bigwig_fpath = Path(bigwig_fpath)
        self.chromosome_dialect = chromosome_dialect
        self.fill_value = fill_value
        with self._open() as bigwig_file:
            self._chromosome_sizes = dict(bigwig_file.chroms())
        logger.info(
            "Opened BigWig %s with %d chromosomes", self.bigwig_fpath, len(self._chromosome_sizes)
        )

    @contextlib.contextmanager
    def _open(self) -> Iterator["pyBigWig.pyBigWig"]:
        """Open the BigWig for one operation; a fresh handle per query is cheap and thread-safe."""
        bigwig_file = pyBigWig.open(str(self.bigwig_fpath))
        if bigwig_file is None:
            raise ValueError(f"Could not open as BigWig file: {self.bigwig_fpath}")
        try:
            yield bigwig_file
        finally:
            bigwig_file.close()

    def _file_chrom(self, chrom: str) -> str:
        """Convert a UCSC chromosome name to the file's name for it.

        Raises
        ------
        KeyError
            If neither the converted nor the original name is in the file.
        """
        converted = utilities.convert_chromosome_name(chrom, dialect=self.chromosome_dialect)
        for candidate in (converted, chrom):
            if candidate in self._chromosome_sizes:
                return candidate
        raise KeyError(
            f"Chromosome {chrom!r} (converted to {converted!r}) not found in "
            f"{self.bigwig_fpath}. Available chromosomes: {list(self._chromosome_sizes)[:10]}..."
        )

    def _query(self, chrom, start, end):
        if start >= end:
            raise ValueError(f"start ({start}) must be less than end ({end})")
        if start < 0:
            raise ValueError(f"start must be non-negative, got {start}")

        file_chrom = self._file_chrom(chrom)
        values = numpy.full(end - start, self.fill_value, dtype=numpy.float64)
        data_end = min(end, self._chromosome_sizes[file_chrom])
        if start < data_end:
            with self._open() as bigwig_file:
                file_values = numpy.array(
                    bigwig_file.values(file_chrom, start, data_end), dtype=numpy.float64
                )
            values[: data_end - start] = numpy.where(
                numpy.isnan(file_values), self.fill_value, file_values
            )
        return pandas.Series(values, index=numpy.arange(start, end))

    @property
    def chromosomes(self) -> dict[str, int]:
        """Chromosome names (in the file's dialect) and sizes in bp."""
        return dict(self._chromosome_sizes)

    def __repr__(self) -> str:
        """Return a short description of this data source."""
        return (
            f"BigWigDataSource(path='{self.bigwig_fpath}', dialect='{self.chromosome_dialect}', "
            f"chromosomes={len(self._chromosome_sizes)})"
        )


class IntervalData:
    HOMER_PEAKFILE_NAMES = ("chrom", "chromStart", "chromEnd", "strand", "normed_tag_count")
    HOMER_ANNOTATEDPEAKS_COLUMN_RENAMER = {
        "Chr": "chrom",
        "Start": "chromStart",
        "End": "chromEnd",
        "Strand": "strand",
    }
    BED_COLUMNS = ("chrom", "chromStart", "chromEnd", "name", "score", "strand")
    REQUIRED_COLUMNS = ["chrom", "chromStart", "chromEnd"]
    FILE_FORMATS = ("bed", "homer", "homer_annotated")

    def __init__(self, interval_data: pandas.DataFrame | str | Path, format: str = "bed"):
        """Genomic intervals, standardized to BED column names in ``self.data``.

        Parameters
        ----------
        interval_data : pandas.DataFrame, str or Path
            Either a DataFrame or a file. A DataFrame with ``chrom``, ``chromStart`` and
            ``chromEnd`` columns is used as is; one without them is read as BED columns
            in order (5 columns means ``name`` was moved to the index). A file is read
            according to ``format``.
        format : str
            File format: ``bed``, ``homer`` (HOMER peak file) or ``homer_annotated``
            (HOMER annotatePeaks.pl output). Ignored for DataFrames.

        Raises
        ------
        ValueError
            If ``format`` is unknown or an unnamed DataFrame has fewer than 3 columns.
        """
        if isinstance(interval_data, pandas.DataFrame):
            data = self._from_frame(interval_data)
        else:
            data = self._from_file(interval_data, format)
        self.data = data.sort_values(["chrom", "chromStart"])

    def query(self, chrom: str, start: int, end: int) -> pandas.DataFrame:
        """Return the intervals on ``chrom`` that overlap the 0-based half-open [start, end).

        Parameters
        ----------
        chrom : str
            Chromosome, as named in the interval data.
        start, end : int
            Region, 0-based half-open.

        Returns
        -------
        pandas.DataFrame
            Overlapping rows of ``self.data``.
        """
        data = self.data
        return data.loc[
            (data["chrom"] == chrom) & (data["chromStart"] < end) & (data["chromEnd"] > start)
        ]

    @classmethod
    def _from_frame(cls, frame: pandas.DataFrame) -> pandas.DataFrame:
        """Standardize a DataFrame, naming BED columns by position if they are unnamed."""
        if set(cls.REQUIRED_COLUMNS).issubset(frame.columns):
            return frame
        logger.info("Treating the interval DataFrame as BED-style columns without names")
        if frame.shape[1] < 3:
            raise ValueError(f"Interval data needs at least 3 columns, got {frame.shape[1]}")
        data = frame.copy()
        if frame.shape[1] == 5:  # name has been moved to the index
            data.columns = ["chrom", "chromStart", "chromEnd", "score", "strand"]
        else:
            data.columns = cls._bed_column_names(frame.shape[1])
            if "name" in data.columns:
                data.index = data["name"]
        data.index.name = "IntervalID"
        return data

    @classmethod
    def _bed_column_names(cls, num_columns: int) -> list[str]:
        """Return BED column names for the first columns and ``column_<i>`` for the rest."""
        return list(cls.BED_COLUMNS[:num_columns]) + [
            f"column_{column_idx}" for column_idx in range(len(cls.BED_COLUMNS), num_columns)
        ]

    @classmethod
    def _from_file(cls, interval_fpath: str | Path, format: str) -> pandas.DataFrame:
        """Read a BED, HOMER peak or HOMER annotated peak file."""
        if format not in cls.FILE_FORMATS:
            raise ValueError(
                f"Unknown interval file format {format!r}; expected {cls.FILE_FORMATS}"
            )
        logger.info("Reading %s intervals from %s", format, interval_fpath)

        if format == "bed":
            data = pandas.read_csv(interval_fpath, sep="\t", comment="#", header=None)
            data.columns = cls._bed_column_names(data.shape[1])
            if "name" in data.columns:
                data = data.set_index("name")
        elif format == "homer":
            data = pandas.read_csv(interval_fpath, sep="\t", index_col=0, comment="#", header=None)
            data.columns = (
                list(cls.HOMER_PEAKFILE_NAMES) + list(data.columns)[len(cls.HOMER_PEAKFILE_NAMES) :]
            )
            data.index.name = "peak_id"
        else:
            data = pandas.read_csv(interval_fpath, index_col=0, sep="\t")
            data.index.name = data.index.name.split(" ")[0]
            data = data.rename(columns=cls.HOMER_ANNOTATEDPEAKS_COLUMN_RENAMER)

        logger.info("Read %d intervals from %s", len(data), interval_fpath)
        return data


GeneModelQueryResult = tuple[dict[str, dict], dict[str, dict], dict[str, dict], dict[str, list]]


class _GeneModels(abc.ABC):
    """Base class for sources of gene models (genes, transcripts, exons, UTRs)."""

    @abc.abstractmethod
    def _query(self, chrom: str, start: int, end: int) -> GeneModelQueryResult:
        """Return gene models overlapping a 0-based half-open region."""

    def query(self, chrom: str, start: int, end: int) -> GeneModelQueryResult:
        """Return gene models overlapping a region.

        Parameters
        ----------
        chrom : str
            Chromosome, UCSC style.
        start, end : int
            Region, 0-based half-open.

        Returns
        -------
        tuple of (dict, dict, dict, dict)
            ``genes``, ``transcripts`` and ``components`` keyed by ID, each record holding
            ``contig``, 0-based half-open ``start``/``end``, ``strand`` and the GFF3
            attributes; genes list their ``transcripts`` and transcripts their
            ``components``. The fourth dict maps gene ``Name`` to gene IDs.
        """
        return self._query(chrom, start, end)


def _parse_gff3_attributes(attribute_column: str) -> dict[str, str]:
    """Parse a GFF3 column-9 string into a dict of %-decoded values.

    Parameters
    ----------
    attribute_column : str
        ``key=value`` pairs separated by ``;``.

    Returns
    -------
    dict of str to str
        Attributes. Empty fields (e.g. from a trailing ``;``) are skipped and only the
        first ``=`` separates a key from its value.
    """
    pairs = (field.split("=", 1) for field in attribute_column.strip().split(";") if "=" in field)
    return {key.strip(): unquote(value) for key, value in pairs}


class Gff3Annotations(_GeneModels):
    def __init__(
        self,
        gff3_fpath: str | Path,
        incoming_chromosome_name_converter: Callable[[str], str] = lambda chrom: (
            utilities.convert_chromosome_name(chrom, dialect="ensembl")
        ),
        feature_sources: tuple[str, ...] = DEFAULT_FEATURE_SOURCES,
        gene_types: tuple[str, ...] = DEFAULT_GENE_TYPES,
        transcript_types: tuple[str, ...] = DEFAULT_TRANSCRIPT_TYPES,
        component_types: tuple[str, ...] = DEFAULT_COMPONENT_TYPES,
    ):
        """Gene models from a bgzipped, tabix-indexed GFF3 file (Ensembl layout).

        Prepare a plain GFF3 with ``utilities.bgzip_gff`` and
        ``utilities.generate_tabix_index``.

        Parameters
        ----------
        gff3_fpath : str or Path
            Bgzipped GFF3 with a ``.tbi`` index next to it.
        incoming_chromosome_name_converter : callable
            Converts a UCSC query name to the file's contig name (default: Ensembl).
        feature_sources : tuple of str
            Column-2 sources to keep.
        gene_types, transcript_types, component_types : tuple of str
            Column-3 feature types treated as genes, transcripts and drawn components.
            A component is linked to a gene only through a parent transcript whose type
            is in ``transcript_types``; others are logged as orphans and not drawn. Add
            the file's transcript types here if exons are missing.
        """
        self.gff3_fpath = Path(gff3_fpath)
        self.tabix_file = pysam.TabixFile(str(gff3_fpath))
        self.incoming_chromosome_name_converter = incoming_chromosome_name_converter
        self.feature_sources = feature_sources
        self.gene_types = gene_types
        self.transcript_types = transcript_types
        self.component_types = component_types
        logger.info("Opened tabix-indexed GFF3 %s", self.gff3_fpath)

    @staticmethod
    def _feature_record(split_line: list[str]) -> dict:
        """Build a feature record with 0-based half-open coordinates and its attributes."""
        record = {
            "contig": split_line[0],
            # GFF3 is 1-based closed; -1 on start gives 0-based half-open (end unchanged).
            "start": int(split_line[3]) - 1,
            "end": int(split_line[4]),
            "strand": split_line[6],
        }
        record.update(_parse_gff3_attributes(split_line[8]))
        return record

    def _query(self, chrom, start, end):
        genes, transcripts, components = {}, {}, {}
        gene_names_to_ids = {}

        # pysam fetch takes 0-based half-open coordinates, the same as the query.
        query_rows = self.tabix_file.fetch(
            self.incoming_chromosome_name_converter(chrom), start, end
        )
        for line in query_rows:
            split_line = line.strip("\n").split("\t")
            source, feature_type = split_line[1], split_line[2]
            if source not in self.feature_sources:
                continue
            record = self._feature_record(split_line)

            if feature_type in self.gene_types:
                record["Name"] = record.get("Name", record["ID"])
                record["transcripts"] = []
                genes[record["ID"]] = record
                gene_names_to_ids.setdefault(record["Name"], []).append(record["ID"])
            elif feature_type in self.transcript_types:
                record["components"] = []
                transcripts[record["ID"]] = record
            elif feature_type in self.component_types:
                record["type"] = feature_type
                components[record.get("exon_id", str(len(components)))] = record

        self._link_children(genes, transcripts, "transcripts", "transcript")
        self._link_children(transcripts, components, "components", "component")
        logger.info(
            "Found %d genes, %d transcripts, %d components in %s:%d-%d",
            len(genes),
            len(transcripts),
            len(components),
            chrom,
            start,
            end,
        )
        return genes, transcripts, components, gene_names_to_ids

    @staticmethod
    def _link_children(
        parents: dict[str, dict], children: dict[str, dict], child_list_key: str, child_kind: str
    ) -> None:
        """Append each child's ID to the ``child_list_key`` list of every parent it names."""
        # ponytail: Parent is split on ',' after %-decoding, so an ID with an encoded %2C
        # mis-splits; split before decoding if such IDs turn up.
        for child_id, child_data in children.items():
            for parent_id in child_data.get("Parent", "").split(","):
                if parent_id in parents:
                    parents[parent_id][child_list_key].append(child_id)
                else:
                    logger.warning(
                        "Orphan %s %s with missing parent %s", child_kind, child_id, parent_id
                    )


class _MatrixData(abc.ABC):
    """Base class for binned contact matrices. Subclasses set ``bin_size``."""

    bin_size: int

    @abc.abstractmethod
    def _query(self, chrom: str, start: int, end: int) -> pandas.DataFrame:
        """Return the matrix block for a region."""

    def query(self, chrom: str, start: int, end: int) -> pandas.DataFrame:
        """Return the square block of bins from ``start`` to ``end``, rounded to bins.

        Parameters
        ----------
        chrom : str
            Chromosome, as named in the matrix source.
        start, end : int
            Region in bp; both ends are rounded to the nearest bin and included.

        Returns
        -------
        pandas.DataFrame
            Matrix block indexed and labeled by bin start.
        """
        return self._query(chrom, start, end)


class HicDataDir(_MatrixData):
    def __init__(self, fpath_template: str, bin_size: int = 10000):
        """Hi-C matrices stored as one TSV file per chromosome.

        Parameters
        ----------
        fpath_template : str
            Path with one ``{}`` placeholder for the chromosome name. Each file is a
            square matrix with bin starts as row and column labels.
        bin_size : int
            Bin size in bp.
        """
        self.fpath_template = fpath_template
        self.bin_size = bin_size

    def _query(self, chrom, start, end):
        chrom_fpath = self.fpath_template.format(chrom)
        chrom_data = pandas.read_csv(chrom_fpath, sep="\t", index_col=0)
        chrom_data.columns = chrom_data.columns.astype(int)
        logger.info("Read a %d x %d matrix from %s", *chrom_data.shape, chrom_fpath)

        rounded_start = utilities.roundto(start, self.bin_size)
        rounded_end = utilities.roundto(end, self.bin_size)
        return chrom_data.loc[rounded_start:rounded_end, rounded_start:rounded_end]


class HicDataDict(_MatrixData):
    def __init__(self, data_dict: dict[str, pandas.DataFrame], bin_size: int):
        """Hi-C matrices held in memory, one DataFrame per chromosome.

        Parameters
        ----------
        data_dict : dict of str to pandas.DataFrame
            Square matrices with bin starts as row and column labels.
        bin_size : int
            Bin size in bp.
        """
        self.data_dict = data_dict
        self.bin_size = bin_size

    def _query(self, chrom, start, end):
        rounded_start = utilities.roundto(start, self.bin_size)
        rounded_end = utilities.roundto(end, self.bin_size)
        return self.data_dict[chrom].loc[rounded_start:rounded_end, rounded_start:rounded_end]

    @staticmethod
    def rename_hic_df(hic_df: pandas.DataFrame) -> pandas.DataFrame:
        """Relabel ``<chrom>-<bin start>`` row and column labels as integer bin starts.

        Parameters
        ----------
        hic_df : pandas.DataFrame
            Matrix with string labels like ``chr1-10000``. Modified in place.

        Returns
        -------
        pandas.DataFrame
            The same DataFrame.
        """
        hic_df.index = [int(name.split("-")[1]) for name in hic_df.index]
        hic_df.columns = [int(name.split("-")[1]) for name in hic_df.columns]
        return hic_df
