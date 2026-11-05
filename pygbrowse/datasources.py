import contextlib
import logging
import os
from urllib.parse import unquote

import numpy
import pandas
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
# DEFAULT_MAXIMUM_TRANSCRIPT_SUPPORT = 5

# ToDo: For each class, allow option of loading into memory or leaving on disk (where applicable)
# ToDo: Add indexing of on-disk csv-like files
# ToDo: Refactor vector data sources to transparently interpolate sparse vectors. Probably will have to drop dict-of-series interface.


class _ChromWrapper:
    def __init__(self, chrom, parent_data_source):
        self.chrom = chrom
        self.parent_data_source = parent_data_source

    def __getitem__(self, key):
        # print(key)
        # ToDo: Add support for step argument
        try:  # See if key is a slice
            query_start = key.start
            query_end = key.stop
        except TypeError:  # if not, treat as a scalar index
            query_start = key
            query_end = key + 1

        return self.parent_data_source.query(
            query_chrom=self.chrom, query_start=query_start, query_end=query_end
        )


class _DataVector:
    def __init__(self, chrom, parent_data_source):
        self.loc = _ChromWrapper(chrom=chrom, parent_data_source=parent_data_source)


class _VectorDataSource:
    # ToDo: Add methods for arithmetic and such, as done for old Pileups class
    def __init__(self, transform=None, smoothing_bandwidth=0):
        self.transform = transform
        if smoothing_bandwidth:
            self.convolution_kernel = utilities.gaussian_kernel(smoothing_bandwidth)
        else:
            self.convolution_kernel = None

    def _query(self, query_chrom, query_start, query_end):
        print("Stub method -- must be overridden by inheritors")

    def query(self, query_chrom, query_start, query_end):
        query_result = self._query(
            query_chrom=query_chrom, query_start=query_start, query_end=query_end
        )

        if self.convolution_kernel is not None:
            query_result = pandas.Series(
                convolve(query_result, self.convolution_kernel, mode="same"),
                index=query_result.index,
            )
        if self.transform:
            query_result = self.transform(query_result)

        return query_result

    def __getitem__(self, key):
        return _DataVector(chrom=key, parent_data_source=self)


class SparseVectors(_VectorDataSource):
    def __init__(self, series_dict, transform=None, convolution_kernel=None):
        self.data = series_dict
        self.transform = transform
        self.convolution_kernel = convolution_kernel

    def _query(self, query_chrom, query_start, query_end):
        this_chrom_vector = self.data[query_chrom]
        # One point either side of the window is kept so that line plots reach the edges.
        start_ipos = max(numpy.searchsorted(this_chrom_vector.index, query_start) - 1, 0)
        end_ipos = numpy.searchsorted(this_chrom_vector.index, query_end) + 1

        return this_chrom_vector.iloc[start_ipos:end_ipos]


class TagDirectory(_VectorDataSource):
    tag_strand_translator = {0: "+", 1: "-"}
    READ_HANDLING_MODES = ("reads", "starts")

    def __init__(
        self,
        tag_directory_path,
        normalize_to=DEFAULT_TAG_COUNT_NORMALIZATION_TARGET,
        transform=None,
        smoothing_bandwidth=0,
        read_handling="reads",
    ):
        super().__init__(transform=transform, smoothing_bandwidth=smoothing_bandwidth)

        if read_handling not in self.READ_HANDLING_MODES:
            raise ValueError(
                f"read_handling must be one of {self.READ_HANDLING_MODES}, got {read_handling!r}"
            )
        self.tag_directory_path = tag_directory_path
        self.read_handling = read_handling
        self.normalization_factor = 1.0

        if normalize_to:
            tag_info_fpath = os.path.join(tag_directory_path, "tagInfo.txt")
            with open(tag_info_fpath) as tag_info_file:
                sizeline = tag_info_file.readlines()[1].strip().split("\t")
            num_tags = int(float(sizeline[2]))
            self.normalization_factor = normalize_to / num_tags

    def _query(self, query_chrom, query_start, query_end):
        coverage = numpy.zeros(query_end - query_start)
        tag_fpath = os.path.join(self.tag_directory_path, f"{query_chrom}.tags.tsv")
        # HOMER tag starts are 1-based: the 0-based query_start is 1-based query_start + 1.
        start_offset = utilities.binary_search_tag_file(tag_fpath, query_start + 1)

        if start_offset is not None:
            with open(tag_fpath, "rb") as tag_file:
                tag_file.seek(start_offset)
                # ponytail: reads starting left of the window are not seen; seeking back by the
                # maximum read length would catch reads that overhang into the window.
                for line in tag_file:
                    fields = line.split(b"\t")
                    # HOMER tag starts are 1-based; -1 gives the 0-based start.
                    read_start = int(fields[1]) - 1
                    if read_start >= query_end:
                        break
                    if self.read_handling == "starts":
                        read_end = read_start + 1
                    else:
                        read_end = read_start + int(fields[4])
                    coverage[
                        max(read_start, query_start) - query_start : min(read_end, query_end)
                        - query_start
                    ] += float(fields[3])

        return pandas.Series(
            coverage * self.normalization_factor, index=numpy.arange(query_start, query_end)
        )


class BigWigDataSource(_VectorDataSource):
    """
    Data source for reading signal data from BigWig files.

    BigWig files are indexed binary files containing per-base signal values
    commonly used for storing coverage, ChIP-seq signal, and other genomic
    continuous data.

    This class provides lazy loading of data - values are read from the BigWig
    file only when queried, making it memory-efficient for large files.

    Parameters
    ----------
    bigwig_path : str
        Path to the BigWig file to read from.
    chromosome_dialect : str, optional
        The chromosome naming convention used in the BigWig file.
        One of 'ucsc' (e.g., 'chr1'), 'ensembl' (e.g., '1'), or 'yeast' (Roman numerals).
        Default is 'ucsc'.
    fill_value : float, optional
        Value to use for positions with no data in the BigWig file.
        Default is 0.0 (not NaN, to ensure smooth plotting).
    transform : callable, optional
        A function to apply to query results. Default is None.
    smoothing_bandwidth : int, optional
        Standard deviation for Gaussian smoothing kernel. 0 means no smoothing.
        Default is 0.

    Examples
    --------
    >>> from pygbrowse import datasources
    >>> bw = datasources.BigWigDataSource('/path/to/signal.bw')
    >>> signal = bw['chr1'].loc[1000000:1001000]
    >>>
    >>> # With chromosome conversion for Ensembl-named BigWig
    >>> bw_ensembl = datasources.BigWigDataSource(
    ...     '/path/to/ensembl_signal.bw',
    ...     chromosome_dialect='ensembl'
    ... )
    """

    def __init__(
        self,
        bigwig_path: str,
        chromosome_dialect: str = "ucsc",
        fill_value: float = 0.0,
        transform=None,
        smoothing_bandwidth: int = 0,
    ):
        """
        Initialize a BigWigDataSource.

        Parameters
        ----------
        bigwig_path : str
            Path to the BigWig file. Must exist and be readable.
        chromosome_dialect : str
            Chromosome naming convention in the BigWig file.
            Valid values: 'ucsc', 'ensembl', 'yeast'
        fill_value : float
            Value to substitute for NaN positions.
        transform : callable or None
            Optional transformation function for query results.
        smoothing_bandwidth : int
            Gaussian smoothing standard deviation. 0 = no smoothing.

        Raises
        ------
        FileNotFoundError
            If bigwig_path does not exist.
        ValueError
            If bigwig_path exists but is not a valid BigWig file.
            If chromosome_dialect is not one of the valid options.
        """
        super(BigWigDataSource, self).__init__(
            transform=transform, smoothing_bandwidth=smoothing_bandwidth
        )

        # Validate file exists before attempting to open
        if not os.path.isfile(bigwig_path):
            raise FileNotFoundError(f"BigWig file not found: {bigwig_path}")

        # Validate chromosome dialect
        valid_dialects = ("ucsc", "ensembl", "yeast")
        if chromosome_dialect not in valid_dialects:
            raise ValueError(
                f"chromosome_dialect must be one of {valid_dialects}, got '{chromosome_dialect}'"
            )

        self.bigwig_path = bigwig_path
        self.chromosome_dialect = chromosome_dialect
        self.fill_value = fill_value

        # Validate file is actually a BigWig by attempting to open it
        self._validate_and_cache_chroms()

    def _validate_and_cache_chroms(self) -> None:
        """
        Validate the BigWig file and cache chromosome information.

        Called during initialization to fail fast if the file is invalid.

        Raises
        ------
        ValueError
            If the file cannot be opened as a BigWig file.
        """
        import pyBigWig

        handle = pyBigWig.open(str(self.bigwig_path))
        if handle is None:
            raise ValueError(f"Could not open as BigWig file: {self.bigwig_path}")

        try:
            self._chromosome_sizes = dict(handle.chroms())
        finally:
            handle.close()

    @contextlib.contextmanager
    def _get_bigwig_handle(self):
        """
        Context manager for safely accessing the BigWig file.

        Opens a new handle for each query to ensure thread safety.
        pyBigWig handles are lightweight so this is acceptable.

        Yields
        ------
        pyBigWig.pyBigWig
            Open BigWig file handle
        """
        import pyBigWig

        handle = pyBigWig.open(str(self.bigwig_path))
        if handle is None:
            raise ValueError(f"Failed to open BigWig file: {self.bigwig_path}")
        try:
            yield handle
        finally:
            handle.close()

    def _convert_query_chrom(self, query_chrom: str) -> str:
        """
        Convert a query chromosome name to the dialect used by the BigWig file.

        The query chromosome name is assumed to be in UCSC format (consistent
        with the rest of pygbrowse).

        Parameters
        ----------
        query_chrom : str
            Chromosome name from query (UCSC format, e.g., 'chr1')

        Returns
        -------
        str
            Chromosome name in the BigWig file's dialect

        Raises
        ------
        KeyError
            If the converted chromosome name is not present in the BigWig file.
        """
        # Convert from UCSC (query format) to the BigWig file's dialect
        converted_chrom = utilities.convert_chromosome_name(
            chrom_string=query_chrom, dialect=self.chromosome_dialect
        )

        # Validate chromosome exists in the file
        if converted_chrom not in self._chromosome_sizes:
            # Try the original name as fallback
            if query_chrom in self._chromosome_sizes:
                return query_chrom
            raise KeyError(
                f"Chromosome '{query_chrom}' (converted to '{converted_chrom}') "
                f"not found in BigWig file. Available chromosomes: "
                f"{list(self._chromosome_sizes.keys())[:10]}..."
            )

        return converted_chrom

    def _query(self, query_chrom: str, query_start: int, query_end: int) -> pandas.Series:
        """
        Query the BigWig file for signal values in a genomic region.

        Parameters
        ----------
        query_chrom : str
            Chromosome name (UCSC format expected, e.g., 'chr1')
        query_start : int
            Start position (0-based, inclusive)
        query_end : int
            End position (0-based, exclusive)

        Returns
        -------
        pandas.Series
            Series with genomic positions as index and signal values as data.
            Positions without data in the BigWig are filled with self.fill_value.

        Raises
        ------
        KeyError
            If chromosome is not found in the BigWig file.
        ValueError
            If query coordinates are invalid (e.g., start >= end).
        """
        # Validate query coordinates
        if query_start >= query_end:
            raise ValueError(
                f"query_start ({query_start}) must be less than query_end ({query_end})"
            )

        if query_start < 0:
            raise ValueError(f"query_start must be non-negative, got {query_start}")

        # Convert chromosome name to BigWig file's dialect
        bigwig_chrom = self._convert_query_chrom(query_chrom)

        # Clamp query_end to chromosome size to avoid pyBigWig errors
        chrom_size = self._chromosome_sizes[bigwig_chrom]
        effective_end = min(query_end, chrom_size)
        effective_start = min(query_start, chrom_size)

        # Handle case where query is entirely beyond chromosome end
        if effective_start >= chrom_size:
            return pandas.Series(
                numpy.full(query_end - query_start, self.fill_value),
                index=numpy.arange(query_start, query_end),
            )

        # Get values from BigWig file
        with self._get_bigwig_handle() as bw:
            values = bw.values(bigwig_chrom, effective_start, effective_end)

        # Convert to numpy array and handle NaN values
        values = numpy.array(values, dtype=numpy.float64)
        nan_mask = numpy.isnan(values)
        if nan_mask.any():
            values[nan_mask] = self.fill_value

        # Build the result series
        result_index = numpy.arange(effective_start, effective_end)
        result = pandas.Series(values, index=result_index)

        # If query extended beyond chromosome, pad with fill_value
        if effective_end < query_end:
            padding_index = numpy.arange(effective_end, query_end)
            padding = pandas.Series(
                numpy.full(len(padding_index), self.fill_value), index=padding_index
            )
            result = pandas.concat([result, padding])

        return result

    @property
    def chromosomes(self) -> dict:
        """
        Return dictionary of chromosome names and sizes from the BigWig file.

        Returns
        -------
        dict
            Keys are chromosome names (in the BigWig file's dialect),
            values are chromosome sizes in base pairs.
        """
        return dict(self._chromosome_sizes)

    def __repr__(self) -> str:
        """Return a string representation of this data source."""
        return (
            f"BigWigDataSource("
            f"path='{self.bigwig_path}', "
            f"dialect='{self.chromosome_dialect}', "
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

    def __init__(self, interval_data, format="bed"):
        """
        Loads genomic interval information in various formats and stores them in a standardized form as a
        pandas.DataFrame in self.data.

        :param:`interval_data` should be a pandas.DataFrame representing BED-formatted genomic data, or,
        alternatively, a filename pointing to one of the following file formats:

        * A BED file
        * A HOMER peak file
        * A HOMER annotated peak file.

        If a filename is passed instead of a DataFrame, :param:`format` should be specified. Allowed values are:
            'bed', 'homer', 'homer_annotated'

        :param interval_data:
        :param format:
        """
        if isinstance(interval_data, pandas.DataFrame):
            data = self._from_frame(interval_data)
        else:
            data = self._from_file(interval_data, format)
        self.data = data.sort_values(["chrom", "chromStart"])

    def query(self, chrom, start, end):
        """Return the intervals on ``chrom`` that overlap the 0-based half-open [start, end)."""
        data = self.data
        return data.loc[
            (data["chrom"] == chrom) & (data["chromStart"] < end) & (data["chromEnd"] > start)
        ]

    @classmethod
    def _from_frame(cls, frame):
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
    def _bed_column_names(cls, num_columns):
        return list(cls.BED_COLUMNS[:num_columns]) + [
            f"column_{column_idx}" for column_idx in range(len(cls.BED_COLUMNS), num_columns)
        ]

    @classmethod
    def _from_file(cls, interval_fpath, format):
        if format not in cls.FILE_FORMATS:
            raise ValueError(f"Unknown interval file format {format!r}; expected {cls.FILE_FORMATS}")
        logger.info("Reading %s intervals from %s", format, interval_fpath)

        if format == "bed":
            data = pandas.read_csv(interval_fpath, sep="\t", comment="#", header=None)
            data.columns = cls._bed_column_names(data.shape[1])
            if "name" in data.columns:
                data = data.set_index("name")
        elif format == "homer":
            data = pandas.read_csv(interval_fpath, sep="\t", index_col=0, comment="#", header=None)
            data.columns = (
                list(cls.HOMER_PEAKFILE_NAMES)
                + list(data.columns)[len(cls.HOMER_PEAKFILE_NAMES) :]
            )
            data.index.name = "peak_id"
        else:
            data = pandas.read_csv(interval_fpath, index_col=0, sep="\t")
            data.index.name = data.index.name.split(" ")[0]
            data = data.rename(columns=cls.HOMER_ANNOTATEDPEAKS_COLUMN_RENAMER)

        logger.info("Read %d intervals from %s", len(data), interval_fpath)
        return data


class _GeneModels:
    def __init__(self):
        """
        Superclass for data sources that describe gene models (gene boundaries, transcript
        boundaries, exons, introns, UTRs, etc.).
        """
        pass

    def _query(self, query_chromosome, query_start, query_end):
        print("Must be overridden by inheritors!")

    def query(self, chromosome, start, end):
        return self._query(query_chromosome=chromosome, query_start=start, query_end=end)


def _parse_gff3_attributes(attribute_column):
    """Parse a GFF3 column-9 string into a dict of %-decoded values.

    Empty fields (e.g. from a trailing ';') are skipped and only the first '=' splits a
    key from its value.
    """
    pairs = (field.split("=", 1) for field in attribute_column.strip().split(";") if "=" in field)
    return {key.strip(): unquote(value) for key, value in pairs}


class Gff3Annotations(_GeneModels):
    def __init__(
        self,
        gff3_filename,
        incoming_chromosome_name_converter=lambda x: utilities.convert_chromosome_name(
            x, dialect="ensembl"
        ),
        feature_sources=DEFAULT_FEATURE_SOURCES,
        gene_types=DEFAULT_GENE_TYPES,
        transcript_types=DEFAULT_TRANSCRIPT_TYPES,
        component_types=DEFAULT_COMPONENT_TYPES,
    ):
        super().__init__()

        self.tabix_file = pysam.TabixFile(str(gff3_filename))
        self.incoming_chromosome_name_converter = incoming_chromosome_name_converter
        self.feature_sources = feature_sources
        self.gene_types = gene_types
        self.transcript_types = transcript_types
        self.component_types = component_types

    def _feature_record(self, split_line):
        fields = _parse_gff3_attributes(split_line[8])
        record = {
            "contig": split_line[0],
            # GFF3 is 1-based closed; -1 on start gives 0-based half-open (end unchanged).
            "start": int(split_line[3]) - 1,
            "end": int(split_line[4]),
            "strand": split_line[6],
        }
        record.update(fields)
        return record

    def _query(self, query_chromosome, query_start, query_end):
        genes, transcripts, components = {}, {}, {}
        gene_names_to_ensembl_ids = {}

        query_rows = self.tabix_file.fetch(
            self.incoming_chromosome_name_converter(query_chromosome), query_start, query_end
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
                gene_names_to_ensembl_ids.setdefault(record["Name"], []).append(record["ID"])
            elif feature_type in self.transcript_types:
                record["components"] = []
                transcripts[record["ID"]] = record
            elif feature_type in self.component_types:
                record["type"] = feature_type
                components[record.get("exon_id", str(len(components)))] = record

        self._link_children(genes, transcripts, "transcripts", "transcript")
        self._link_children(transcripts, components, "components", "component")
        return genes, transcripts, components, gene_names_to_ensembl_ids

    @staticmethod
    def _link_children(parents, children, child_list_key, child_kind):
        # ponytail: Parent is split on ',' after %-decoding, so an ID with an encoded %2C
        # mis-splits; split before decoding if such IDs turn up.
        for child_id, child_data in children.items():
            for parent_id in child_data.get("Parent", "").split(","):
                if parent_id in parents:
                    parents[parent_id][child_list_key].append(child_id)
                else:
                    logger.warning("Orphan %s %s with missing parent %s", child_kind, child_id, parent_id)


class _MatrixData:
    def __init__(self):
        pass

    def _query(self):
        print("Must be overridden by inheritors")

    def query(self, chrom, start, end):
        return self._query(query_chrom=chrom, query_start=start, query_end=end)


class HicDataDir(_MatrixData):
    def __init__(self, fname_template, bin_size=10000):
        self.fname_template = fname_template
        self.bin_size = bin_size

    def _query(self, query_chrom, query_start, query_end):
        this_chrom_fname = self.fname_template.format(query_chrom)
        this_chrom_data = pandas.read_csv(this_chrom_fname, sep="\t", index_col=0)
        this_chrom_data.columns = this_chrom_data.columns.astype(int)

        rounded_start = utilities.roundto(query_start, self.bin_size)
        rounded_end = utilities.roundto(query_end, self.bin_size)

        return this_chrom_data.loc[rounded_start:rounded_end, rounded_start:rounded_end]


class HicDataDict(_MatrixData):
    def __init__(self, data_dict, bin_size):
        # self.data_dict = {chrom:self.rename_hic_df(data_dict[chrom]) for chrom in data_dict}
        self.data_dict = data_dict
        self.bin_size = bin_size

    def _query(self, query_chrom, query_start, query_end):
        rounded_start = utilities.roundto(query_start, self.bin_size)
        rounded_end = utilities.roundto(query_end, self.bin_size)

        return self.data_dict[query_chrom].loc[rounded_start:rounded_end, rounded_start:rounded_end]

    @staticmethod
    def rename_hic_df(hic_df):
        hic_df.index = [int(name.split("-")[1]) for name in hic_df.index]
        hic_df.columns = [int(name.split("-")[1]) for name in hic_df.columns]
        return hic_df
