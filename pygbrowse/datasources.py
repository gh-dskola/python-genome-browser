import contextlib
import os

import numpy
import pandas
import pysam
from scipy.signal import convolve

from . import utilities
from .utilities import log_print

DEFAULT_TAG_COUNT_NORMALIZATION_TARGET = 10000000
DEFAULT_FEATURE_SOURCES = ('ensembl', 'havana', 'ensembl_havana')
DEFAULT_GENE_TYPES = (
'gene', 'RNA', 'mt_gene', 'lincRNA_gene', 'miRNA_gene', 'ncRNA_gene', 'rRNA_gene', 'snRNA_gene', 'snoRNA_gene',
'processed_transcript')
DEFAULT_TRANSCRIPT_TYPES = ('mRNA', 'transcript', 'lincRNA', 'lnc_RNA', 'miRNA', 'ncRNA', 'snRNA', 'snoRNA')
DEFAULT_COMPONENT_TYPES = ('CDS', 'three_prime_UTR', 'five_prime_UTR')
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

        return self.parent_data_source.query(query_chrom=self.chrom, query_start=query_start, query_end=query_end)


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
        print('Stub method -- must be overridden by inheritors')

    def query(self, query_chrom, query_start, query_end):
        query_result = self._query(query_chrom=query_chrom, query_start=query_start, query_end=query_end)
        
        if self.convolution_kernel is not None:
            query_result = pandas.Series(convolve(query_result, self.convolution_kernel, mode='same'), index=query_result.index)
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
        start_ipos = numpy.searchsorted(this_chrom_vector.keys(), query_start) - 1
        end_ipos = numpy.searchsorted(this_chrom_vector.keys(), query_end) + 1
        
        return this_chrom_vector.iloc[start_ipos:end_ipos]
        

class TagDirectory(_VectorDataSource):
    tag_strand_translator = {0: '+', 1: '-'}

    def __init__(self, tag_directory_path, normalize_to=DEFAULT_TAG_COUNT_NORMALIZATION_TARGET, transform=None,
                 smoothing_bandwidth=0):
        super(TagDirectory, self).__init__(transform=transform, smoothing_bandwidth=smoothing_bandwidth)

        self.tag_directory_path = tag_directory_path

        if normalize_to:
            # extract total tag count from tagInfo.txt
            tag_info_fname = os.path.join(tag_directory_path, 'tagInfo.txt')
            with open(tag_info_fname, 'rt') as tag_info_file:
                sizeline = tag_info_file.readlines()[1].strip().split('\t')
            num_tags = int(float(sizeline[2]))

            self.normalization_factor = normalize_to / num_tags

    def _query(self, query_chrom, query_start, query_end, read_handling='reads'):
        # ToDo: Add argument validation to all functions and methods with string parameters
        # ToDo: Add verbosity-based logging output
        # ToDo; Compare performance with memory-mapped pandas DataFrames
        query_result = pandas.Series(numpy.zeros(query_end - query_start), index=numpy.arange(query_start, query_end))

        tag_filename = os.path.join(self.tag_directory_path, '{}.tags.tsv'.format(query_chrom))
        start_offset = utilities.binary_search_tag_file(tag_filename=tag_filename, search_target=query_start + 1)

        done = False
        with open(tag_filename, 'rt') as tag_file:
            tag_file.seek(start_offset)
            # print(start_offset)
            while not done:
                line_fields = tag_file.readline().strip().split('\t')
                # print(line_fields)
                if len(line_fields) > 1:
                    # chrom = line_fields[0]
                    read_start = int(line_fields[1]) - 1
                    # strand = self.tag_strand_translator[int(line_fields[2])]
                    depth = float(line_fields[3])

                    if read_handling == 'starts':
                        assert read_start > query_start
                        if read_start < query_end:
                            query_result.loc[read_start] += depth
                        else:
                            done = True

                    elif read_handling == 'reads':
                        # ToDo: Hard to do this in a streaming fashion because we don't know how far upstream to seek to capture left-overhanging reads.
                        read_len = int(line_fields[4])
                        if query_start < read_start <= query_end or query_start < read_start + read_len <= query_end:
                            # print(max(read_start, query_start), min(read_start + read_len,
                            #                                         query_end))
                            query_result.loc[max(read_start, query_start):min(read_start + read_len,
                                                                              query_end)] += depth  # trim to visible vector
                        else:
                            done = True

        query_result *= self.normalization_factor

        return query_result


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
        chromosome_dialect: str = 'ucsc',
        fill_value: float = 0.0,
        transform=None,
        smoothing_bandwidth: int = 0
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
            transform=transform,
            smoothing_bandwidth=smoothing_bandwidth
        )

        # Validate file exists before attempting to open
        if not os.path.isfile(bigwig_path):
            raise FileNotFoundError(f"BigWig file not found: {bigwig_path}")

        # Validate chromosome dialect
        valid_dialects = ('ucsc', 'ensembl', 'yeast')
        if chromosome_dialect not in valid_dialects:
            raise ValueError(
                f"chromosome_dialect must be one of {valid_dialects}, "
                f"got '{chromosome_dialect}'"
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

        handle = pyBigWig.open(self.bigwig_path)
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

        handle = pyBigWig.open(self.bigwig_path)
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
            chrom_string=query_chrom,
            dialect=self.chromosome_dialect
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

    def _query(
        self,
        query_chrom: str,
        query_start: int,
        query_end: int
    ) -> pandas.Series:
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
                f"query_start ({query_start}) must be less than "
                f"query_end ({query_end})"
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
                index=numpy.arange(query_start, query_end)
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
                numpy.full(len(padding_index), self.fill_value),
                index=padding_index
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
    # HOMER_PEAKFILE_HEADER_ROW = 39
    # HOMER_PEAKFILE_COLUMN_RENAMER = {'chr': 'chrom', 'start': 'chromStart', 'end': 'chromEnd'}
    HOMER_PEAKFILE_NAMES = ('chrom', 'chromStart', 'chromEnd', 'strand', 'normed_tag_count')
    HOMER_ANNOTATEDPEAKS_COLUMN_RENAMER = {'Chr': 'chrom', 'Start': 'chromStart', 'End': 'chromEnd', 'Strand': 'strand'}

    def __init__(self, interval_data, format='bed'):
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
        try:
            _ = interval_data.loc[:, ['chrom', 'chromStart', 'chromEnd', 'strand']]

        except KeyError:  # maybe it's a BED DataFrame without column names?
            log_print('Guessing this is a BED-style DataFrame without column names')

            assert interval_data.shape[1] >= 3, 'Not enough columns (got {})!'.format(interval_data.shape[1])

            if interval_data.shape[1] >= 6:  # assume name is still separate column
                self.data = interval_data.copy()
                self.data.columns = ['chrom', 'chromStart', 'chromEnd', 'name',
                                     'score', 'strand'] + list(self.data.columns)[6:]
                self.data.index = self.data['name']

            elif interval_data.shape[1] == 5:  # assume name has been made the index and deleted from the columns
                self.data = interval_data.copy()
                self.data.columns = ['chrom', 'chromStart', 'chromEnd', 'score',
                                     'strand']
            else:
                self.data = interval_data.copy()
                self.data.columns = ['chrom', 'chromStart', 'chromEnd', 'score',
                                     'strand'][:interval_data.shape[1] - 5]

            self.data.index.name = 'IntervalID'

        except (AttributeError,):  # guessing it's a filename string
            log_print('Guessing {} is a filename'.format(interval_data))
            # if format == 'auto':
            #     extension = filename.split('.')[-1]
            #     if extension.lower() == 'bed':
            #         format = 'bed'
            #     elif extension.lower() == 'homer':
            #         # ToDo: Add more sophisticated methods of detecting formats since, e.g. .txt can refer to many.
            #         format = 'homer'

            if format == 'bed':
                self.data = pandas.read_csv(interval_data, sep='\t', index_col=3, comment='#', header=None,
                                            names=['chrom', 'chromStart', 'chromEnd', 'score', 'strand'])
            elif format == 'homer':
                self.data = pandas.read_csv(interval_data, sep='\t', index_col=0, comment='#', header=None)
                self.data.columns = list(self.HOMER_PEAKFILE_NAMES) + list(self.data.columns)[len(self.HOMER_PEAKFILE_NAMES):]
                self.data.index.name = 'peak_id'
                # self.data = self.data.rename(columns=self.HOMER_PEAKFILE_COLUMN_RENAMER)

            elif format == 'homer_annotated':
                self.data = pandas.read_csv(interval_data, index_col=0, sep='\t')
                self.data.index.name = self.data.index.name.split(' ')[0]
                self.data = self.data.rename(columns=self.HOMER_ANNOTATEDPEAKS_COLUMN_RENAMER)

        else:  # seems to be a properly-formatted DataFrame so just  store it
            self.data = interval_data

        self.data = self.data.sort_values(['chrom', 'chromStart'])


class _GeneModels():
    def __init__(self):
        """
        Superclass for data sources that describe gene models (gene boundaries, transcript
        boundaries, exons, introns, UTRs, etc.).        
        """
        pass

    def _query(self, query_chromosome, query_start, query_end):
        print('Must be overridden by inheritors!')

    def query(self, chromosome, start, end):
        return self._query(query_chromosome=chromosome, query_start=start, query_end=end)


from pygbrowse.datasources import _GeneModels

class Gff3Annotations(_GeneModels):
    def __init__(self,
                 gff3_filename,
                 incoming_chromosome_name_converter=lambda x: utilities.convert_chromosome_name(x, dialect='ensembl'),
                 outgoing_chromosome_name_converter=lambda x: utilities.convert_chromosome_name(x, dialect='ucsc'),
                 feature_sources=DEFAULT_FEATURE_SOURCES,
                 gene_types=DEFAULT_GENE_TYPES,
                 transcript_types=DEFAULT_TRANSCRIPT_TYPES,
                 component_types=DEFAULT_COMPONENT_TYPES,
                 # maximum_transcript_support=DEFAULT_MAXIMUM_TRANSCRIPT_SUPPORT
                 ):

        super(Gff3Annotations, self).__init__()

        self.tabix_file = pysam.TabixFile(gff3_filename)
        self.incoming_chromosome_name_converter = incoming_chromosome_name_converter
        self.outgoing_chromosome_name_converter = outgoing_chromosome_name_converter
        self.feature_sources = feature_sources
        self.gene_types = gene_types
        self.transcript_types = transcript_types
        self.component_types = component_types
        # self.maximum_transcript_support = maximum_transcript_support

    def _query(self, query_chromosome, query_start, query_end):
        gene_names_to_ensembl_ids = {}
        genes = {}
        transcripts = {}
        components = {}
        component_num = 0  # serial index for components without IDs

        query_rows = self.tabix_file.fetch(self.incoming_chromosome_name_converter(query_chromosome), query_start,
                                           query_end)

        for line in query_rows:
            split_line = line.strip('\n').split('\t')
            source, feature_type = split_line[1], split_line[2]

            if source in self.feature_sources:
                contig = split_line[0]
                start = int(split_line[3])
                end = int(split_line[4])
                strand = split_line[6]

                fields = dict(field_value_pair.split('=') for field_value_pair in split_line[8].split(';'))
#                 print(line)

                if feature_type in self.gene_types:
                    ensembl_id = fields['ID']
                    gene_name = fields['Name']
                    #                     assert ensembl_id not in genes, 'Duplicate entry for gene {} on line {}'.format(ensembl_id,
                    #                                                                                                     line_num)

                    genes[ensembl_id] = {'contig': contig,
                                         'start': start - 1,  # convert 1-based to 0-based
                                         'end': end,
                                         'strand': strand,
                                         'transcripts': []}

                    genes[ensembl_id].update(fields)
                    if gene_name not in gene_names_to_ensembl_ids:
                        gene_names_to_ensembl_ids[gene_name] = []
                    gene_names_to_ensembl_ids[gene_name].append(ensembl_id)
                    # print('\t added gene {}'.format(ensembl_id))

                elif feature_type in self.transcript_types:
#                     print('\ttranscript has gene parent {}. {}'.format(parent, parent in genes))
                    # try:
                    #     transcript_support_level = int(fields['transcript_support_level'].split(' ')[0])
                    # except ValueError:
                    #     passed_support_filter = False
                    # else:
                    #     passed_support_filter = transcript_support_level < self.maximum_transcript_support
                    ensembl_id = fields['ID']
                    transcripts[ensembl_id] = {'contig': contig,
                                               'start': start - 1,  # convert 1-based to 0-based
                                               'end': end,
                                               'strand': strand,
                                               'components': []}
                    transcripts[ensembl_id].update(fields)

                    # print('\t added transcript {} with parent {}'.format(ensembl_id, parent))


                elif feature_type in self.component_types:
#                     print('\tcomponent has transcript parent {}. {}'.format(parent, parent in transcripts))
                    if 'exon_id' in fields:
                        ensembl_id = fields['exon_id']
                    else:
                        ensembl_id = str(component_num)
                        component_num += 1

                    components[ensembl_id] = {'contig': contig,
                                              'start': start - 1,  # convert 1-based to 0-based
                                              'end': end,
                                              'strand': strand,
                                              'type': feature_type}
                    components[ensembl_id].update(fields)
                    
        for transcript_id, transcript_data in transcripts.items():
            parent = transcript_data['Parent']
            if parent in genes:
                genes[parent]['transcripts'].append(transcript_id)
            else:
                print('orphan transcript {} with missing parent {}!'.format(transcript_id, parent))
                
        for component_id, component_data in components.items():
            parent = component_data['Parent']
            if parent in transcripts:
                transcripts[parent]['components'].append(component_id)
            else:
                print('orphan component {} with missing parent {}!'.format(component_id, parent))

        return genes, transcripts, components, gene_names_to_ensembl_ids


        
class _MatrixData:
    def __init__(self):
        pass
    
    def _query(self):
        print('Must be overridden by inheritors')
        
    def query(self, chrom, start, end):
        return self._query(query_chrom=chrom, query_start=start, query_end=end)

        
class HicDataDir(_MatrixData):
    def __init__(self, 
                 fname_template='/home/dskola/projects/coupled_peaks/hic/c57_hic_corrs_{}.tsv',
                 binsize=10000):
        self.fname_template = fname_template
        self.binsize = binsize
        
    def _query(self, query_chrom, query_start, query_end):
        this_chrom_fname = self.fname_template.format(query_chrom)
        this_chrom_data = pandas.read_csv(this_chrom_fname, sep='\t', index_col=0)
        
        rounded_start = utilities.roundto(query_start, binsize)
        rounded_end = utilities.roundto(query_end, binsize)

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
        hic_df.index = [int(name.split('-')[1]) for name in hic_df.index]
        hic_df.columns = [int(name.split('-')[1]) for name in hic_df.columns] 
        return hic_df