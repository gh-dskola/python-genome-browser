# ToDo:

1. Gene Models
    1. Ensure compatibility with all versions of Ensembl GFF3
    2. Add GTF support
    3. Add RefSeq support
    4. Auto-scale heights of gene models
2. Vector data
    1. Add support for stranded data
    2. Add support for text Wig / bedGraph files
    3. Add BAM support
    4. `TagDirectory`: include reads that start left of the window and overhang into it
3. Interaction data
    1. Add arbitrary interactions (not confined to bins). Plot "arbitrary bins" (start-end)
    2. Ensure current parser is compatible with BED-PE
4. Matrix data
    1. Use tabix to search CSV matrices on disk and read only applicable rows
    2. Plot diagonal cells
        * Generate image, then transform, instead of transforming data and generating image. We can
        probably do this transparently by passing a Transform object to the axes.
5. General
    1. Shrink logo.
    2. Make HiC bins diagonal instead of square.
    3. Add parameter for subplots to share yaxis limits.
    4. Allow for querying entire chromosome by passing 0s as start and end arguments
    5. Add automatic windowing function to reduce complexity of plots (and corresponding PDF file sizes).
