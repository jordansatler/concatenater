# Concatenater

This script will concatenate a set of loci. The loci need to be in fasta format (files ending in .fa, .fas, .fasta) and are not to be interleaved. By default, all loci are concatenated. If you would like to filter loci that have some minimum number of individuals sampled, you can add this threshold in the concatenate_loci function. It is currently commented out but you can uncomment and set your sampling threshold (as a percent of missing taxa).

usage:  
```python
    python concatenater.py /path/to/fasta/files
```

Output will be two files. One file will be the concatenated data set in phylip format, and the other file will show the locus partitions in the concatenated file.
