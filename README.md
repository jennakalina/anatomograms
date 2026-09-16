# anatomograms

All code and data associated with Drosophila anatomograms.

## Directories & References

**source_info.csv** contains source data information on each anatomogram, including both the paper and the URL to get the data that was used.

**scripts** contains the scripts used to map paper annotations to anatomogram terms, subset by cell type, and calculate various average expression, percent expression, and CPM counts.

**processedData** contains zips of all files used by the anatomogram web tool. There is one CSV file for each data source that contains all genes and their expression metrics for each graph and anatomical term. Expanding the zip shows the different subdirectories for different anatomograms (if a data source maps to multiple graphs) and each file is named by the anatomical term it maps to.
