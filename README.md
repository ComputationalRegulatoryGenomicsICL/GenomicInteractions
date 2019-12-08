# GenomicInteractions

*GenomicInteractions* is an R/Bioconductor package for manipulating and investigating chromatin interaction data. It provides simple annotation of genomic features which interfaces with existing Bioconductor packages, as well as numerous plotting functions and summary statistics.  

## Installation

To install this package, start R and enter:

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicInteractions")
```

To install the latest version from this repository using the devtools package:

```
library(devtools)
install_github("ComputationalRegulatoryGenomicsICL/GenomicInteractions")
```

## Bioconductor links and documentation

The [Bioconductor page](http://www.bioconductor.org/packages/release/bioc/html/GenomicInteractions.html) contains documentation for the current release.

To view example Hi-C and ChIA-PET analysis for using this package on your system, start R and enter:

```
browseVignettes("GenomicInteractions")
```

## Citation 

From within R, enter `citation("GenomicInteractions")`:

```
  Harmston, N., Ing-Simmons, E., Perry, M., Baresic, A. and
  Lenhard, B. (2015). GenomicInteractions: An R/Bioconductor
  package for manipulating and investigating chromatin interaction
  data. BMC Genomics 16:963.
  
A BibTeX entry for LaTeX users is

  @Article{,
    title = {GenomicInteractions: An R/Bioconductor package for manipulating and investigating chromatin interaction data},
    author = {Nathan Harmston and Elizabeth Ing-Simmons and Malcolm Perry and Anja Bare{\v{s}}i{\'{c}} and Boris Lenhard},
    year = {2015},
    journal = {BMC Genomics},
    volume = {16},
    issue = {963},
    doi = {10.1186/s12864-015-2140-x},
    url = {https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2140-x},
  }
```

## Contributing

We welcome contributions in any form: suggestions, issues, bugfixes. Pull requests should be made to the development branch.

