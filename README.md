##GenomicInteractions

*GenomicInteractions* is an R/Bioconductor package for manipulating and investigating chromatin interaction data. It provides simple annotation of genomic features which interfaces with existing Bioconductor packages, as well as numerous plotting functions and summary statistics.  

###Installation

To install this package, start R and enter:

```
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicInteractions")
```

To install the latest version from this repository using the devtools package:

```
library(devtools)
install_github("ComputationalRegulatoryGenomicsICL/GenomicInteractions")
```

###Bioconductor links and documentation

The bioconductor page contains documentation for the released package: [Bioconductior Release 3.0](http://www.bioconductor.org/packages/release/bioc/html/GenomicInteractions.html)

To view example Hi-C and ChIA-PET analysis for using this package on your system, start R and enter:

```
browseVignettes("GenomicInteractions")
```

###Citation 

From within R, enter `citation("GenomicInteractions")`:

```
Harmston, N., Ing-Simmons, E., Perry, M., Baresic, A., and Lenhard, B. (2014). GenomicInteractions: R package for handling genomic interaction data. R package version 1.0.0.
```

###Contributing

We welcome contributions in any form: suggestions, issues, bugfixes. Pull requests should be made to the development branch.

