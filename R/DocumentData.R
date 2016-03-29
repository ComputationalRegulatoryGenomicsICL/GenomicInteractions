#' Example HiC dataset
#' 
#' This dataset contains HiC data from Seitan et al. 2013. The data was 
#' analysed using HOMER (Heinz et al. 2010) at a resolution of 100kb to 
#' find significant interactions. This example dataset has been filtered
#' to retain only interactions on chromosomes 14 and 15 with a FDR < 0.1.
#' The data has also been annotated for overlaps with Refseq promoters.
#' See the HiC analysis vignette (\code{vignettes(GenomicInteractions)}) for 
#' more information on how this dataset was created.
#' 
#' @name hic_example_data
#' @docType data
#' @keywords datasets
#' @usage data(hic_example_data)
#' @format A GenomicInteractions object with length 8171.
#' @return GenomicInteractions object
#' @references Seitan, V. C. et al. Cohesin-based chromatin interactions enable 
#' regulated gene expression within pre-existing architectural compartments. 
#' Genome Res. 23, 2066-77 (2013).
#' 
#' Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of 
#' Lineage-Determining Transcription Factors Prime cis-Regulatory Elements 
#' Required for Macrophage and B Cell Identities. Mol Cell 2010 
#' May 28;38(4):576-589.
NULL

#' Mouse Refseq promoters from chr 14-15
#' 
#' This dataset contains a subset of the promoters from the Refseq annotation
#' for mouse genome build mm9. See the HiC analysis vignette 
#' (\code{vignettes(GenomicInteractions)}) for more information on how this 
#' dataset was created.
#' 
#' @name mm9_refseq_promoters
#' @docType data
#' @keywords datasets
#' @usage data(mm9_refseq_promoters)
#' @format A GRanges object with length 2441.
#' @format A GRanges object
#' @return A GRanges object.
NULL

#' Human Refseq transcripts from chr 17-18
#' 
#' This dataset contains a subset of the transcripts from the Refseq annotation
#' for mouse genome build hg19 See the ChIA-PET analysis vignette 
#' (\code{vignettes(GenomicInteractions)}) for more information on how this 
#' dataset was created.
#' 
#' @name hg19.refseq.transcripts
#' @docType data
#' @keywords datasets
#' @usage data(hg19.refseq.transcripts)
#' @format A GRanges object with length 2441.
#' @format A GRanges object
#' @return A GRanges object.
NULL

#' Putative enhancers from mouse thymus data
#' 
#' This dataset contains a set of mouse thymus enhancers derived from ChIP-seq data
#' from mouse thymus, as described in Shen et al. 2012. See the HiC analysis vignette 
#' for more details. (\code{vignettes(GenomicInteractions)}). 
#' 
#' @name thymus_enh
#' @docType data
#' @keywords datasets
#' @usage data("thymus_enhancers")
#' @format A GRanges object
#' @return A GRanges object
#' @references Shen, Y et al. A map of cis-regulatory sequences in the mouse genome. Nature (2012).
NULL
