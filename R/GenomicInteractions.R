## Definition of GenomicInteractions

#' A S4 class to represent interactions between genomic regions.
#' 
#'  @slot experiment_name Character. Experiment name.
#'  @slot description Character. Longer description of experiment. 
#'  @slot genome_name Character. Genome version for experiment data, should correspond to a BSgenome data package.
#'  @slot anchor_one,anchor_two GRanges. Set of anchors of interactions.
#'  @slot counts Numeric. Counts of reads supporting each interaction.
#'  @slot normalised_counts Numeric. Normalised counts of reads supporting each interaction.
#'  @slot pvalue Numeric. P-values for individual interactions.
#'  @slot fdr Numeric. FDRs for individual interactions. 
#'  
#' This class is used to store information on which genomic regions are interacting with each other, the number of counts supporting each interaction,
#' and associated p-values and false-discovery rates (FDR). Objects of this class contain information of the genomic coordinates of the interacting
#' regions and the strength of these interactions, and associated metadata such as the name of the dataset and a brief description of the dataset. 
#' Interacting regions are stored as a pair of GenomicRanges: each set of anchor regions is stored as a separate GenomicRanges object, accessed by 
#' \code{getAnchorOne} and \code{getAnchorTwo}. 
#' 
#' @examples
#' 
#' showClass("GenomicInteractions")
#' 
#' library(BSgenome.Mmusculus.UCSC.mm9)
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5), seqlengths=seqlengths(Mmusculus))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5), seqlengths=seqlengths(Mmusculus))
#' test <- new("GenomicInteractions", experiment_name="test", description="this is a test", 
#'                  genome_name="BSgenome.Mmusculus.UCSC.mm9", anchor_one = anchor.one, 
#'                  anchor_two = anchor.two, counts=as.integer(c(2,1,2,3)), pvalue=c(0.1, 0.3, 0.1, 0.08))
#'  
#' @import GenomicRanges
#' @export GenomicInteractions
setClass("GenomicInteractions", 
    representation(experiment_name = "character", 
                   description = "character", 
                   genome_name = "character", 
                   anchor_one ="GRanges", 
                   anchor_two="GRanges", 
                   counts="integer", 
                   normalised_counts = "numeric", 
                   pvalue="numeric", 
                   fdr="numeric"),
    prototype(experiment_name = character(), 
              description = character(), 
              genome_name = character(), 
              anchor_one = GRanges(), 
              anchor_two = GRanges(), 
              counts = integer(), 
              normalised_counts = numeric(), 
              pvalue = numeric(), 
              fdr = numeric() ),
    validity = function(object){ 
        if(!(object@genome_name %in% suppressWarnings(suppressMessages(BSgenome::available.genomes())))){
            return("'genome_name' must be a name of one of the genome packages available in BSgenome! See 'available.genomes()'")
        }else if(length(object@anchor_one) == 0 ){
            return("anchor one cannot be of length 0")
        }else if(length(object@anchor_two) == 0 ){
            return("anchor two cannot be of length 0")
        }else if(length(object@anchor_one) != length(object@anchor_two)){
            return("length of anchor one and anchor two do not match")
        }else if(length(object@counts) == 0 & length(object@normalised_counts) == 0){
            return("counts or normalised counts must contain information")
        }else if(object@experiment_name == ""){
            return("experiment_name cannot be empty")
        }else if( (length(object@anchor_one) != length(object@counts)) & length(object@counts ) > 0 ){
            return("Anchor length and length of counts do not match")
        }else if( (length(object@anchor_one) != length(object@normalised_counts)) & length(object@normalised_counts ) > 0 ){
            return("Anchor length and length of normalised counts do not match")
        }else if( (length(object@anchor_one) != length(object@pvalue)) & length(object@pvalue) > 0 ){
            return("Anchor length and length of pvalue do not match") 
        }else if( (length(object@anchor_one) != length(object@fdr)) & length(object@fdr) > 0 ){
            return("Anchor length and length of FDR do not match")
        }else{
            return(TRUE)}}
)

#' Function to create GenomicInteraction objects.
#' 
#' Function to create GenomicInteraction objects from a variety of files. The resulting objects contain information
#' on which genomic regions are interacting with each other, and the number of counts supporting each interaction. 
#' It is also possible to store information on associated p-values and false-discovery rates (FDR). 
#' It is possible to create GenomicInteractions objects for various datasets including Hi-C and ChIA-PET. It is possible 
#' to read interactions from a variety of files including BAM files, bed files (BED12 and BEDPE) and from the output 
#' from standard processing pipelines, such as HOMER and ChIA-PET tool. GenomicInteractions objects can also be created 
#' using calls of the form \code{new("GenomicInteractions", ...)}. For hiclib, it expects the directory in which the files 
#' extracted using h5dictToTxt.py from the hdf5 file are located, where as for all of the other file types it expects the full
#' filename. 
#' 
#' @param fn Filename or, if type="hiclib", folder
#' @param type One of "chiapet.tool", "chiapet.encode", "bed12", "bedpe", "hiclib", "homer", "bam". 
#' @param experiment_name Experiment name.
#' @param description Description of experiment.
#' @param gname Genome name to use for constructing the GenomicInteractions object.
#' @return a GenomicInteractions object
#' 
#' @importFrom Rsamtools scanBamFlag ScanBamParam scanBam bamFlagAsBitMatrix
#' @importFrom IRanges IRanges
#' @import data.table
#' @importFrom stringr str_split
#' @importFrom rtracklayer import.bed
#' 
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' 
#' k562.rep1 = GenomicInteractions(file.path(system.file(package="GenomicInteractions"), "extdata", "k562.rep1.cluster.pet3+.txt"), 
#'                                      type="chiapet.tool", 
#'                                      experiment_name="k562", 
#'                                      description="k562 pol2 8wg16", 
#'                                      gname="BSgenome.Hsapiens.UCSC.hg19")
#'                                      
#' k562.rep1
#' 
#' @export

GenomicInteractions = function(fn, type, experiment_name, description, gname){
    genome = .loadGenome(gname)
	if(type == "chiapet.tool"){
	    dat = read.table(fn, header=TRUE, stringsAsFactors=FALSE, sep="\t")
        anchor1 = GRanges(dat[,"chrom_left"], 
                            IRanges(as.integer(dat[,"start_left"])+1, as.integer(dat[,"end_left"])), 
                            seqlengths=seqlengths(genome))
        anchor2 = GRanges(dat[,"chrom_right"], 
                            IRanges(as.integer(dat[,"start_right"])+1, as.integer(dat[,"end_right"])), 
                            seqlengths=seqlengths(genome))
        counts = as.integer(dat[,"pet.counts.between.left.and.right.anchors"])
        p.value = as.numeric( dat[, "p.value"])
        fdr = as.numeric( dat[, "FDR"])
    }else if(type == "chiapet.encode"){
        dat = .processChiapetName(unique(import.bed(fn)$name))
        anchor1 = GRanges(dat[,"chrom.left."], 
                          IRanges(as.integer(dat[,"start.left."]), as.integer(dat[,"end.left."])), 
                          seqlengths=seqlengths(genome))
        anchor2 = GRanges(dat[,"chrom.right."], 
                          IRanges(as.integer(dat[,"start.right."]), as.integer(dat[,"end.right."])), 
                          seqlengths=seqlengths(genome))
        counts = as.integer(dat[,"counts"])
        p.value = numeric(0) 
        fdr = numeric(0) 
    }else if(type == "bed12"){
        bedfile = import.bed(fn)
        dat = .processChiapetName(unique(bedfile$name))
        anchor1 = GRanges(dat[,"chrom.left."], 
                          IRanges(as.integer(dat[,"start.left."]), as.integer(dat[,"end.left."])), 
                          seqlengths=seqlengths(genome))
        anchor2 = GRanges(dat[,"chrom.right."], 
                          IRanges(as.integer(dat[,"start.right."]), as.integer(dat[,"end.right."])), 
                          seqlengths=seqlengths(genome))
        counts = as.integer(dat[,"counts"])
        p.value = numeric(0)
        fdr = numeric(0) 
    }else if(type == "bedpe"){
        dat = read.table(fn, stringsAsFactors=FALSE, sep="\t")
        anchor1 = GRanges(dat[,1], 
                          IRanges(dat[,2]+1, dat[,3]), strand=ifelse(ncol(dat) >= 10, dat[,9], "*"), 
                          seqlengths=seqlengths(genome))
        anchor2 = GRanges(dat[,4], 
                          IRanges(dat[,5]+1, dat[,6]), strand=ifelse(ncol(dat) >= 10, dat[,10], "*"), 
                          seqlengths=seqlengths(genome))
        counts = as.integer(rep(1, length(anchor1))) 
        p.value = numeric(0)
        fdr = numeric(0)
    }else if(type == "hiclib"){ 
	    dat = .importHicLib(fn, genome)
        anchor1 = GRanges(dat$chrm1, 
                          IRanges(dat$mid1-round(dat$fraglength1/2), dat$mid1 + round(dat$fraglength1/2)), 
                          seqlengths=seqlengths(genome), fragid=dat$fragid1)
		anchor2 = GRanges(dat$chrm2, 
                          IRanges(dat$mid2-round(dat$fraglength2/2), dat$mid2 + round(dat$fraglength2/2)), 
                          seqlengths=seqlengths(genome), fragid=dat$fragid2)
        counts = dat$N
		p.value = numeric(0) 
		fdr = numeric(0)
    }else if(type == "homer"){
        dat = .importHomer(fn)
        anchor1 = GRanges(dat$chr.1., 
                          IRanges(dat$start.1., dat$end.1.), 
                          seqlengths=seqlengths(genome))
        anchor2 = GRanges(dat$chr.2., 
                          IRanges(dat$start.2., dat$end.2.), 
                          seqlengths=seqlengths(genome))
        counts = as.integer(dat$Interaction.Reads)
        p.value = exp(as.numeric(dat$LogP))
        fdr = as.numeric(dat$FDR.Benjamini.)
	}else if(type == "bam"){
        dat = .readBam(fn, genome)
        anchor1 = dat[[1]]
        anchor2 = dat[[2]]
        counts = as.integer(rep(1, length(anchor1))) 
        p.value = numeric(0)
        fdr = numeric(0)
	}else{
        stop("type is not one of \"chiapet.tool\", \"chiapet.encode\", \"bed12\", \"bedpe\", \"hiclib\", \"homer\", \"bam\"")
	}
    giobject = new("GenomicInteractions", 
                 experiment_name = experiment_name, 
                 description = description, 
                 genome_name = gname, 
                 anchor_one=anchor1, 
                 anchor_two=anchor2, 
                 counts=counts, 
                 pvalue=p.value, 
                 fdr=fdr)
    return(giobject)
}
