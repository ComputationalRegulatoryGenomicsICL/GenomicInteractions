
#' Summarise Interactions between defined anchors
#'
#' Calculate the number of of paired-end reads mapping between a defined set of anchors.
#' This function will ignore counts present in the input data.
#'
#' @return A GenomicInteractions object with annotated counts between anchors
#' @docType methods
#' @rdname countsBetweenAnchors-methods
#' @export
setGeneric("countsBetweenAnchors",function(x, y, ...){standardGeneric ("countsBetweenAnchors")})

#' @param x A GenomicInteractions object
#' @param y A GenomicRanges object
#' @param ignore_overlaps Allow overlapping anchors. Use this when you have overlapping anchors
#'                        but be careful with multi-mapping. The "within" option can help with this.
#' @param ... Extra parameters to pass to findOverlaps
#' @import GenomicRanges
#' @rdname countsBetweenAnchors-methods
#' @docType methods
#' @export
setMethod("countsBetweenAnchors", list("GenomicInteractions", "GRanges"), function(x, y, ignore_overlaps=FALSE, ...) {
    #check anchors are unique
    if (ignore_overlaps == FALSE && any(countOverlaps(y, y) > 1)) stop("anchors are not unique")
    one = overlapsAny(anchorOne(x), y, ...)
    two = overlapsAny(anchorTwo(x), y, ...)
    x.valid = x[one & two]
    overlaps = findOverlaps(sort(x.valid), y, select="first", ...) # select produces matrix not Hits
    interactions = paste(overlaps[[1]], overlaps[[2]], sep=":")
    tabulated = table(interactions)

    pairs_list = strsplit(names(tabulated), ":")
    pairs_one = as.integer(sapply(pairs_list, function(x) x[1]))
    pairs_two = as.integer(sapply(pairs_list, function(x) x[2]))

    anchor_one = y[pairs_one]
    anchor_two = y[pairs_two]
    counts = as.integer(tabulated)

    final_counts = GenomicInteractions(
                       anchor_one=anchor_one,
                       anchor_two=anchor_two,
                       experiment_name = name(x),
                       description = description(x),
                       counts=counts)

    return(sort(final_counts))
})

#' Remove all but one occurences of a duplicated interaction
#'
#' Removes all but the first occurence of a duplicated interaction (defined as
#' having identical coordinates for both anchors). N.B. this does not summarise
#' the total counts of all the duplicates. It is designed for removing potential
#' PCR duplicates after reading in .bam files.
#'
#' @param GIObject A GenomicInteractions object.
#' @return A GenomicInteractions object that is a subset of the input object.
#' @import GenomicRanges
#' @export

removeDups <- function(GIObject){
    dat <- data.frame(Chr1 = seqnames(anchorOne(GIObject)),
                      Start1 = start(anchorOne(GIObject)),
                      Chr2 = seqnames(anchorTwo(GIObject)),
                      Start2 = start(anchorTwo(GIObject))
    )
    idx <- which(!duplicated(dat))
    reads_removed <- length(GIObject) - length(idx)
    percent_removed <- signif(100*reads_removed / length(GIObject), 3)
    message(paste0("Removing ", reads_removed, " duplicate PETs (", percent_removed, "%)"))
    return(GIObject[idx])
}

#' Tests whether anchors have the same strand.
#'
#' This is designed for processing .bam files.
#'
#' @param GIObject A GenomicInteractions object
#' @return A logical vector denoting with TRUE if both anchors of an interaction
#'  are on the same strand and FALSE otherwise.

sameStrand <- function(GIObject){
    return(strand(anchorOne(GIObject))==strand(anchorTwo(GIObject)))
}

#' Get self ligation threshold with SD method from Heidari et al
#' 
#' This function calculates a self ligation threshold according to 
#' the method published in Heidari et al., Genome Research, 2014. 
#' Briefly, paired reads are divided into in evenly sized bins. For
#' each bin, the log2 ratio of reads that are aligned to opposite strand
#' vs to the same strand is calculated. Twice the standard deviation of 
#' this ratio at high distances is used a cutoff to determine which bins 
#' are likely to contain mostly self-liagted reads.
#' 
#' @param GIObject a GenomicInteractions object of paired end reads
#' @param bins Number of evenly sized bins to use.
#' @param distance_th The threshold, in base pairs, to use as a cutoff to 
#' pick which bins to use to determine the standard deviation.
#' @param plot TRUE by default. Whether to plot the log2ratio of opposite
#' to same strand reads vs distance.
#'
#' @import dplyr
#' @import ggplot2
#' @export
#' @return The cutoff in base pairs below which an interaction is likely to be a self ligation.

get_self_ligation_threshold <- function(GIObject, bins=100, distance_th=400000, plot=TRUE){
    #get df
    stranded_df <- data.frame(Distance=calculateDistances(GIObject), SameStrand=sameStrand(GIObject))
    stranded_cis_df <- stranded_df[complete.cases(stranded_df),]
    stranded_cis_df <- stranded_cis_df[order(stranded_cis_df$Distance),]

    #bin data
    bin_n <- nrow(stranded_cis_df )/bins

    cuts <- 1:bins * bin_n
    breaks <- stranded_cis_df$Distance[cuts]

    stranded_cis_df$Bin <- as.numeric(as.character(cut(stranded_cis_df$Distance, breaks=breaks, labels=breaks[1:length(breaks)-1], include.lowest = TRUE)))
    byBin <- group_by(stranded_cis_df, Bin)

    #summarise by bin
    sum_byBin <- summarise(byBin, Total=n(), SameStrand=sum(SameStrand))
    sum_byBin <- mutate(sum_byBin, OppStrand=Total-SameStrand,
                        log2Ratio=log2((OppStrand+1)/(SameStrand+1))) #pseudocount to avoid NaN errors
    sum_byBin <- mutate(sum_byBin, OppPercent=100*OppStrand/Total, SamePercent=100*SameStrand/Total)

    #get cutoff of log2ratio
    sum_byBin %>%
        filter(Bin > distance_th) %>%
        select(log2Ratio) %>%
        unlist() %>%
        mean() -> longrange_mean_log2

    sum_byBin %>%
        filter(Bin > distance_th) %>%
        select(log2Ratio) %>%
        unlist() %>%
        sd() -> longrange_sd_log2

    lower <- longrange_mean_log2 - 2*longrange_sd_log2
    upper <- longrange_mean_log2 + 2*longrange_sd_log2
    bp_cutoff <- min(sum_byBin[sum_byBin$log2Ratio > lower & sum_byBin$log2Ratio < upper,"Bin"])

    if (plot){
        print(ggplot(sum_byBin, aes(x=Bin, y=log2Ratio)) + geom_line() + geom_point() +
                  geom_hline(aes_string(yintercept=lower)) +
                  geom_hline(aes_string(yintercept=upper)) +
                  coord_cartesian(xlim=c(0, 20000)) + geom_vline(xintercept=bp_cutoff, linetype="dashed") +
                  xlab("Distance (bp)") + ylab("log2 ratio opposite strand pairs / same strand pairs")
        )
    }
    return(bp_cutoff)
}

#' get self ligation threshold with binomial test
#' 
#' This function calculates a self ligation threshold according to 
#' a method based on that of Heidari et al., Genome Research, 2014. 
#' Briefly, paired reads are divided into in evenly spaced bins. For
#' each bin, the number of reads that are aligned to opposite strand
#' vs to the same strand is calculated. A binomial test is used to test
#' if this is significantly different from the 50:50 ratio expected by 
#' chance if all reads are real interactions. 
#' 
#' @param GIObject a GenomicInteractions object of paired end reads
#' @param bin.size Bin size in base pairs.
#' @param max.distance The maximum distance to consider between reads. 
#' Reads further apart than this distance should be very unlikely to be
#'  self ligations.
#' @param p.cutoff P value cut off for a significant difference from 50:50. Default: 0.05
#' @param adjust Method to use to adjust p values. Default: fdr. See `help(p.adjust)` for 
#' accepted values. Can also be NA for no adjustment.
#' @param plot TRUE by default. Whether to plot the percentage of reads 
#' on opposite strands vs difference and the binomial test p value vs distance.
#'
#' @export
#' @return The cutoff in base pairs below which an interaction is likely to be a self ligation.


get_binom_ligation_threshold = function(GIObject, max.distance=20000, bin.size=500, p.cutoff=0.05, adjust="fdr", plot=TRUE){

    #make data frame
    stranded_df <- data.frame(Distance=calculateDistances(GIObject), SameStrand=sameStrand(GIObject))
    stranded_cis_df <- stranded_df[complete.cases(stranded_df),]
    stranded_cis_df <- stranded_cis_df[order(stranded_cis_df$Distance),]
    stranded_cis_df = stranded_cis_df[ stranded_cis_df$Distance < max.distance, ]

    #bin data
    bins = cut(stranded_cis_df$Distance, breaks=seq(0, max.distance, by=bin.size), include.lowest=TRUE)
    stranded_cis_df$Bin = bins
    byBin <- group_by(stranded_cis_df, Bin)
    sum_byBin <- summarise(byBin, Total=n(), SameStrand=sum(SameStrand))

    #get and adjust p values
    sum_byBin$p.value <- sapply(1:nrow(sum_byBin), function(x){binom.test(sum_byBin$SameStrand[x], sum_byBin$Total[x])$p.value})

    if(!is.na(adjust)){
        sum_byBin$p.value <- p.adjust(sum_byBin$p.value, method=adjust)
    }

    #get cutoff
    bp_cutoff <- seq(0, max.distance, by=bin.size)[min(which(sum_byBin$p.value > p.cutoff))]

    if (plot){
        #data for plotting
        sum_byBin <- mutate(sum_byBin, OppStrand=Total-SameStrand,
                            Bin=bin.size*as.numeric((Bin)),
                            OppPercent=100*OppStrand/Total)

        #plot % opposite strand reads and cutoff
        print(ggplot(sum_byBin, aes(x=Bin, y=OppPercent)) + geom_line() + geom_point() +
                  coord_cartesian(xlim=c(0, max.distance)) + geom_vline(xintercept=bp_cutoff, linetype="dashed") +
                  ylab("Opposite Strand Percentage")
        )
        #plot p values, p value cutoff and distance cutoff
        print(ggplot(sum_byBin, aes(x=Bin, y=p.value)) + geom_line() + geom_point() +
                  coord_cartesian(xlim=c(0, max.distance)) + geom_hline(xintercept=p.cutoff, linetype="dashed") +
                  geom_vline(xintercept=bp_cutoff, linetype="dashed") +
                  ylab("p value") + xlab("Distance (bp)")
        )
    }
    #return distance cutoff
    return(bp_cutoff)
}


