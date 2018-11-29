
.valid.counts <- function(object){
    counts <- assays(object)$counts
    if(length(dim(counts))!=3){
        return("The counts assays must have 3 dimensions")
    }
    NULL
}

.valid.design <- function(object){
    design <- colData(object)
    cns <- colnames(design)
    if(!all(c("sample", "group") %in% cns)){
        return("design dataframe must have sample and group columns")
    }
    if(!all(table(design$sample)==2) | length(unique(design$group))!=2){
        return("design must be paired")
    }
    NULL
}

.valid.lengths <- function(object){
    lengths <- rowData(object)
    cns <- colnames(lengths)
    if(!all(c("iLen", "sLen") %in% cns)){
        return("lengths dataframe must have iLen and sLen columns")
    }
    if(!is.integer(as.matrix(lengths))){
        return("lengths must be integer")
    }
    NULL
}

.valid.PDseDataSet <- function(object){
    c(.valid.counts(object),
      .valid.design(object),
      .valid.lengths(object))
}

#' @rdname PDseDataSet
#' @import SummarizedExperiment
#' @export
setClass("PDseDataSet",
         contains = "SummarizedExperiment",
         validity = .valid.PDseDataSet)

#' PDseDataSet object and constuctor
#'
#' `PDseDataSet` is a subclass of `SummarizedExperiment`. It can used to store inclusion and skipping splicing counts for pair designed samples.
#'
#' @param counts The counts of splicing events, including inclusion and skipping counts in 3 dimensions for each sample.
#' @param design The paired design data.frame, including sample column for sample ids and group column for design factors.
#' @param lengths Two columns iLen and sLen for the effective lengths of inclusion and skipping isoforms.
#' @return A PDseDataSet object
#' @rdname PDseDataSet
#' @importFrom methods new
#' @export
PDseDataSet <- function(counts, design, lengths){
    se <- SummarizedExperiment(assays = list(counts = counts),
                               rowData = DataFrame(lengths),
                               colData = DataFrame(design))
    new("PDseDataSet", se)
}

.appNA <- function(datL){
    ssize <- max(lengths(datL[[1]]))
    for(ct in c("I1", "S1", "I2", "S2")){
        datL[[ct]] <- lapply(datL[[ct]], function(x)c(x, rep(NA, ssize - length(x))))
    }
    return(datL)
}

#' PDseDataSet from rMATs/PAIRADISE Mat format
#'
#' The Mat format should have 7 columns, arranged as follows:
#' Column 1 contains the ID of the alternative splicing events.
#' Column 2 contains counts of isoform 1 corresponding to the first group.
#' Column 3 contains counts of isoform 2 corresponding to the first group.
#' Column 4 contains counts of isoform 1 corresponding to the second group.
#' Column 5 contains counts of isoform 2 corresponding to the second group.
#' Column 6 contains the effective length of isoform 1.
#' Column 7 contains the effective length of isoform 2.
#' Replicates in columns 2-5 should be separated by commas, e.g. "1623,432,6" for three replicates and the replicate order should be consistent for each column to ensure pairs are matched correctly.
#' @param dat The Mat format dataframe.
#' @return A PDseDataSet object
#' @import abind
#' @rdname PDseDataSetFromMat
#' @export
PDseDataSetFromMat <- function(dat){
    datL <- clean.data(dat)
    datL <- .appNA(datL)
    lens <- data.frame(iLen = as.integer(datL$length_I),
                       sLen = as.integer(datL$length_S))

    iCounts <- cbind(do.call(rbind, datL$I1),
                     do.call(rbind, datL$I2))
    sCounts <- cbind(do.call(rbind, datL$S1),
                     do.call(rbind, datL$S2))
    counts <- abind(iCounts, sCounts, along = 3)
    ##design
    design <- DataFrame(sample=rep(paste0("S", seq(max(datL$M))), 2),
                        group=rep(c("T", "N"), each=max(datL$M)))
    ##counts
    ids <- paste(design$sample, design$group, sep=".")
    rownames(counts) <- rownames(lens) <- datL$exonList
    colnames(counts) <- ids
    
    PDseDataSet(counts, design, lens)
}

#' PDseDataSet counts
#' @param object A PDseDataSet object
#' @return A counts matrix 
counts <- function(object){
    stopifnot(is(object, "PDseDataSet"))
    assays(object)$counts
}

setGeneric("counts")
