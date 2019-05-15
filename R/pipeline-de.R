
##library(tximport)
##library(GenomicFeatures)
##library(DESeq2)
##library(AnnotationDbi)


#' Run deseq2
#' 
#' @param pipeline A pipeline object
#' @param ... other arguments
#' @export
deseq2 <- function (pipeline, ...)
{
    UseMethod("deseq2", pipeline)
}

#' Run salmon for genes expression
#' 
#' @param pipeline A pipeline object
#' @rdname deseq2
#' @export
deseq2 <- function (pipeline, design=~condition)
{
    last.res <- getResultsByClass(pipeline, "rnaseq_quantification")
    if (is.null(last.res))
    {
        warning("Have you called salmon()? Doing nothing! No deseq2 call!")
        return(pipeline)
    }

    logging("Running DESeqDataSetFromMatrix()", .module="deseq2")
    metadata <- getFilter(pipeline)
    inter <- intersect(colnames(last.res$counts), row.names(metadata))
    counts <- apply(last.res$counts, c(1, 2), as.integer)
    dds <- DESeq2::DESeqDataSetFromMatrix(
                       counts,
                       metadata[inter,, drop=FALSE],
                       design)

    ## TODO: custom filter because too low expression
    dds <- dds[ apply(DESeq2::counts(dds), 1,
                      function (x) sum(x > 5) > ncol(dds)/2), ]
    logging("Running estimateSizeFactors()", .module="deseq2")
    ## TODO: estimateSizeFactorsForMatrix. Necessary ??!?
    dds <- DESeq2::estimateSizeFactors(dds)
    ## logging("Running estimateSizeFactorsForMatrix()", .module="deseq2")
    ## dds <- DESeq2::estimateSizeFactorsForMatrix(
    ##                    counts(dds)/
    ##                    last.res$effective_length[, inter, drop=FALSE])

    logging("Running DESeq()", .module="deseq2")
    dds <- DESeq2::DESeq(dds, parallel=T)

    pipeline$results <- append(list(dds), pipeline$results)
    
    return(pipeline)
}


#' Return the results
#' 
#' @param pipeline A pipeline object
#' @param c1 Condition 1
#' @param c2 Condition 2
#' @param condition The column name of sampleTable passed to newSamples()
#'     to use as contrast
#' @param isoforms True if isoforms result (default: NULL)
#' @export
deseq2Results  <- function (pipeline, c1, c2, condition="condition", isoforms=NULL, ...)
{
    dds <- getResultsByClass(pipeline, "DESeqDataSet")
    isoforms <- startsWith(row.names(counts(dds))[1], "ENST")
    res <- DESeq2::results(dds, contrast=c(condition, c1, c2), ...)
    if (isoforms)
    {
        pipeline <- tx2genes(pipeline)

        ##rownames(res)  <- sapply(strsplit(rownames(res),'\\.'),'[',1)
        res <- merge(as.data.frame(res), pipeline$tx2gene,
                     by.x=0, by.y="TXNAME", all.x=T)
        rownames(res) <- res$Row.names
        res <- res[-1]
    }
    return (res)
}

## pipeline() %>% metadata("data/metadata.tsv") %>%
##     options(
##         salmon.index="~/data/ensembl/GRCh38/Homo_sapiens.GRCh38.cdna.all.salmon.index",
##         gtf="~/data/ensembl/GRCh38/Homo_sapiens.GRCh38.95.gtf.gz") %>%
##     salmonGenes()
