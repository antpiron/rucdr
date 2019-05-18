
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

    txi  <- last.res$txi
    
    logging("Running DESeqDataSetFromMatrix()", .module="deseq2")
    metadata <- getFilter(pipeline)
    inter <- intersect(colnames(Counts(txi)), row.names(metadata))

    txi.filtered <- lapply(txi, function(x) if(is.matrix(x)) return(x[,inter]) else return(x))
    logging("Running DESeqDataSetFromTximport()", .module="deseq2")
    dds <- DESeq2::DESeqDataSetFromTximport(txi.filtered,
                                            metadata[inter,, drop=FALSE],
                                            design)
 
    logging("Running estimateSizeFactors()", .module="deseq2")
    ## TODO: estimateSizeFactorsForMatrix. Necessary ??!?
    dds <- DESeq2::estimateSizeFactors(dds)
    ## TODO: custom filter because too low expression
    idx <- rowSums(DESeq2::counts(dds, normalized=TRUE) >= 5 ) >= ncol(DESeq2::counts(dds))/2
    dds <- dds[idx, ]

    logging("Running DESeq()", .module="deseq2")
    BiocParallel::register(BiocParallel::MulticoreParam(pipeline$option$nthreads))
    dds <- DESeq2::DESeq(dds, parallel=TRUE)

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
deseq2Results  <- function (pipeline, c1, c2,
                            condition="condition", name=NULL, ...)
{
    logging("Running deseq2Results()", .module="deseq2")
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

    pipeline <- pushResults(pipeline, res, name=name)
    
    return (pipeline)
}

## pipeline() %>% metadata("data/metadata.tsv") %>%
##     options(
##         salmon.index="~/data/ensembl/GRCh38/Homo_sapiens.GRCh38.cdna.all.salmon.index",
##         gtf="~/data/ensembl/GRCh38/Homo_sapiens.GRCh38.95.gtf.gz") %>%
##     salmonGenes()
