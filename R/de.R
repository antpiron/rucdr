
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
    if (is.null(pipeline$salmon))
    {
        warning("Have you called salmon()? Doing nothing! No deseq2 call!")
        return(pipeline)
    }

    if (is.null(pipeline$salmon$txi.isoforms))
    {
        warning("No pipeline$salmon$txi.isoforms. Doing nothing! No deseq2 call!")
        return(pipeline)
    }

    metadata <- pipeline$metadata[
                             ! is.na(pipeline$metadata$salmon.quant.sf),]
    ## TODO: is metadata corresponding to txi.isoforms columns?
    dds <- DESeq2::DESeqDataSetFromTximport(pipeline$salmon$txi.isoforms,
                                            metadata,
                                            design)
    ## TODO: custom filter because too low expression 
    dds <- dds[ apply(DESeq2::counts(dds), 1,
                      function (x) sum(x > 5) > 4), ]
    dds <- DESeq2::estimateSizeFactors(dds)
    dds <- DESeq2::DESeq(dds, parallel=T)
    pipeline$dseq2 <- list(dds.isoforms=dds)
    
    return(pipeline)
}


#' Return the results
#' 
#' @param pipeline A pipeline object
#' @param c1 Condition 1
#' @param c2 Condition 2
#' @param condition The column name of sampleTable passed to newSamples()
#'     to use as contrast
#' @param isoforms True if isoforms result (default: True)
#' @export
deseq2Results  <- function (pipeline, c1, c2, condition="condition", isoforms=T)
{
    dds <- if (isoforms) pipeline$dseq2$dds.isoforms else pipeline$dseq2$dds.genes
    res <- DESeq2::results(dds, contrast=c(condition, c1, c2))
    if (isoforms)
    {
        rownames(res)  <- sapply(strsplit(rownames(res),'\\.'),'[',1)
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
