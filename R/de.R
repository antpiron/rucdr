
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
    if (is.null(pipeline$salmon) || is.null(pipeline$salmon$tlast.txi) )
    {
        warning("Have you called salmon()? Doing nothing! No deseq2 call!")
        return(pipeline)
    }

    txi <- if ( "isoforms" == pipeline$salmon$tlast.txi ) pipeline$salmon$txi.isoforms else pipeline$salmon$txi.genes
    
    if (is.null(txi))
    {
        warning(paste0("No txi for", pipeline$salmon$tlast.txi ,
                       ". Doing nothing! No deseq2 call!"))
        return(pipeline)
    }

    metadata <- pipeline$metadata.selection
    ## TODO: is metadata corresponding to txi.isoforms columns?
    dds <- DESeq2::DESeqDataSetFromTximport(txi,
                                            metadata,
                                            design)
    ## from filter()
    ## print("=========")
    ## print(colnames(dds))
    ## print(setdiff(pipeline$metadata.selection$id, colnames(dds)))
    if (! is.null(pipeline$metadata.selection) )
        dds <- dds[,pipeline$metadata.selection$id]
    ## TODO: custom filter because too low expression
    ## print("=========")
    ## print(colnames(dds))
    dds <- dds[ apply(DESeq2::counts(dds), 1,
                      function (x) sum(x > 5) > ncol(dds)/2), ]
    dds <- DESeq2::estimateSizeFactors(dds)
    ## print("=========")
    ## print(colnames(dds))
    ## print(dim(dds))
    dds <- DESeq2::DESeq(dds, parallel=T)
    
    if (is.null(pipeline$dseq2))
        pipeline$dseq2 <- list()
        
    if ("isoforms" == pipeline$salmon$tlast.txi)
        pipeline$dseq2$dds.isoforms  <- dds
    else
        pipeline$dseq2$dds.genes <- dds
    
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
    isoforms <- if (is.null(isoforms)) "isoforms" == pipeline$salmon$tlast.txi else isoforms
    dds <- if (isoforms) pipeline$dseq2$dds.isoforms else pipeline$dseq2$dds.genes
    res <- DESeq2::results(dds, contrast=c(condition, c1, c2), ...)
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
