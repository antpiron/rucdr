library(magrittr)



#' Run salmon over the fastq files
#' 
#' @param pipeline A pipeline object
#' @param ... other arguments
#' @export
salmon <- function (pipeline, ...)
{
    UseMethod("salmon", pipeline)
}

#' Run salmon over the fastq files
#' 
#' @param pipeline A pipeline object
#' @rdname salmon
#' @export
salmon.pipeline <- function (pipeline)
{
    if (is.null(pipeline$metadata))
    {
        warning("No metadata imported! Nothing to do!")
        return(pipeline)
    }

    pipeline <- pipeline %>%
        filter( (exists("quant.sf.fn") && is.not.empty(quant.sf.fn) &&
                 file.exists(as.character(quant.sf.fn))) ||
                is.not.empty(fastq1))
    metadata_filtered <-  getFilter(pipeline)

    quant <- salmonQuant(metadata_filtered,
                         index=pipeline$option$salmon.index,
                         outputdir=file.path(pipeline$option$output.dir,
                                             "salmon"),
                         njobs=pipeline$option$njobs,
                         nthreads=pipeline$option$nthreads)
    ## TODO: fix filtering. Set metadata.
    getFilter(pipeline)$quant.sf.fn <- quant$quant.sf.fn
    if (is.null(pipeline$metadata$quant.sf.fn))
        pipeline$metadata$quant.sf.fn <- NA
    pipeline$metadata[quant$id, "quant.sf.fn"] <- quant$quant.sf.fn
    quant <- quant[! is.na(quant$quant.sf.fn),]
   
    res <- rnaseq(quant)

    pipeline <- pushResults(pipeline, res)
    
    return(pipeline)
}

#' Run salmon for genes expression
#' 
#' @param pipeline A pipeline object
#' @param ... other arguments
#' @export
salmonGenes <- function (pipeline, ...)
{
    UseMethod("salmonGenes", pipeline)
}

#' Run salmon for genes expression
#' 
#' @param pipeline A pipeline object
#' @rdname salmonGenes
#' @export
salmonGenes.pipeline <- function(pipeline)
{
    logging("salmonGenes.pipeline(): start.",
            .module="salmon")
    isoforms <- getResultsByClass(pipeline, .class = "salmon_isoforms")
    if (is.null(isoforms))
    {
        logging("salmonGenes.pipeline(): No isoforms. Running salmon().",
                .module="salmon")
        pipeline <- salmon(pipeline)
        isoforms <- getResultsByClass(pipeline,
                                      .class = "salmon_isoforms")
    }
    if (is.null(isoforms))
        stop("salmonGenes.pipeline(): isoforms not found.")

    txi <- tximport::summarizeToGene(
                         isoforms$txi,
                         pipeline$tx2gene,
                         ignoreTxVersion=T)
    
    res <- structure(list(
        txi   = txi
    ), class  = c("salmon_genes", "rnaseq_quantification") )

    
    logging("salmonGenes.pipeline(): pushing results.",
            .module="salmon")
    pipeline <- pushResults(pipeline, res)
    
    return(pipeline)
}
