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

    metadata_filtered <- pipeline %>%
        filter( (is.not.empty(quant.sf.fn) &&
                 file.exists(as.character(quant.sf.fn))) ||
                is.not.empty(fastq1)) %>%
        getFilter(pipeline)

    rnaseq(salmonQuant(metadata_filtered))

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
    if (is.null(pipeline$salmon))
        pipeline <- salmon(pipeline)

    if (is.null(pipeline$salmon$txi.isoforms))
    {
        warning("No pipeline$salmon$txi.isoforms")
        return(pipeline)
    }
    
    pipeline$salmon$txi.genes  <- tximport::summarizeToGene(
                                                pipeline$salmon$txi.isoforms,
                                                pipeline$tx2gene,
                                                ignoreTxVersion=T)
    pipeline$salmon$tlast.txi <- "genes"

    return(pipeline)
}
