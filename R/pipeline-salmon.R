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
    quant <- quant[! is.na(quant$quant.sf.fn),]
    ## TODO: set quant.sf.fn in metadata
    
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

    pipeline <- tx2genes(pipeline)

    rn <- gsub('\\.[0-9]*$', '',
               row.names(isoforms$counts))
    toGenesID <- pipeline$tx2gene[rn,"GENEID"]
    ## print(head(toGenesID))

    sumCounts <- function (counts)
    {
        counts <- aggregate(counts, data.frame(toGenesID), sum)
        row.names(counts) <- counts$toGenesID
        counts[,-1]
    }

    logging("salmonGenes.pipeline(): Computing counts.",
            .module="salmon")
    counts <- sumCounts(isoforms$counts)
    logging("salmonGenes.pipeline(): Computing tpm.",
            .module="salmon")
    tpm <- sumCounts(isoforms$tpm)
    ## print(head(genes))

    meanLength <- function(.length)
    {
        nr <- nrow(isoforms$counts)
        nc <- ncol(isoforms$counts)
        mat <- matrix(rep(1:nr, nc), nrow=nr)
        colnames(mat) <- colnames(isoforms$counts)
        Length <- aggregate(mat,
                            data.frame(toGenesID),
                            function(x)
                            {
                                weighted.mean(.length[x,],
                                              isoforms$counts[x,])
                            })
        row.names(Length) <- Length$toGenesID
        Length[,-1] 
    }

    logging("salmonGenes.pipeline(): Computing length.",
            .module="salmon")
    Length <- meanLength(isoforms$length)
    logging("salmonGenes.pipeline(): Computing effective_length.",
            .module="salmon")
    effective_length <- meanLength(isoforms$effective_length)
    
    res <- structure(list(
        data   = isoforms$data,
        counts = counts,
        tpm    = tpm,
        rpkm   = matrix(double()),
        length = Length,
        effective_length = effective_length
    ), class  = c("salmon_genes", "rnaseq_quantification") )

    
    pipeline <- pushResults(pipeline, res)
    
    return(pipeline)
}
