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
    ## TODO: test this
    getFilter(pipeline)$quant.sf.fn <- quant$quant.sf.fn
    quant <- quant[! is.na(quant$quant.sf.fn),]
   
    res <- rnaseq(quant, nthreads=pipeline$option$nthreads)

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

    logging("salmonGenes.pipeline(): generating counts/tpm/length/effective_length.",
            .module="salmon")
    ToDoList <- list(isoforms$counts, isoforms$tpm,
                     isoforms$length, isoforms$effective_length)
    listCounts <- parallel::mclapply(1:4,
                                function(i)
                                    if (i <= 2)
                                        sumCounts(ToDoList[[i]])
                                    else
                                        meanLength(ToDoList[[i]]),
                                mc.cores=pipeline$option$nthreads)
    logging("salmonGenes.pipeline(): generating counts/tpm/length/effective_length finished.",
            .module="salmon")
   
    res <- structure(list(
        data   = isoforms$data,
        counts = listCounts[[1]],
        tpm    = listCounts[[2]],
        rpkm   = matrix(double()),
        length = listCounts[[3]],
        effective_length = listCounts[[4]]
    ), class  = c("salmon_genes", "rnaseq_quantification") )

    
    logging("salmonGenes.pipeline(): pushing results.",
            .module="salmon")
    pipeline <- pushResults(pipeline, res)
    
    return(pipeline)
}
