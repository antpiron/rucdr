


salmon_process <- function (pipeline, sample)
{
    paired <- ! is.null(sample$fastq2) &&  ! is.na(sample$fastq2) && "" != sample$fastq2
    salmon.output.dir <- file.path(pipeline$option$output.dir, "salmon",
                                   sample$name)
    salmon.output.quant.sf <- file.path(salmon.output.dir, "quant.sf")
    if (file.exists(salmon.output.quant.sf))
    {
        return(salmon.output.quant.sf)
    }
    
    tmp.dir <- tempfile(pattern="salmon-")
    dir.create(tmp.dir)

    fastq.param <- if (paired) {
                       c("-1", as.character(sample$fastq1),
                         "-2", as.character(sample$fastq2))
                   } else {
                       c("-r", as.character(sample$fastq1)) }

    param <- c("quant", "--seqBias", "--gcBias",
               "--no-version-check", "--validateMappings",
               "-p", pipeline$option$nthreads,
               "-i", as.character(pipeline$option$salmon.index),
               "-l", "A",
               fastq.param,
               "-o", tmp.dir)

    ret <- system2("salmon", param, wait = TRUE)

    if (0 != ret)
    {
        warning( paste0("Salmon returned ", ret, ". Check ", tmp.dir,".") )
        return(NA)
    }
    file.rename(tmp.dir, salmon.output.dir)
    
    return(salmon.output.quant.sf)
}

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
    if (is.null(pipeline$metadata$salmon.quant.sf))
    {
        if (is.null(pipeline$metadata$fastq1))
        {
            warning("No fastq1 column in metadata")
            return(pipeline)
        }
        indexes <- which(!is.na(pipeline$metadata$fastq1))
        filenames <- parallel::mclapply(
                                   indexes,
                                   function (i)
                                   {
                                       salmon_process(pipeline,
                                                      as.list(pipeline$metadata[i,]))
                                   },
                                   mc.cores=pipeline$option$njobs)
        
        filenames <- unlist(filenames)
    }
    else
    {
        indexes <- which(!is.na(pipeline$metadata$salmon.quant.sf))
        filenames <- pipeline$metadata$salmon.quant.sf[indexes]
        names(filenames) <- pipeline$metadata$name[indexes]
    }

    names(filenames) <- pipeline$metadata[indexes,"name"]
    
    pipeline$metadata[,"salmon.quant.sf"] <- NA
    pipeline$metadata[indexes,"salmon.quant.sf"] <- filenames
    
    filenames <- filenames[! is.na(filenames) ]
    
    pipeline$salmon = list()
    if (length(filenames) > 0)
        pipeline$salmon$txi.isoforms <- tximport::tximport(
                                                      as.character(filenames),
                                                      type = "salmon", txOut = TRUE,
                                                      ignoreTxVersion=T)
    else
        pipeline$salmon$txi.isoforms <- NULL

    pipeline <- tx2genes(pipeline)

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
        
    return(pipeline)
}
