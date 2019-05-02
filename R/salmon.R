


salmon_process <- function (pipeline, sample)
{
    paired <- ! is.null(sample$fastq2) &&  ! is.na(sample$fastq2) && "" != sample$fastq2
    if (! file.exists(sample$fastq1) )
    {
        warning( paste0("salmon_process(): file ",
                        sample$fastq1, " does not exist") )
        return(NA)

    }
    if (paired && ! file.exists(sample$fastq2) )
    {
        warning( paste0("salmon_process(): file ",
                        sample$fastq2, " does not exist") )
        return(NA)

    }
    salmon.output.dir <- file.path(pipeline$option$output.dir, "salmon",
                                   sample$id)
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
    tmp.quant.sf = file.path(tmp.dir, "quant.sf")
    if (! file.exists(tmp.quant.sf) ||
        1000 != length(readLines(file(tmp.quant.sf), n = 1000)))
    {
        warning( paste0("Salmon do not seem to have created the file ",
                        tmp.quant.sf, ".") )
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
    if (is.null(pipeline$metadata$salmon.quant.sf) &&
        is.null(pipeline$metadata$fastq1))
    {
        warning("No fastq1 nor salmon.quant.sf column in metadata")
        return(pipeline)
    }

    indexes <- 1:nrow(pipeline$metadata)
    filenames <- parallel::mclapply(
                               indexes,
                               function (i)
                               {
                                    if ( is.null(pipeline$metadata$fastq1) ||
                                        is.null(pipeline$metadata$fastq1[i]))
                                   {
                                       return(as.character(
                                           pipeline$metadata$salmon.quant.sf[i]))
                                   }
                                   
                                   as.character(salmon_process(
                                       pipeline,
                                       as.list(pipeline$metadata[i,])))
                               },
                               mc.cores=pipeline$option$njobs)
    

    filenames <- as.character(filenames)
    names(filenames) <- pipeline$metadata$id[indexes]
    
    pipeline$metadata[,"salmon.quant.sf"] <- NA
    pipeline$metadata[indexes,"salmon.quant.sf"] <- filenames
    
    filenames <- filenames[! is.na(filenames) ]
    
    pipeline$salmon = list()
    ##print(filenames)
    if (length(filenames) > 0)
    {
        pipeline$salmon$txi.isoforms <- tximport::tximport(
                                                      filenames,
                                                      type = "salmon", txOut = TRUE,
                                                      ignoreTxVersion=T)
        pipeline$salmon$tlast.txi <- "isoforms"
    }
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
    pipeline$salmon$tlast.txi <- "genes"

    return(pipeline)
}
