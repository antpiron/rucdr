
pipeline <- function (metadata, ...)
{
    UseMethod("pipeline", metadata)
}


#' Initialize the pipeline
#' 
#' @param sampleTable A data frame with at least 3 columns
#'     "fastq1", "name" and "condition"
#' @export
pipeline.metadata <- function (metadata)
{
    ncores <- parallel::detectCores(logical = TRUE)
    njobs <- 1
    nthreads  <- ncores
    
    if (ncores > 8)
    {
        l  <- c(8,6,4,2,1)
        nthreads <- l[which(sapply(l, function (i) 0 == (ncores %% i)))][1]
        njobs=ncores %/% nthreads
    }
    
    structure(list(option=list(
                       output.dir="output/",
                       njobs=njobs,
                       nthreads=nthreads),
                   metadata=metadata),
              class = "pipeline")
}

pipeline.data.frame <- function (metadata)
{
    if (is.null(metadata$id)) { stop("Mandatory `id` column missing.") }

    
    pipeline(structure(metadata,
                       class = c("metadata", "data.frame")))
}
