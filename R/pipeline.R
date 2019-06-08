
#' @export
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
                   metadata=metadata,
                   results=list()),
              class = "pipeline")
}

#' @export
pipeline.data.frame <- function (metadata)
{
    if (is.null(metadata$id)) { stop("Mandatory `id` column missing.") }

    
    pipeline(structure(metadata,
                       class = c("metadata", "data.frame")))
}


#' @export
getResultsByClass <- function (pipeline, ...)
{
    UseMethod("getResultsByClass", pipeline)
}

#' @export
getResultsByClass.pipeline <- function (pipeline,
                                        .class="rnaseq_quantification")
{
    pos <- Position(function (x) is(x, .class),
                    pipeline$results)

    if (is.na(pos))
        NULL
    else
       pipeline$results[[pos]] 
    
}

#' @export
pushResults <- function (pipeline, ...)
{
    UseMethod("pushResults", pipeline)
}

#' @export
pushResults.pipeline <- function (pipeline, res, name=NULL)
{
    pipeline$results <- append(list(res), pipeline$results)
    if ( ! is.null(name) )
        names(pipeline$results)[1] <- name

    return(pipeline)
}


#' @export
getResultsByName <- function (pipeline, ...)
{
    UseMethod("getResultsByName", pipeline)
}

#' @export
getResultsByName.pipeline <- function (pipeline,
                                       name="1")
{
    pos <- Position(function (x) x == name,
                    names(pipeline$results))

    if (is.na(pos))
        NULL
    else
       pipeline$results[[pos]]   
}


