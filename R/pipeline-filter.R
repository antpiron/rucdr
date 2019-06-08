#' Return the filtered metadata
#' 
#' @param pipeline A pipeline object
#' @param ... other arguments
#' @export
getFilter <- function (pipeline, ...)
{
    UseMethod("getFilter", pipeline)
}

#' @export
getFilter.pipeline <- function (pipeline, ...)
{
    if ( is.null(pipeline$metadata.filtered) )
        pipeline$metadata
    else
        pipeline$metadata.filtered
}

#' @export
"getFilter<-" <- function(pipeline, value)
{
    
    if ( is.null(pipeline$metadata.filtered) )
        pipeline$metadata <- value
    else
        pipeline$metadata.filtered <- value
    pipeline
}

#' @export
unFilter <- function (pipeline, ...)
{
    UseMethod("unFilter", pipeline)
}

#' @export
unFilter.pipeline <- function (pipeline)
{
    pipeline$metadata.filtered <- NULL
}

#' Filter the metadata rows
#' 
#' @param pipeline A pipeline object
#' @param ... other arguments
#' @export
filter <- function (pipeline, ...)
{
    UseMethod("filter", pipeline)
}

#' Filter the metadata using dplyr filter
#' 
#' @param pipeline A pipeline object
#' @param unfiltered 
#' @rdname filter
#' @export
filter.pipeline <- function (pipeline, ..., unfiltered=FALSE)
{
    metadata.filtered <- if (unfiltered) pipeline$metadata else getFilter(pipeline)
    pipeline$metadata.filtered <- subset(metadata.filtered, ...)
    

    return(pipeline)
}
