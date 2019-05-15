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
#' @rdname filter
#' @export
filter.pipeline <- function (pipeline, ...)
{
    metadata.filtered <- getFilter(pipeline)
    pipeline$metadata.filtered <- subset(metadata.filtered, ...)
    

    return(pipeline)
}
