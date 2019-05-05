

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
    pipeline$metadata.selection <- dplyr::filter(pipeline$metadata, ...)

    return(pipeline)
}
