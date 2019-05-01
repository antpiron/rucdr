

#' Set options
#' 
#' @param pipeline A pipeline object
#' @param ... other arguments
#' @export
options <- function (pipeline, ...)
{
    UseMethod("options", pipeline)
}

#' Set options
#' 
#' @param pipeline A pipeline object
#' @param ... Options
#' @rdname options
#' @export
options.pipeline <- function (pipeline, ...)
{
    pipeline$option <- modifyList(pipeline$option, list(...))
    
    return(pipeline)
}
