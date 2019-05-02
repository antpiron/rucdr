
#' Load the metadata files
#' 
#' @param pipeline A pipeline object
#' @param ... other arguments
#' @export
metadata <- function (pipeline, ...)
{
    UseMethod("metadata", pipeline)
}


#' Load the metadata files
#' 
#' @param pipeline A pipeline object
#' @param filename A tab separated file 
#' @rdname fastq
#' @export
metadata.pipeline <- function (pipeline, filename="metadata.tsv")
{
    metadata <- read.table(filename, header = TRUE,
                              sep="\t", quote="")
        
    if (is.null(metadata$id)) { stop("Mandatory `id` column missing.") }

    row.names(metadata) <- metadata$id
    pipeline$metadata <- metadata    

    return(pipeline)
}
