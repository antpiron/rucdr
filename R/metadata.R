
#' Load the metadata files
#' 
#' @param pipeline A pipeline object
#' @param ... other arguments
#' @export
metadata <- function (files, ...)
{
    UseMethod("metadata", files)
}


#' Load the metadata files
#' 
#' @param pipeline A pipeline object
#' @param filename A tab separated file 
#' @rdname metadata
#' @export
metadata.character <- function (files)
{
    metadata <- Reduce(function (a, b)
    {
        df <- read.table(b, header = TRUE,
                         sep="\t", quote="")
        concat.data.frame(a, df)
    },
    files, data.frame(id="-1", stringsAsFactors=F))
    metadata <- metadata[-1,]
        
    if (is.null(metadata$id)) { stop("Mandatory `id` column missing.") }

    row.names(metadata) <- metadata$id

    structure(metadata,
              class = c("metadata", "data.frame"))
}

metadata.factor <- function(files)
{
    metadata(as.character(files))
}
