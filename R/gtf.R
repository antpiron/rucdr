
#' Load a gtf file
#' 
#' @param pipeline A pipeline object
#' @param ... other arguments
#' @export
GTF <- function (pipeline, ...)
{
    UseMethod("GTF", pipeline)
}

#' Load a gtf file. A cache is created to improve the loading
#' time (cache name `filename`.sqlite).
#'
#' @param pipeline A pipeline object
#' @param filename The gtf filename
#' @rdname GTF
#' @export
GTF.pipeline <- function (pipeline)
{
    filename <- pipeline$option$gtf
    print(filename)
    sqlite <- paste0(filename, ".sqlite")
    if (file.exists(sqlite))
        pipeline$txdb <- AnnotationDbi::loadDb(sqlite)
    else
    {
        pipeline$txdb <- GenomicFeatures::makeTxDbFromGFF(filename)
        AnnotationDbi::saveDb(pipeline$txdb, sqlite)
    }

    return(pipeline)
}

tx2genes <- function (pipeline, ...)
{
    UseMethod("tx2genes", pipeline)
}

tx2genes.pipeline <- function (pipeline)
{
    if (is.null(pipeline$txdb))
        pipeline <- GTF(pipeline)

    if (! is.null(pipeline$tx2gene) )
        return(pipeline)   
    
    k <- biomaRt::keys(pipeline$txdb, keytype = "TXNAME")
    pipeline$tx2gene <- biomaRt::select(pipeline$txdb, k,
                                        "GENEID", "TXNAME")
    row.names(pipeline$tx2gene) <- pipeline$tx2gene$TXNAME

    return(pipeline)
}

