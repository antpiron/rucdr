

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
    pipeline$metadata.selection <- as.data.frame(dplyr::filter(pipeline$metadata, ...))
    row.names(pipeline$metadata.selection) <- pipeline$metadata.selection$id
    
    pipeline$metadata.selection <- pipeline$metadata.selection[
                                                ! is.na(pipeline$metadata.selection$salmon.quant.sf),]

    cols <- apply( pipeline$metadata.selection, 1, is.factor)
    pipeline$metadata.selection[, cols]  <-  lapply(pipeline$metadata.selection[, cols], factor)
    ## print("=========================")
    ## print(pipeline$metadata.selection$salmon.quant.sf)
    filenames <- pipeline$metadata.selection[
                              file.exists(
                                  as.character(pipeline$metadata.selection$salmon.quant.sf)),"salmon.quant.sf"]
    if( length(filenames) > 0)
    {

        pipeline$salmon$txi.isoforms <- tximport::tximport(
                                                      filenames,
                                                      type = "salmon", txOut = TRUE,
                                                      ignoreTxVersion=T)
        
        if ( "genes" == pipeline$salmon$tlast.txi)
        {
            pipeline <- salmonGenes(pipeline)
        }
    }


    return(pipeline)
}
