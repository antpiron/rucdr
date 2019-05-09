library(dplyr)


checkAction <- function(flag, .message="", .stop=F, .warning=T)
{
    if (! flag)
    {
        if (.stop)
            stop(.message)

        if (.warning)
            warning(message)
    }
    
    return(flag)
}

checkColumns <- function (df, ...)
{
    UseMethod("checkColumns", df)
}

checkColumns.data.frame <- function(df, mandatory, .message="", ...)
{
    checkAction(
        all(colnames(df) %in% mandatory),
        paste0(
            .message,
            "Expected column in data frame missing, mandatory: ",
            paste0(mandatory, collapse = ", "), "."),
        ...)
}

lmerge <- function (data, on, col, col.names)
{
    counts <- Reduce(
        function (df, quant.sf)
        {
            nr <- quant.sf[,c(on, col)]
            full_join(df, nr, by = on)
        },
        data, data.frame(Name=character(), stringsAsFactors=F))   
    
    row.names(counts) <- counts[,on]
    counts <- counts[,-which(colnames(counts) == on )]
    if (! is.null(col.names) )
        colnames(counts) <- col.names
    
    as.matrix(counts)
}
