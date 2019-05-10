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

checkFile <- function(file, ...)
{
    checkAction(! is.na(file) || ! is.null(file) ||  file.exists(file), ...)
}

checkCharacter <- function(str, ...)
{
    checkAction(is.character(str), ...)
}


checkColumns <- function (df, ...)
{
    UseMethod("checkColumns", df)
}


checkColumns.data.frame <- function(df, mandatory, .message="", ...)
{
    checkAction(
        all(mandatory %in% colnames(df)),
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


logging <- function (message, .level=0, .module=NULL)
{
    flag <- is.null(.module) || ! exists("logmodules") || is.null(logmodules) ||
        (exists("logmodules") && .module %in% logmodules)
    
    if ( exists("loglevel") && .level < loglevel && flag)
        message(paste0("Logging: ", message))
}


file.move <- function(src, dst)
{
    ## fail because files can be on different filesystems and R sucks
    ## file.rename(tmp.dir, salmon.output.dir)
    ## TODO: Do something portable
    ret <- system2("mv", c(src, dst), wait = TRUE)
    if (0 != ret)
        stop(paste0("runSalmon(): mv ", src, " ", dst,
                    " returned ", ret))
}
