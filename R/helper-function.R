

checkColumns <- function (df, ...)
{
    UseMethod("checkColumns", df)
}

checkColumns.data.frame <- function(df, mandatory, .stop=F)
{
    flag <- all(colnames(df) %in% mandatory)
    if (.stop && ! flag)
        stop(paste0("Expected column in data frame missing, mandatory: ",
                    paste0(mandatory, collapse = ", "), "."))
    
    return(flag)
}
