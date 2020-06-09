library(dplyr)



checkAction <- function(flag, .message="", .stop=F, .warning=T)
{
    if (! flag)
    {
        if (.stop)
            stop(.message)

        if (.warning)
            warning(.message)
    }
    
    return(flag)
}

checkFile <- function(file, ...)
{
    checkAction( is.character(file) && file.exists(file), ... )
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
    init <- data.frame()
    init[,on] <- character(0)
    init <- as.data.frame(init)
    counts <- Reduce(
        function (df, quant.sf)
        {
            nr <- quant.sf[,c(on, col), drop=FALSE]
            dplyr::full_join(df, nr, by = on)
        },
        data, init)   
    
    row.names(counts) <- counts[,on]
    counts <- counts[,-which(colnames(counts) == on ), drop=FALSE]
    if (! is.null(col.names) )
        colnames(counts) <- col.names
    
    as.matrix(counts)
}

#' @export
logging <- function (message, .level=0, .module=NULL)
{
    flag <- is.null(.module) || ! exists("logmodules") || is.null(logmodules) ||
        (exists("logmodules") && .module %in% logmodules)
    
    if ( exists("loglevel") && .level < loglevel && flag)
        message(paste0("Logging: ", message, "\n"))
}

#' @export
file.move <- function(src, dst)
{
    ## fail because files can be on different filesystems and R sucks
    ## file.rename(tmp.dir, salmon.output.dir)
    ## TODO: Do something portable
    ## dir.create(dirname(dst), recursive=T)
    ret <- system2("mv", c(src, dst), wait = TRUE)
    if (0 != ret)
        stop(paste0("file.move(): mv ", src, " ", dst,
                    " returned ", ret))
}

concat.data.frame <- function(d1, d2) {
  d1.names <- colnames(d1)
  d2.names <- colnames(d2)

  d2.add <- setdiff(d1.names, d2.names)
  d1.add <- setdiff(d2.names, d1.names)

  if(length(d2.add) > 0) {
      d2[, d2.add] <- NA
  }

  if(length(d1.add) > 0) {
      d1[,d1.add] <- NA
    }

  return(rbind(d1, d2))
}

# TODO: Test
is.not.empty <- function (x)
{
    !is.na(x) && !is.null(x) && x != ""
}

#' @export
import <- function (filename, paths = c())
{
    if ( startsWith(filename, '/') )
        source(filename)
    else
    {
        wd <- tryCatch({dirname(sys.frame(1)$ofile)},
                       error=function (e) {file.path(".")})
        path <- file.path(wd, filename)
        
        if (! file.exists(path) )
        {
            paths <- c(paths, strsplit(Sys.getenv("R_IMPORT_DIR"), ':')[[1]])
            for ( cpath in paths )
            {
                path <- file.path(cpath, filename)
                if (file.exists(path) )
                {
                    source(path)
                    break
                }
            }
        }
    }
}

