



salmon.process <- function (sample, nthreads=4,
                            index="", output.dir="")
{
    fastq1 <- sample["fastq1"]
    fastq2 <- sample["fastq2"]
    name <- sample["name"]
    paired <- ! is.na(fastq2)

    tmp.dir <- tempfile(pattern="dir-")
    dir.create(tmp.dir)
    
    param <- c("quant", "--seqBias", "--gcBias",
               "--no-version-check", "--validateMappings",
               "-p", nthreads, "-i", index, "-l", "A",
               ifelse(paired,
                      c("-1", fastq1, "-2", fastq2),
                      c("-r", fastq1)),
               "-o", tmp.dir)


    ret <- system2("salmon", param, wait = TRUE)

    if (0 != ret)
    {
        return(NA)
    }
    mydir <- file.path(output.dir,name)
    file.rename(tmp.dir, mydir)
    
    return(file.path(mydir,"quant.sf"))
}


#' @export
salmon <- function (sampleTable, output.dir="output/salmon/",
                    njobs=1, nthreads=4)
{
    parallel::mclapply(1:nrow(sampleTable),
                       function (i)
                       {
                           salmon.process(sampleTable[i,])
                       },
                       mc.cores=njobs)

    return(sampleTable)
}
