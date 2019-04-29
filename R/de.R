
library(tximport)
library(GenomicFeatures)
library(DESeq2)

#' @export
loadGTF <- function (filename)
{
    sqlite <- paste0(filename, ".sqlite")
    if (file.exists(sqlite))
        txdb <- loadDb(sqlite)
    else
    {
        txdb <- makeTxDbFromGFF(filename)
        saveDb(txdb,sqlite)
    }

    return(txdb)
}

#'
#' @param sampleTable A data frame with at least 3 columns
#'     "filename", "name" and "condition"
#' @export
newSamples <- function (sampleTable=
                            setNames(
                                data.frame(matrix(ncol = 3, nrow = 0)),
                                c("filename", "name", "condition")))
{
    names(sample) <- sampleTable$name
    
    txi.isoforms <- tximport(sampleTable$filename,
                             type = "salmon", txOut = TRUE,
                             ignoreTxVersion=T)
    s <-  structure(list(txdb=txdb,
                         samplesTable=sampleTable,
                         txi=txi.isoforms,
                         ),
                    class = "samples")

    return(s)
}

#' @param sample A list of samples as returned by newSamples()
#' @param txdb A genomic database (gtf)
#' @export
summarizeSamples <- function(samples, txdb)
{
    k <- keys(txdb, keytype = "TXNAME")
    tx2gene <- select(txdb, k, "GENEID", "TXNAME")

    samples$txi  <- summarizeToGene(samples$txi, tx2gene,
                                    ignoreTxVersion=T)
    
    return(samples)
}

#' @export
differential.expression <- function (samples, design=~condition)
{
    dds <- DESeqDataSetFromTximport(samples$txi,
                                    samples$sampleTable,
                                    design)
    ## filter too low expression
    dds <- dds[ apply(counts(dds),1,function (x) sum(x > 5) > 4), ]
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds, parallel=T)
    
    return(dds)
}
