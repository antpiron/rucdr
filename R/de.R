
library(tximport)
library(GenomicFeatures)
library(DESeq2)

#' Load a gtf file into a tx handle. A cache is created to improve the loading
#' time (cache name `filename`.sqlite).
#'
#' @param filename The gtf filename
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
                                c("filename", "name", "condition")),
                        txdb)
{
    names(sampleTable) <- sampleTable$name
    
    txi.isoforms <- tximport(sampleTable$filename,
                             type = "salmon", txOut = TRUE,
                             ignoreTxVersion=T)
    k <- keys(txdb, keytype = "TXNAME")
    tx2gene <- biomaRt::select(txdb, k, "GENEID", "TXNAME")

    s <-  structure(list(txdb=txdb,
                         samplesTable=sampleTable,
                         txi=txi.isoforms,
                         isoforms=TRUE,
                         tx2gene=tx2gene,
                         dds=NULL
                         ),
                    class = "samples")

    return(s)
}

#' @param sample A list of samples as returned by newSamples()
#' @param txdb A genomic database (gtf)
#' @export
summarizeSamples <- function(samples)
{
    samples$txi  <- summarizeToGene(samples$txi, samples$tx2gene,
                                    ignoreTxVersion=T)
    samples$isoforms=FALSE
        
    return(samples)
}

#' @export
differential.expression <- function (samples, design=~condition)
{
    samples$dds <- DESeqDataSetFromTximport(samples$txi,
                                    samples$sampleTable,
                                    design)
    ## filter too low expression (TODO: Adapt)
    samples$dds <- dds[ apply(counts(dds), 1,
                              function (x) sum(x > 5) > 4), ]
    samples$dds <- estimateSizeFactors(dds)
    samples$dds <- DESeq(dds, parallel=T)
    
    return(samples)
}




#' @export
condition.samples  <- function (samples, c1, c2)
{
    res <- results(samples$dds, contrast=c("condition", c1, c2))
    if (samples$isoforms)
        {
            rownames(res)  <- sapply(strsplit(rownames(res),'\\.'),'[',1)
            res <- merge(as.data.frame(res), samples$txtogene,
                         by.x=0, by.y="TXNAME", all.x=T)
            rownames(res) <- res$Row.names
            res <- res[-1]
        }
    return (res)
}
