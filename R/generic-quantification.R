library(tibble)
library(dplyr)


#' Run quantification
#' 
#' @param pipeline A pipeline object
#' @param ... other arguments
#' @export
rnaseq <- function (df, ...)
{
    UseMethod("rnaseq", df)
}

#' Base object for rnaseq quantification
#' 
#' @param input A data frame describing the input files.
#'    `input` should have at least the columns `quant.sf.fn` and `id`.
#' @return An `rnaseq_quantification` S3 object with members
#'    $counts, $tpm, $rpkm. 
#' @export
rnaseq.data.frame <- function (input)
{
    checkColumns(input, c("quant.sf.fn", "id"), .stop=T)

    data <- lapply(input$quant.sf.fn,
                   function (fn)
                   {
                       rt <- read.table(as.character(fn),
                                        sep= "\t", header = TRUE,
                                        quote="", fill=T, check.names=F,
                                        stringsAsFactors=F)
                       rt
                   })
    names(data) <- input$id

    counts <- lmerge(data, on="Name", col="NumReads", input$id)
    tpm <- lmerge(data, on="Name", col="TPM", input$id)
    
    structure(list(
        data   = data,
        counts = counts,
        tpm    = tpm,
        rpkm   = matrix(double())
        ), class  = c("rnaseq_quantification") )
}

print.rnaseq_quantification <- function (rnaseq)
{
    print(paste0("Data: ", paste0(names(rnaseq$data), collapse=" ")))
    print("Counts:")
    print(head(rnaseq$counts))
    print("TPM:")
    print(head(rnaseq$tpm))
    print("RPKM:")
    print(head(rnaseq$rpkm))          
}

## df <- data.frame(quant.sf.fn=c("/home/apiron/git/bioinfo/ph.d/tools/tiger/output/ensembl-GRCh38-95/ULB-T2D261010/salmon/quant.sf", "/home/apiron/git/bioinfo/ph.d/tools/tiger/output/ensembl-GRCh38-95/ULB-T2D310311/salmon/quant.sf"), id=c("ULB-T2D261010","ULB-T2D310311"))


