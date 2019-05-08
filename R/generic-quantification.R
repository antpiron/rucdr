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

    columns_to_matrix <- function (column)
    {
        counts <- Reduce(
            function (df, quant.sf)
            {
                nr <- quant.sf[,c("Name", column)]
                full_join(df, nr, by = "Name")
            },
            data, data.frame(Name=character(), stringsAsFactors=F))   
        
        row.names(counts) <- counts$Name
        counts <- subset(counts, select=-c(Name))
        colnames(counts) <- input$id

        as.matrix(counts)
    }

    counts <- columns_to_matrix("NumReads")
    tpm <- columns_to_matrix("TPM")
    
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


