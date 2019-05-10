

runSalmon <- function (fq1, fq2, id,
                       outputdir=file.path("output", "salmon", id),
                       nthreads=4,
                       index="")
{
    logging(paste0("runSalmon() for ", id), .module="salmon")

    paired <- ! is.null(fq2) &&  ! is.na(fq2) && "" != fq2
    
    if (! checkFile(fq1,
                    .message=paste0("runSalmon(): fastq1 ",
                                    fq1, " does not exist")) )
        return(NA)
    
    if (paired && ! checkFile(fq1,
                              .message=paste0("runSalmon(): fastq2 ",
                                              fq2, " does not exist")) )
        return(NA)
        
    logging(paste0("runSalmon(): ", id, " is paired."),
            .module="salmon")
    output.quant.sf <- file.path(outputdir, "quant.sf")

    if (file.exists(output.quant.sf))
    {
        logging(paste0('runSalmon(): ', output.quant.sf,
                       " already exists."),
                .module="salmon")
        return(output.quant.sf)
    }
    
    tmp.dir <- tempfile(pattern="salmon-")
    dir.create(tmp.dir)

    fastq.param <- if (paired) {
                       c("-1", as.character(fq1),
                         "-2", as.character(fq2))
                   } else {
                       c("-r", as.character(fq1)) }

    param <- c("quant", "--seqBias", "--gcBias",
               "--no-version-check", "--validateMappings",
               "-p", nthreads,
               "-i", as.character(index),
               "-l", "A",
               fastq.param,
               "-o", tmp.dir)

    ret <- system2("salmon", param, wait = TRUE)

    if (0 != ret)
    {
        message( paste0("'runSalmon(): Salmon returned ",
                        ret, ". Check ", tmp.dir,".") )
        return(NA)
    }
    tmp.quant.sf = file.path(tmp.dir, "quant.sf")
    if (! file.exists(tmp.quant.sf) ||
        1000 != length(readLines(file(tmp.quant.sf), n = 1000)))
    {
        message( paste0("Salmon do not seem to have created the file ",
                        tmp.quant.sf, ".") )
        return(NA)
    }
    logging(paste0("'runSalmon(): renaming ", tmp.dir,
                   " to ", outputdir, "."),
            .module="salmon")
    file.move(tmp.dir, outputdir)
   
    return(output.quant.sf)
}


salmonQuant <- function (metadata, index="", outputdir="output/salmon",
                         njobs=1, nthreads=4)
{
    checkColumns(metadata, c("fastq1", "id"), .stop=T)
    checkCharacter(metadata$fastq1,
                   .message="salmonQuant(): metadata$fastq1 is not characters.",
                   .stop=T)
    checkCharacter(metadata$id,
                   .message="salmonQuant(): metadata$id is not characters.",
                   .stop=T)

    indexes <- 1:nrow(metadata)
    filenames <- parallel::mclapply(
                               indexes,
                               function (i)
                               {                                  
                                   as.character(runSalmon(
                                       metadata$fastq1, metadata$fastq2,
                                       metadata$id,
                                       nthreads=nthreads,
                                       index=index
                                   ))
                               },
                               mc.cores=njobs,
                               mc.silent=FALSE)
    

    filenames <- as.character(filenames)
    names(filenames) <- metadata$id[indexes]
    
    metadata[,"salmon.quant.sf"] <- NA
    metadata[indexes,"salmon.quant.sf"] <- filenames

    structure(metadata,
              class  = c("salmon", "data.frame") )
}



#' Base object for rnaseq quantification
#' 
#' @param input A data frame describing the input files.
#'    `input` should have at least the columns `quant.sf.fn` and `id`.
#' @return An `rnaseq_quantification` S3 object with members
#'    $counts, $tpm, $rpkm. 
#' @export
rnaseq.salmon <- function (input)
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
    ## TODO: RPKM https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/#highlighter_305629
    
    structure(list(
        data   = data,
        counts = counts,
        tpm    = tpm,
        rpkm   = matrix(double())
        ), class  = c("rnaseq_quantification") )
}



## df <- data.frame(fastq1=c("/home/apiron/git/bioinfo/ph.d/tools/omics/fadista/SRR957910_1.fastq.gz"), fastq2=c("/home/apiron/git/bioinfo/ph.d/tools/omics/fadista/SRR957910_2.fastq.gz"), id=c("SRR957910"),stringsAsFactors=F)
## salmonQuant(df, index="~/data/ensembl/GRCh38/Homo_sapiens.GRCh38.cdna.all.salmon.index")
