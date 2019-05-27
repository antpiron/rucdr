

#' @export
runSalmon <- function (fq1, fq2, id,
                       outputdir=file.path("output", "salmon", id),
                       nthreads=4,
                       index="")
{
    logging(paste0("runSalmon() for ", id), .module="salmon")

    output.quant.sf <- file.path(outputdir, "quant.sf")
    logging(paste0("runSalmon()  output.quant.sf = ", output.quant.sf),
            .module="salmon")
    if (file.exists(output.quant.sf))
    {
        logging(paste0('runSalmon(): ', output.quant.sf,
                       " already exists."),
                .module="salmon")
        return(output.quant.sf)
    }
    
    paired <- ! is.null(fq2) &&  ! is.na(fq2) && "" != fq2
    if (paired)
        logging(paste0("runSalmon(): ", id, " is paired."),
                .module="salmon")
    
    if (! checkFile(fq1,
                    .message=paste0("runSalmon(): fastq1 ",
                                    fq1, " does not exist")) )
        return(NA)
    
    if (paired && ! checkFile(fq2,
                              .message=paste0("runSalmon(): fastq2 ",
                                              fq2, " does not exist")) )
        return(NA)
           
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

    logging(paste0("runSalmon() salmon running for ", id),
            .module="salmon")
    ret <- system2("salmon", param, wait = TRUE)
    logging(paste0("runSalmon() salmon finished for ", id),
            .module="salmon")

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
        message( paste0(
            "Salmon do not seem to have created the file properly ",
            tmp.quant.sf, ". Inexistent or truncated file.") )
        return(NA)
    }
    logging(paste0("'runSalmon(): renaming ", tmp.dir,
                   " to ", outputdir, "."),
            .module="salmon")
    ##dir.create(outputdir, recursive=T)
    file.move(tmp.dir, outputdir)
   
    return(output.quant.sf)
}


#' @export
salmonQuant <- function (metadata, index="", outputdir="output/salmon",
                         njobs=1, nthreads=4)
{
    checkColumns(metadata, c("fastq1", "id"), .stop=T)
    checkCharacter(
        metadata$fastq1,
        .message="salmonQuant(): metadata$fastq1 is not characters.",
        .stop=T)
    checkCharacter(
        metadata$id,
        .message="salmonQuant(): metadata$id is not characters.",
        .stop=T)

    logging("salmonQuant() parallel quantification", .module="salmon")
    indexes <- 1:nrow(metadata)
    if (! dir.exists(outputdir))
        dir.create(outputdir, recursive=T)
    filenames <- parallel::mclapply(
                               indexes,
                               function (i)
                               {                                  
                                   as.character(runSalmon(
                                       metadata$fastq1[i],
                                       metadata$fastq2[i],
                                       metadata$id[i],
                                       outputdir=file.path(outputdir,
                                                           metadata$id[i]),
                                       nthreads=nthreads,
                                       index=index
                                   ))
                               },
                               mc.cores=njobs,
                               mc.silent=FALSE)
    
    logging("salmonQuant() parallel quantification finished",
            .module="salmon")

    filenames <- as.character(filenames)
    names(filenames) <- metadata$id[indexes]
    
    metadata[,"quant.sf.fn"] <- NA
    metadata[indexes,"quant.sf.fn"] <- filenames

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
rnaseq.salmon <- function (input, nthreads=4)
{
    checkColumns(input, c("quant.sf.fn", "id"), .stop=T)
    ## print(input)
    logging("rnaseq.salmon(): Loading quant.sf in txi object.",
            .module="salmon")
    filenames <- input$quant.sf.fn
    names(filenames) <- input$id
    txi.tx <- tximport::tximport(filenames, type = "salmon",
                                 txOut = TRUE, ignoreTxVersion=T)
     
    structure(list(
        txi = txi.tx
        ), class  = c("salmon_isoforms", "rnaseq_quantification") )
}



## df <- data.frame(fastq1=c("/home/apiron/git/bioinfo/ph.d/tools/omics/fadista/SRR957910_1.fastq.gz"), fastq2=c("/home/apiron/git/bioinfo/ph.d/tools/omics/fadista/SRR957910_2.fastq.gz"), id=c("SRR957910"),stringsAsFactors=F)
## salmonQuant(df, index="~/data/ensembl/GRCh38/Homo_sapiens.GRCh38.cdna.all.salmon.index")
