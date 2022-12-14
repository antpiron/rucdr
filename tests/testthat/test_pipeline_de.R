context("DE")
library(rucdr)
library(magrittr)


loglevel <- 1
set.seed(1)
nsamples <- 10
dir.create(file.path("output", "salmon"), recursive=T)
sampleTable <- data.frame(quant.sf.fn=file.path("output", "salmon",
                                                paste0("S",1:nsamples),
                                                "quant.sf"),
                          fastq1=file.path("output", "salmon",
                                           paste0("S",1:nsamples,
                                                  "_R1.fastq.gz")),
                          id=paste0("S",1:nsamples),
                          condition=c(rep("ctl", nsamples/2),
                                      rep("pal", nsamples/2)),
                          stringsAsFactors=F)
row.names(sampleTable)  <- sampleTable$id
transcripts <- c("ENST00000456328", "ENST00000450305",
                 "ENST00000473358", "ENST00000469289",
                 "ENST00000607096", "ENST00000606857")
apply(sampleTable, 1,
       function (entry)
       {
           quant.sf <- data.frame(
               Name=transcripts,
               Length=rep(100,length(transcripts)),
               EffectiveLength=rep(100,length(transcripts)))
           quant.sf$TPM <- rnorm(length(transcripts), 100, 10)
           if ("pal" == entry[["condition"]])
           { quant.sf$TPM[1] <- rnorm(1, 1000, 10) }
           quant.sf$TPM[quant.sf$TPM < 0] <- 0
           quant.sf$TPM  <- quant.sf$TPM * 1E6 / sum(quant.sf$TPM)
           quant.sf$NumReads  <- quant.sf$TPM * 100

           dir=file.path("output", "salmon", entry[["id"]])
           dir.create(dir, recursive=T)
           
           write.table(quant.sf,
                       file=entry[["quant.sf.fn"]],
                       row.names=FALSE, quote=FALSE, sep="\t")
       })

on.exit(unlink("output", recursive=T),
        add = TRUE)

pl <- sampleTable %>% pipeline()
## pl$metadata <- sampleTable

pl <- pl %>% options(gtf="data/test.gtf") %>%  salmon() 

##print(pl$metadata)
##print("===== metadata.selection")
##print(pl$metadata.selection)
##print(pl$salmon)

isoforms <- getResultsByClass(pl, .class = "salmon_isoforms")

test_that("salmon()", {
    expect_s3_class(pl, "pipeline")
    expect_true(!is.null(isoforms))
    expect_equal(dim(isoforms$txi$counts),
                 c(length(transcripts), nrow(sampleTable)))
    expect_equal(colnames(isoforms$txi$counts),
                 as.character(sampleTable$id))
})

## print("====== before deseq2")
## print(isoforms$counts)

pl <- pl %>% deseq2()
test_that("deseq2()", {
    expect_s3_class(pl, "pipeline")
})

pl <- pl %>% deseq2Results("ctl", "pal", name="test")
res <- pl %>% getResultsByName("test")
test_that("deseq2Results()", {
    expect_s3_class(res, "data.frame")
    expect_true(res["ENST00000456328",]$padj < 0.05)
    expect_true(all(res[transcripts[! ("ENST00000456328" == transcripts)],]$padj > 0.05))
})

## test filter
## TODO: better testing
pl <- pl %>% filter(! id %in% c("S1", "S6")) %>% salmon()

pl <- pl %>% deseq2()
test_that("deseq2()", {
    expect_s3_class(pl, "pipeline")
})

pl <- pl %>% deseq2Results("ctl", "pal", name="test2")
res <- pl %>% getResultsByName("test2")
test_that("deseq2Results()", {
    expect_s3_class(res, "data.frame")
    expect_true(res["ENST00000456328",]$padj < 0.05)
    expect_true(all(res[transcripts[! ("ENST00000456328" == transcripts)],]$padj > 0.05))
})

## Summarize by genes
pl <- sampleTable %>% pipeline() %>%
    options(gtf="data/test.gtf") %>%
    salmonGenes()
genes <- getResultsByClass(pl, .class = "salmon_genes")
test_that("salmonGenes()", {
    expect_s3_class(pl, "pipeline")
    expect_true(!is.null(genes))
    expect_s3_class(genes, "salmon_genes")
    ## expect_equal(dim(isoforms$counts),
    ##              c(length(transcripts), nrow(sampleTable)))
    expect_equal(colnames(genes$txi$counts),
                 as.character(sampleTable$id))
})



