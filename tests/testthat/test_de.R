context("DE")
library(rucdr)
library(magrittr)



set.seed(1)
nsamples <- 10
sampleTable <- data.frame(salmon.quant.sf=file.path("output",
                                                    paste0("S",1:nsamples),
                                                    "quant.sf"),
                          name=paste0("S",1:nsamples),
                          condition=c(rep("ctl", nsamples/2),
                                      rep("pal", nsamples/2)))
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

           dir=file.path("output", entry[["name"]])
           dir.create(dir, recursive=T)
           
           write.table(quant.sf,
                       file=entry[["salmon.quant.sf"]],
                       row.names=FALSE, quote=FALSE, sep="\t")
       })

on.exit(unlink("output", recursive=T),
        add = TRUE)

pl <- pipeline()
pl$metadata <- sampleTable

pl <- pl %>% options(gtf="data/test.gtf") %>%  salmon() 

## print(pl$metadata)
## print(pl$salmon)

test_that("salmon()", {
    expect_s3_class(pl, "pipeline")
    expect_true(!is.null(pl$salmon$txi.isoforms))
    expect_equal(dim(pl$salmon$txi.isoforms$counts),
                 c(length(transcripts), nrow(sampleTable)))
    expect_equal(colnames(pl$salmon$txi.isoforms$counts),
                 as.character(sampleTable$name))
})


pl <- pl %>% deseq2()
test_that("deseq2()", {
    expect_s3_class(pl, "pipeline")
})

res <- pl %>% deseq2Results("ctl", "pal")
test_that("deseq2Results()", {
    expect_s3_class(res, "data.frame")
    expect_true(res["ENST00000456328",]$padj < 0.05)
    expect_true(all(res[transcripts[! ("ENST00000456328" == transcripts)],]$padj > 0.05))
})


## TODO: test summarizeSamples
