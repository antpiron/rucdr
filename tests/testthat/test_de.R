context("DE")
library(rucdr)
library(magrittr)


return()


set.seed(1)
nsamples <- 10
sampleTable <- data.frame(filename=paste0("data/S",1:nsamples,
                                          ".quant.sf"),
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
           write.table(quant.sf,	file=entry[["filename"]],
                       row.names=FALSE, quote=FALSE, sep="\t")
       })

on.exit(sapply(as.character(sampleTable$filename),
               function (f) file.remove(f)),
        add = TRUE)

samples <- newSamples(sampleTable, txdb)

##print(samples$txi)

test_that("newSamples()", {
    expect_s3_class(samples, "samples")
    expect_true(all(samples$sampleTable$name == sampleTable$name))
    expect_equal(dim(samples$txi$counts),
                 c(length(transcripts), nrow(sampleTable)))
})

de <- differential.expression(samples)
test_that("differential.expression()", {
    expect_s3_class(de, "samples")
})

res <- differential.expression.results(de, "ctl", "pal")
test_that("differential.expression.reults()", {
    expect_s3_class(res, "data.frame")
    expect_true(res["ENST00000456328",]$padj < 0.05)
    expect_true(all(res[transcripts[! ("ENST00000456328" == transcripts)],]$padj > 0.05))
})


res.pipe <- sampleTable %>% 
    newSamples(suppressWarnings(loadGTF("data/test.gtf"))) %>%
    differential.expression %>%
    differential.expression.results("ctl", "pal")

print(res.pipe)
test_that("differential.expression.results() with pipe", {
    expect_s3_class(res.pipe, "data.frame")
    expect_true(res.pipe["ENST00000456328",]$padj < 0.05)
    expect_true(all(res.pipe[transcripts[! ("ENST00000456328" == transcripts)],]$padj > 0.05))
})

## TODO: test summarizeSamples
