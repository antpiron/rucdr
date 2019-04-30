context("DE")
library(rucdr)

txdb <- suppressWarnings(loadGTF("data/test.gtf"))
on.exit(file.remove("data/test.gtf.sqlite"))

test_select <- function ()
{
    k <- biomaRt::keys(txdb, keytype = "TXNAME")
    tx2gene <- biomaRt::select(txdb, k, "GENEID", "TXNAME")
    check <- tx2gene$GENEID[tx2gene$TXNAME == "ENST00000456328"] ==
        "ENSG00000223972"
    return(check)
}


test_that("txdb is a S4 txdb", {
    expect_type(txdb, "S4")
    expect_true(test_select())
})

sampleTable <- data.frame(filename=paste0("data/S",1:10, ".quant.sf"),
                          name=paste0("S",1:10),
                          condition=c(rep("ctl", 5), rep("pal", 5)))
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
           quant.sf$TPM <- c(rnorm(length(transcripts), 100, 50))
           if ("pal" == entry[["condition"]])
           { quant.sf$TPM[1] <- rnorm(1, 1000, 200) }
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

res <- condition.samples(de, "ctl", "pal")
test_that("differential.expression()", {
    expect_s3_class(res, "data.frame")
    expect_true(res["ENST00000456328",]$padj < 0.01)
    expect_true(all(res[transcripts[! ("ENST00000456328" == transcripts)],]$padj > 0.05))
})

