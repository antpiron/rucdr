context("Pipeline metadata")
library(rucdr)
library(magrittr)

df <- data.frame(id=c("S1","S2"),
                 fastq1=c("S1_R1.tar.gz", "S2_R2.tar.gz"),
                 fastq2=c("S2_R1.tar.gz", ""))
write.table(df, file="data/metadata.tsv", sep="\t", quote=F)
on.exit(file.remove("data/metadata.tsv"))


pl <- pipeline() %>% metadata("data/metadata.tsv")

test_that("pl is a pipeline()", {
    expect_s3_class(pl, "pipeline")
    expect_true(! is.null(pl$metadata))
    expect_equal(as.character(pl$metadata[1,"id"]), "S1")
})

