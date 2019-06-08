context("Pipeline metadata")
library(rucdr)
library(magrittr)


test_that("pl is a pipeline()", {
    df <- data.frame(id=c("S1","S2"),
                     fastq1=c("S1_R1.tar.gz", "S2_R2.tar.gz"),
                     fastq2=c("S2_R1.tar.gz", ""))
    write.table(df, file="data/metadata.tsv", sep="\t", quote=F)
    on.exit(file.remove("data/metadata.tsv"))
    
    metadata <- metadata("data/metadata.tsv")
    expect_s3_class(metadata, "metadata")
    expect_true(! is.null(metadata))
    expect_equal(as.character(metadata[1,"id"]), "S1")
})

