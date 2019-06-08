context("Pipeline GTF")
library(rucdr)
library(magrittr)


pl <- data.frame(id=character()) %>%
    pipeline() %>% options(gtf="data/test.gtf") %>%  GTF()
on.exit(file.remove("data/test.gtf.sqlite"))

test_select <- function (txdb)
{
    k <- biomaRt::keys(txdb, keytype = "TXNAME")
    tx2gene <- biomaRt::select(txdb, k, "GENEID", "TXNAME")
    check <- tx2gene$GENEID[tx2gene$TXNAME == "ENST00000456328"] ==
        "ENSG00000223972"
    return(check)
}


test_that("pl is a pipeline()", {
    expect_s3_class(pl, "pipeline")
    expect_true(! is.null(pl$option))
    expect_true(! is.null(pl$option$gtf))
    expect_true(! is.null(pl$txdb))
    expect_type(pl$txdb, "S4")
    expect_true(test_select(pl$txdb))
})


## test cache
pl <- data.frame(id=character()) %>%
    pipeline() %>% options(gtf="data/test.gtf") %>%  GTF()
test_that("pl is a pipeline()", {
    expect_s3_class(pl, "pipeline")
    expect_true(! is.null(pl$option))
    expect_true(! is.null(pl$option$gtf))
    expect_true(! is.null(pl$txdb))
    expect_type(pl$txdb, "S4")
    expect_true(test_select(pl$txdb))
})

