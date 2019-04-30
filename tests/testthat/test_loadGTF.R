context("Load GTF")
library(rucdr)

txdb <- loadGTF("data/test.gtf")
on.exit(file.remove("data/test.gtf.sqlite"))

test_select <- function ()
{
    k <- keys(txdb, keytype = "TXNAME")
    tx2gene <- select(txdb, k, "GENEID", "TXNAME")
    check <- tx2gene$GENEID[tx2gene$TXNAME == "ENST00000456328"] ==
        "ENSG00000223972"
    return(check)
}


test_that("txdb is a S4 txdb", {
    expect_type(txdb, "S4")
    expect_true(test_select())
})

