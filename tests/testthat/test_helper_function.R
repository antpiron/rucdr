context("helper-function")
library(rucdr)




test_that("checkFile() is fine", {
    suppressWarnings({
        expect_false(checkFile(NA))
        expect_false(checkFile(NULL))
        expect_false(checkFile("jkkkk"))
    })
})

test_that("checkCharacter() is fine", {
    suppressWarnings({
        expect_false(checkCharacter(NA))
        expect_false(checkCharacter(NULL))
        expect_false(checkCharacter(factor("jjjj")))
        expect_true(checkCharacter("jkkkk"))
    })
})

test_that("checkColumns() is fine", {
    df <- data.frame(a=1, b=2, c=3)
    expect_true(checkColumns(df, c("a","b")))
    suppressWarnings({
        expect_false(checkColumns(df, c("a","b","d")))
    })
})

test_that("lmerge() is fine", {
    data <- list(
        data.frame(a="1", b="2", c="3", stringsAsFactors=F),
        data.frame(a="1", b="2", c="3", stringsAsFactors=F)
    )
    res <- lmerge(data, on="a", col="b", col.names=c("one", "two"))
    expect_true(all(res == "2"))
    expect_equal(colnames(res), c("one", "two"))
})

test_that("file.move() is fine", {
    suppressWarnings({file.remove(c("data/test", "data/test1"))})
    x <- data.frame()
    write.table(x, file='data/test', col.names=FALSE)
    on.exit(suppressWarnings(
    {file.remove(c("data/test", "data/test1"))}))
    file.move("data/test", "data/test1")
    expect_true(file.exists("data/test1"))
})

test_that("concat.data.frame() is fine", {
    res <- concat.data.frame(data.frame(id=1:2),
                             data.frame(id=3, plop="a"))

    expect_equal(res$id, 1:3)
    expect_equal(colnames(res), c("id", "plop"))
    expect_true(is.na(res$plop[1]))
})
