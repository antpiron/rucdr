
distance.spearman <- function (a, b)
{
    1-cor(a, b,
          method="spearman",
          use="pairwise.complete.obs")
}

#' Build an heatmap plot from expression samples
#'
#' @param data A #expression x #sample matrix.
#' @param distance.fn A function(a,b) returning the distance between a and b (default: 1-cor(method="spearman")).
#' @param mc.cores The number of cores.
#' @return an heatmap plot.
#' @examples
#' expression.heatmap(matrix(rnorm(64), nrow=8))
#' @export
expression.heatmap <- function (data,
                                distance.fn=distance.spearman,
                                mc.cores=4, ...)
{
    nsamples <- ncol(data)
    distMAT <- mat.or.vec(nsamples, nsamples)
    for (i in 1:(nsamples-1))
    {
        res <- unlist(
            parallel::mclapply((i+1):nsamples,
                               function (j) 
                               {
                                   distance.spearman(data[,i],
                                                     data[,j])
                               },
                               mc.cores=mc.cores)
        )
        distMAT[i,(i+1):nsamples] <- res
        distMAT[(i+1):nsamples,i] <- res
    }
    rownames(distMAT)=colnames(data)
    colnames(distMAT)=colnames(data)

    heatmap(distMAT, scale="none", ...)
}

