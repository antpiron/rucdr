
library(reshape2)
library(ggplot2)

defaultStepSize <- function (list1, list2)
{
    ceiling(min(sqrt( c(length(list1), length(list2)) )))
}	

numericListOverlap <- function(sample1, sample2, stepsize, alternative)
{
    n <- length(sample1)
    overlap <- function(a, b)
    {
        count <- as.integer(sum(as.numeric(sample1[1:a] %in% sample2[1:b])))
        the.mean <- a*b/n
        
        switch(alternative,
               enrichment={
                   pval <- phyper(q=count-1,
                                  m=a, n=n-a+1, k=b,
                                  lower.tail=FALSE)
                   fdr <- the.mean / count
                   signs <- 1L
               },
               two.sided={
                   signs <- sign(count - the.mean)
                   symmetric <- 2*the.mean - count
                   lower <- min(count, symmetric)
                   upper <- max(count, symmetric)

                   ## TODO: is this right? Think so.
                   fdr <- the.mean / upper
                   
                   pval <- phyper(q=lower,
                                  m=a, n=n-a+1, k=b,
                                  lower.tail=TRUE) +
                       phyper(q=upper,
                              m=a, n=n-a+1, k=b,
                              lower.tail=FALSE)
               })
        
        return( c(counts=count,
                  pval=as.numeric(pval),
                  signs=as.integer(signs),
                  fdr=as.numeric(fdr)) )    
    }
  
    logging(paste0("n = ", n, " ; stepsize = ", stepsize),
            .module="RRHO")
    indexes <- expand.grid(i=seq(1,n,by=stepsize),
                          j=seq(1,n,by=stepsize))
    overlaps <- apply(indexes, 1,
                      function(x) overlap(x['i'], x['j']))
    
    nrows<- sqrt(ncol(overlaps))
    matrix.counts <- matrix(overlaps['counts',], ncol=nrows)
    matrix.pvals <- matrix(overlaps['pval',], ncol=nrows)
    matrix.signs <- matrix(overlaps['signs',], ncol=nrows)  
    matrix.fdr <- matrix(overlaps['fdr',], ncol=nrows)  
  
    return(list(counts=matrix.counts,
                pval=matrix.pvals,
                signs=matrix.signs,
                fdr=matrix.fdr))  
}

#' @export RRHO
RRHO <- function(list1, list2, 
                 stepsize=NULL, 
                 alternative="two.sided")
{
    checkAction((alternative=='two.sided' || alternative=='enrichment'),
                .message="Alternative should be 'two.sided' or 'enrichment'.",
                .stop=T)

    if (is.null(stepsize))
        stepsize=defaultStepSize(list1, list2)

    ## Keep the common ids
    inter <- intersect(names(list1), names(list2))
    list1 <- list1[inter]
    list2 <- list2[inter]
        
    ## Order lists
    list1  <- list1[order(list1[inter],decreasing=TRUE)]
    list2  <- list2[order(list2[inter],decreasing=TRUE)]
    nlist1 <- length(list1)
    nlist2 <- length(list2)
    
    
    hypermat <- numericListOverlap(names(list1), names(list2),
                                   stepsize, alternative)
    
    hypermat$padj <- matrix(p.adjust(hypermat$pval, method="BY"),
                            nrow=nrow(hypermat$pval))

    return(structure(hypermat,
                     class="rrho"))
}

#' @export
plot <- function (rrho, ...)
{
    UseMethod("plot", rrho)
}

#' @export
plot.rrho <- function (rrho, min.pval=1e-12,
                       colors=c('darkblue', 'darkgreen',
                                'darkorange', 'darkred'))
{
    max.log <- -log(min.pval)
    no.zero <- apply(rrho$padj, 1:2,
                     function (p) min(p + min.pval, 1) )
    signed.log.pval <- -log(no.zero) * rrho$signs 

    b <- c(-max.log, 0, max.log)

    ggplot(data = melt(signed.log.pval),
           aes(x=Var1, y=Var2, fill=value)) + 
        geom_tile() +
        scale_fill_gradientn(colors = colors, breaks = b,
                             labels = format(b),
                             limits=b[c(1,length(colors))],
                             name="-log p.val") + 
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),)
}



## n <- 100
## list1 <- 1:n
## names(list1) <- 1:n
## list2 <- 1:n
## names(list2) <- 1:n
## rrho <- RRHO(list1,list2)
## plot(rrho)
