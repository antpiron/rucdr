
library(reshape2)
library(ggplot2)
library(ggrepel)

defaultStepSize <- function (list1, list2)
{
    ceiling(min(sqrt( c(length(list1), length(list2)) )))
}	

numericListOverlap <- function(sample1, sample2, stepsize, alternative, tol)
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
                   fdr <- the.mean / (count+1e-12)
                   signs <- 1L
               },
               two.sided={
                   signs <- sign(count - the.mean)
                   symmetric <- 2*the.mean - count
                   lower <- min(count, symmetric)
                   upper <- max(count, symmetric)

                   ## if(signs < 0){
                   ##     lower<- count 
                   ##     upper<- 2*the.mean - count 
                   ## } else{
                   ##     lower<- 2*the.mean - count 
                   ##     upper<- count 
                   ## }

                   ## TODO: is this right? Think so.
                   fdr <- the.mean / (upper+1e-12)

                   ## TODO: escape zeroes, cheat from RRHO orignal package, anything better?
                   ## hypergeom is expecting integer not double!?!!!!!!!
                   lower <- min(lower+tol, the.mean)
                   upper <- max(upper-tol, the.mean)
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
                 alternative="two.sided", tol=0.000001)
{
    checkAction((alternative=='two.sided' || alternative=='enrichment'),
                .message="Alternative should be 'two.sided' or 'enrichment'.",
                .stop=T)


    ## Keep the common ids
    inter <- intersect(names(list1), names(list2))
    list1 <- list1[inter]
    list2 <- list2[inter]
        
    ## Order lists
    list1  <- list1[order(list1[inter],decreasing=TRUE)]
    list2  <- list2[order(list2[inter],decreasing=TRUE)]
    nlist1 <- length(list1)
    nlist2 <- length(list2)

    if (is.null(stepsize))
        stepsize=defaultStepSize(list1, list2)
    
    hypermat <- numericListOverlap(names(list1), names(list2),
                                   stepsize, alternative, tol=tol)
    
    hypermat$padj <- matrix(p.adjust(hypermat$pval, method="BY"),
                            nrow=nrow(hypermat$pval))

    return(structure(c(hypermat, list(list1=list1, list2=list2),
                       nlist1=nlist1, nlist2=nlist2,
                       stepsize=stepsize),
                     class="rrho"))
}


getDown.index  <- function(l)
{
    ind <- which(l < 0)
    if (length(ind) > 0)
        ind[1]
    else
        length(l)+1
}

intersect.clean <- function(mat, x.bigger=TRUE, y.bigger=TRUE)
{
    op <- if (y.bigger) `<` else `>`
    sign.x <- sign(0.5-x.bigger)
    sign.y <- sign(0.5-y.bigger)
    
    rec <- function(mat)
    {
        if (nrow(mat) > 0)
        {
            higher <- intersect.clean(mat[op(mat[1,2], mat[,2]),, drop=FALSE])
            return(rbind(mat[1,], higher))
        }
        return(mat)
    }
    mat <- mat[order(sign.x*mat[,1], sign.y*mat[,2], decreasing=FALSE),,
                 drop=FALSE]
    
    return(rec(mat))
}

getUPDOWN.ij <-  function(rrho)
{
    if (is.null(rrho$best$updown))
    {
        ## -1 because half-brained started the R indexes at 1
        x.ind <- (getDown.index(rrho$list1) - 1) %/% rrho$stepsize
        y.len <- ncol(rrho$padj) - 1
        y.ind <- (getDown.index(rrho$list2)-1) %/% rrho$stepsize + 1
        logging(paste0("getUPDOWN.ij(): x = ", x.ind, " ; y = ", y.ind),
                .module="RRHO")
        ## start from 2 because phyper does not make any sense for 1
        UP.mat <- ifelse(rrho$sign[2:x.ind,y.ind:y.len,drop=FALSE] > 0, 1,
                         rrho$padj[2:x.ind,y.ind:y.len,drop=FALSE])
        updown <- t(t(which(UP.mat == min(UP.mat), arr.ind=TRUE)) + c(1,y.ind-1))
        updown <- intersect.clean(updown, y.bigger=FALSE)
        rrho$best$updown <- updown
    }
    
    return(rrho$best$updown)
}

#' @export
getUPDOWN <- function (rrho, ...)
{
    UseMethod("getUPDOWN", rrho)
}

#' @export
getUPDOWN.rrho <- function(rrho)
{
    min.ind  <- getUPDOWN.ij(rrho)
    logging(paste0("min.ind = ", min.ind[1,1], ", ", min.ind[1,2] ),
            .module="RRHO")
    as.vector(unique(unlist(apply(min.ind, 1,
            function (x)
                intersect(names(rrho$list1)[1:((x[1]-1) * rrho$stepsize) + 1],
                          names(rrho$list2)[((x[2]-1) * rrho$stepsize + 1):length(rrho$list2)])))))
}

getDOWNUP.ij <-  function(rrho)
{
    if (is.null(rrho$best$downup))
    {
        ## -1 because half-brained started the R indexes at 1
        x.ind <- (getDown.index(rrho$list1) - 1) %/% rrho$stepsize + 1
        x.len <- nrow(rrho$padj) - 1
        y.ind <- (getDown.index(rrho$list2)-1) %/% rrho$stepsize
        logging(paste0("getUPDOWN.ij(): x = ", x.ind, " ; y = ", y.ind),
                .module="RRHO")
        ## start from 2 because phyper does not make any sense for 1
        UP.mat <- ifelse(rrho$sign[x.ind:x.len,2:y.ind,drop=FALSE] > 0, 1,
                         rrho$padj[x.ind:x.len,2:y.ind,drop=FALSE])
        downup <- t(t(which(UP.mat == min(UP.mat), arr.ind=TRUE)) + c(x.ind-1,1))
        downup <- intersect.clean(downup, x.bigger=FALSE)
        rrho$best$downup <- downup
    }
    
    return(rrho$best$downup)
}

#' @export
getDOWNUP <- function (rrho, ...)
{
    UseMethod("getDOWNUP", rrho)
}

#' @export
getDOWNUP.rrho <- function(rrho)
{
    min.ind <- getDOWNUP.ij(rrho)
    logging(paste0("min.ind = ", min.ind[1,1], ", ", min.ind[1,2] ),
            .module="RRHO")
    as.vector(unique(unlist(apply(min.ind, 1,
            function (x)
                intersect(names(rrho$list1)[((x[1]-1) * rrho$stepsize + 1):length(rrho$list1)],
                           names(rrho$list2)[1:((x[2]-1) * rrho$stepsize) + 1])))))
}

getUPUP.ij <-  function(rrho)
{
    if (is.null(rrho$best$upup))
    {
        ## -1 because half-brained started the R indexes at 1
        x.ind <- (getDown.index(rrho$list1) - 1) %/% rrho$stepsize
        y.ind <- (getDown.index(rrho$list2) - 1) %/% rrho$stepsize
        logging(paste0("getUPUP.ij(): x = ", x.ind, " ; y = ", y.ind),
                .module="RRHO")
        ## start from 2 because phyper does not make any sense for 1
        UP.mat <- ifelse(rrho$sign[2:x.ind,2:y.ind,drop=FALSE] < 0, 1,
                         rrho$padj[2:x.ind,2:y.ind,drop=FALSE])
        upup <- which(UP.mat == min(UP.mat), arr.ind=TRUE) + 1
        upup <- intersect.clean(upup)
        rrho$best$upup <- upup
    }
    
    return(rrho$best$upup)
}

#' @export
getUPUP <- function (rrho, ...)
{
    UseMethod("getUPUP", rrho)
}

#' @export
getUPUP.rrho <- function(rrho)
{
    min.ind  <- getUPUP.ij(rrho)
    logging(paste0("min.ind = ", min.ind[1,1], ", ", min.ind[1,2] ),
            .module="RRHO")
    as.vector(unique(unlist(apply(min.ind, 1,
            function (x)
                intersect(names(rrho$list1)[1:((x[1]-1) * rrho$stepsize) + 1],
                          names(rrho$list2)[1:((x[2]-1) * rrho$stepsize) + 1])))))
}

getDOWNDOWN.ij <-  function(rrho)
{
    if (is.null(rrho$best$downdown))
    {
        ## TODO: remove last line, not needed if len(list)+1
        ## is not divisible by stepsize
        x.len <- nrow(rrho$padj) - 1
        ## -1 because half-brained started the R indexes at 1
        x.ind <- (getDown.index(rrho$list1)-1) %/% rrho$stepsize + 1
        y.len <- ncol(rrho$padj) - 1
        y.ind <- (getDown.index(rrho$list2)-1) %/% rrho$stepsize + 1
        logging(paste0("getDOWNDOWN.ij(): x = ", x.ind, " ; y = ", y.ind),
                .module="RRHO")
        ## 
        DOWN.mat <- ifelse(rrho$sign[x.ind:x.len,y.ind:y.len,drop=FALSE] < 0, 1,
                           rrho$padj[x.ind:x.len,y.ind:y.len,drop=FALSE])
        downdown <- t(t(which(DOWN.mat == min(DOWN.mat), arr.ind=TRUE)) +
                                c(x.ind-1,y.ind-1))
        downdown <- intersect.clean(downdown, x.bigger=FALSE, y.bigger=FALSE)
        rrho$best$downdown <- downdown
    }
    return(rrho$best$downdown)
}

#' @export
getDOWNDOWN <- function (rrho, ...)
{
    UseMethod("getDOWNDOWN", rrho)
}

#' @export
getDOWNDOWN.rrho <- function(rrho)
{
    min.ind <- getDOWNDOWN.ij(rrho)
    logging(paste0("min.ind = ", min.ind[1,1], ", ", min.ind[1,2] ),
            .module="RRHO")
    ##print(rrho$nlist1)
    ##print(min.ind)
    result <- as.vector(unique(unlist(apply(min.ind, 1,
                             function (x)
                                 intersect(names(rrho$list1)[((x[1]-1) * rrho$stepsize + 1):length(rrho$list1)],
                                           names(rrho$list2)[((x[2]-1) * rrho$stepsize + 1):length(rrho$list2)])))))

    ##print("END")
    
    return(result)
}

#' @export
getEnrichment <- function (rrho, ...)
{
    UseMethod("getEnrichment", rrho)
}

#' @export
getEnrichment.rrho <- function (rrho, direction="upup")
{
    res <- switch(direction,
                  upup={
                      getUPUP(rrho)
                  },
                  downdown={
                      getDOWNDOWN(rrho)
                  },
                  updown={
                      getUPDOWN(rrho)
                  },
                  downup={
                      getDOWNUP(rrho)
                  })
    return(res)
}




#' @export
plot <- function (rrho, ...)
{
    UseMethod("plot", rrho)
}

my.colors <- c('darkblue', 'darkorchid4', 'darkgreen', 'lightpink4',
##                'tan',
               'darkorange', 'chocolate4', 'brown4',  'darkred')
## my.colors <- c('darkblue',"cyan4","cyan3","cyan2","cyan1","cyan",
##                "palegreen4","palegreen3","palegreen2","palegreen1",
##                "yellow",
##                "orange1","orange2","orange3","orange4",
##                "red", "red1", "red2", "red3", "red4", 'darkred')


## my.colors <- c(sapply(0:255,
##                       function (x)
##                           paste0("#",format(as.hexmode(255-x),width=2),
##                                   format(as.hexmode(x),width=2), "00")),
##               sapply(0:255,
##                       function (x)
##                           paste0("#00",format(as.hexmode(255-x),width=2),
##                                  format(as.hexmode(x),width=2))))


signed.lpval <- function (pval, signs, .log=log2)
{
    no.zero <- ifelse(pval == 0, 10, pval)
    min.pval <- min(no.zero) / 2
    no.zero <- ifelse(no.zero == 10, min.pval, no.zero)
    return(-.log(no.zero) *signs)
}

#' @export
plot.rrho <- function (rrho,
                       colors=my.colors,
                       labels=c("",""),
                       .log=log2, repel.force=150)
{
    signed.log.pval <- signed.lpval(rrho$padj, rrho$signs, .log)
    max.log <- max(abs(signed.log.pval))
    min.log <- - max.log

    b <- c(min.log, 0, max.log)
    text_up <- grid::textGrob("up", gp=grid::gpar(fontsize=13, fontface="bold"))
    text_down <- grid::textGrob("down", gp=grid::gpar(fontsize=13, fontface="bold"))
    text_up_rot <- grid::textGrob("up", rot=90,
                               gp=grid::gpar(fontsize=13, fontface="bold"))
    text_down_rot <- grid::textGrob("down", rot=90,
                              gp=grid::gpar(fontsize=13, fontface="bold"))
    vperc <- ncol(rrho$pval) / 10
    hperc <- nrow(rrho$pval) / 10
    upup <- as.data.frame(getUPUP.ij(rrho))
    upup$padj <- apply(upup,1,function (x) rrho$padj[x[1],x[2]])
    downdown <-  as.data.frame(getDOWNDOWN.ij(rrho))
    downdown$padj <- apply(downdown,1,function (x) rrho$padj[x[1],x[2]])
    updown <-  as.data.frame(getUPDOWN.ij(rrho))
    updown$padj <- apply(updown,1,function (x) rrho$padj[x[1],x[2]])
    downup <-  as.data.frame(getDOWNUP.ij(rrho))
    downup$padj <- apply(downup,1,function (x) rrho$padj[x[1],x[2]])
    len.colors <- length(colors)
    half.len.colors <- len.colors %/% 2
    colors.values <- seq(0, len.colors) /  len.colors
    ##colors.values <- ifelse(colors.values > 0.5,
    ##                        2^(colors.values+0.5)-2,
    ##                        log2(colors.values+1))
    ##print(colors.values)

    gg <- ggplot2::ggplot() + 
        ggplot2::geom_tile(data = melt(signed.log.pval),
                           aes(x=Var1, y=Var2, fill=value)) +
        ggplot2::scale_fill_gradientn(colors = colors, breaks = b,
                                      labels = format(b),
                                      limits=b[c(1,3)],
                                      ##limits=b[c(1,length(colors))],
                                      name="-log p.val",
                                      values=colors.values) + 
        ggplot2::theme(##axis.title.x=element_blank(),
                     axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank(),
                     ##axis.title.y=element_blank(),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank()) +
        ggplot2::xlab(labels[1]) +  ggplot2::ylab(labels[2]) +
        ggplot2::annotation_custom(text_up,
                                   xmin=hperc,xmax=hperc,
                                   ymin=-hperc/4,ymax=-hperc/4) +
        ggplot2::annotation_custom(text_down,
                                   xmin=nrow(rrho$pval)+hperc,
                                   xmax=nrow(rrho$pval)-3*hperc,
                                   ymin=-hperc/4,ymax=-hperc/4) +
        ggplot2::annotation_custom(text_up_rot,
                                   xmin=-vperc/4,xmax=-vperc/4,
                                   ymin=vperc,ymax=vperc) +
        ggplot2::annotation_custom(text_down_rot,
                                   ymin=ncol(rrho$pval)+vperc,
                                   ymax=ncol(rrho$pval)-3*vperc,
                                   xmin=-vperc/4,xmax=-vperc/4) 
    x.ind <- which(rrho$list1 < 0)
    if ( length(x.ind) > 0 )
        gg  <- gg +
            ggplot2::geom_vline(aes(xintercept = (x.ind[1]-1) / rrho$stepsize + 1), 
                       linetype = "dotted", colour = "gray10",size = 1)
    y.ind <- which(rrho$list2 < 0)
    if ( length(y.ind) > 0 )
        gg  <- gg +
            ggplot2::geom_hline(aes(yintercept = (y.ind[1]-1) / rrho$stepsize + 1), 
                                linetype = "dotted", colour = "gray10",size = 1)

    gg <- gg +
        ggrepel::geom_text_repel(data=upup,
                           aes(x=row, y=col,
                               label=formatC(padj,
                                             format = "e", digits = 1)),
                           hjust=1, vjust=1,
                           color="black",
                           force = repel.force, show.legend = FALSE) +
        ggrepel::geom_text_repel(data=downdown,
                           aes(x=row, y=col,
                               label=formatC(padj,
                                             format = "e", digits = 1)),
                           hjust=0, vjust=0,
                           color="black",
                           force = repel.force, show.legend = FALSE) +
        ggrepel::geom_text_repel(data=updown,
                           aes(x=row, y=col,
                               label=formatC(padj,
                                             format = "e", digits = 1)),
                           hjust=1, vjust=0,
                           color="black",
                           force = repel.force, show.legend = FALSE) +
        ggrepel::geom_text_repel(data=downup,
                           aes(x=row, y=col,
                               label=formatC(padj,
                                             format = "e", digits = 1)),
                           hjust=0, vjust=1,
                           color="black",
                           force = repel.force, show.legend = FALSE) +
        ggplot2::geom_point(data=upup,
                            aes(x=row, y=col)) +
        ggplot2::geom_point(data=downdown,
                            aes(x=row, y=col)) +
        ggplot2::geom_point(data=updown,
                            aes(x=row, y=col)) +
        ggplot2::geom_point(data=downup,
                            aes(x=row, y=col))

    return(gg)
}



myColorRamp <- function(colors, values)
{
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

#' @export
plot3d <- function (rrho,
                    labels=c("x","y"),
                    .log=log2)
{
    padj <- rrho$padj[2:nrow(rrho$padj),2:ncol(rrho$padj)]
    signs <- rrho$signs[2:nrow(rrho$signs),2:ncol(rrho$signs)]
    signed.log.pval <- signed.lpval(padj, signs, .log)
    max.log <- max(abs(signed.log.pval))
    min.log <- - max.log
    melted <- melt(signed.log.pval)
    colnames(melted) <- c(labels, "-log pval")
    ## print(head(melted))
    colors <- myColorRamp(c("blue", "red"), melted[,3]) 
    scatterplot3d::scatterplot3d(melt(signed.log.pval),
                                 color=colors,
                                 xlab=labels[1],
                                 ylab=labels[2],
                                 zlab="-log pval")
}

## n <- 100
## list1 <- 1:n
## names(list1) <- 1:n
## list2 <- 1:n
## names(list2) <- 1:n
## rrho <- RRHO(list1,list2)
## plot(rrho)
