require(survival)
require(multtest)

aCGH.test.old <-
    function(formula, aCGH.obj, test = c("t.test", "coxph"),
             grouping = NULL, p.adjust.method = "fdr", subset = NULL)
{

#    if(missing(formula) || !inherits(formula, "formula"))
#        stop("formula missing or invalid")
#        m <- match.call(expand.dots = FALSE)
#    if(is.matrix(eval(m$data, parent.frame())))
#        m$data <- as.data.frame(data)
#    m[[1]] <- as.name("model.frame")
#    m$... <- NULL
#    mf <- eval(m, parent.frame())
#    DNAME <- paste(names(mf), collapse = " by ")
#    names(mf) <- NULL
#    response <- attr(attr(mf, "terms"), "response")
#    g <- factor(mf[[-response]])
#    if(nlevels(g) != 2)
#        stop("grouping factor must have exactly 2 levels")
#    DATA <- split(mf[[response]], g)
#    browser()

    pheno <- phenotype(aCGH.obj)
    if (!is.null(subset))
        pheno <- pheno[ subset, ]
    test <- match.arg(test)
    switch(test,
           t.test = t.test(formula, pheno),
           coxph = coxph(formula, pheno)
           )
#    do.call(test, formula, data = phenotype(aCGH.obj))
#    invisible()
    
}

aCGH.test <-
    function(frml, aCGH.obj, test = c("survdiff", "coxph"),
             grouping = NULL, p.adjust.method = "fdr", subset = NULL)
{

    l2r <- log2.ratios.imputed(aCGH.obj)
    if (!is.null(subset))
        l2r <- l2r[ subset, ]
    test <- match.arg(test)
    pheno <- phenotype(aCGH.obj)
    resT <- 
        sapply(1:nrow(l2r),
               function(i) {
                   
                   if (i %% 100 == 0)
                       print(i)
                   clone <- l2r[ i, ]
                   switch(test,
                          survdiff = {
                              
                              survdiff.fit <-
                                  try(survdiff(as.formula(frml),
                                               data = pheno))
                              if (inherits(survdiff.fit, "try-error"))
                                  c(0, 1)
                              else
                              {
                                  
                                  etmp <- 
                                      if (is.matrix(survdiff.fit$obs))
                                          apply(survdiff.fit$exp,
                                                1,
                                                sum)
                                      else
                                          survdiff.fit$exp
                                  c(survdiff.fit$chisq,
                                    1 - pchisq(survdiff$chisq,
                                               sum(etmp > 0) - 1)
                                    )
                                  
                              }
                              
                          },
                          coxph = {

                              coxph.fit <-
                                  try(coxph(as.formula(frml),
                                            data = pheno))
                              if (inherits(coxph.fit, "try-error"))
                                  c(0, 1)
                              else
                              {
                                  
                                  logtest <-
                                      -2 * (coxph.fit$loglik[1] -
                                            coxph.fit$loglik[2])
                                  beta <- coxph.fit$coef
                                  df <- length(beta[!is.na(beta)])
                                  c(logtest, 1 - pchisq(logtest, df))
                                  
                              }
                              
                          }
                          )
                   
               }
               )
    rawp <- resT[ 2, ]
    adjp <- p.adjust(rawp, p.adjust.method)

    data.frame(index = 1:ncol(resT),
               teststat = resT[ 1, ],
               rawp = rawp,
               adjp = adjp
               )[ order(adjp), ]

}

threshold.func <-
    function(dat, thresAbs)
{
    
    out <- matrix(0, nrow = nrow(dat), ncol = ncol(dat))
    ##if the same threshold for all samples
    if (length(thresAbs) == 1)
        thresAbs <- rep(thresAbs, ncol(dat))
    if (length(thresAbs) != ncol(dat))
        stop("Error: number of threshold is not the same as number of\
samples")
##    for (i in 1:ncol(dat))
##    {
##        tmp <- dat[,i]
##        tmp[dat[,i] >=thresAbs[i] & !is.na(dat[,i])] <- 1
##        tmp[dat[,i] <=-thresAbs[i] & !is.na(dat[,i])] <- -1
##        tmp[dat[,i] > -thresAbs[i] & dat[,i] < thresAbs[i] & !is.na(dat[,i])] <- 0
        
##        out[,i] <- tmp
##    }
##    out
    sapply(1:ncol(dat),
           function(i) {
               
               tmp <- rep(0, ncol(dat))
               na.col <- is.na(dat[ ,i ])
               col <- dat[ ,i ][!na.col]
               tmp[!na.col] <-
                   ifelse(col > thresAbs[i],
                          1,
                          ifelse(col < thresAbs[i], -1, 0)
                          )
               
               tmp
               
           }
           )
    
}

changeProp.func <-
    function(dat = data.screen.norm.thres, colMatr)
{
    
    out <- matrix(0, nrow = nrow(dat), ncol = nrow(colMatr))
    for (i in 1:nrow(colMatr))
        for (j in 1:nrow(dat))
        {
            
            vec <- dat[j, colMatr[ i, ] == 1]
            out[j, i] <- 
                if (lengthLoss.na(vec) < lengthGain.na(vec))
                    propGain.na(vec)
                else
                    propLoss.na(vec)
            
        }
    out
    
}

changeProp.overall.func <-
    function(dat)
    apply(dat,
          1,
          function(vec) {
               
              if (lengthLoss.na(vec) < lengthGain.na(vec))
                  propGain.na(vec)
              else
                  propLoss.na(vec)
              
          }
          )

table.bac.func <-
    function(dat, colMatr)
{

    nr <- 
        if (nrow(colMatr) > 1)
            6 * (nrow(colMatr) + 1)
        else
            6
    out <- matrix(0, nrow=nrow(dat), ncol=nr)	
    
    ##number of samples:

    sample.ind <- which(colMatr[ 1, ] == 1)
    if (nrow(colMatr) > 1)
        for (j in 2:nrow(colMatr))
            sample.ind <- c(sample.ind, which(colMatr[ j, ] == 1))

    ##all samples 
    sample.ind <- unique(sample.ind)
    
    dat.all <- dat[ ,sample.ind ]
    
    for (i in 1:nrow(dat))
    {
        
        len <- length(sample.ind)
        vec <- dat.all[ i, ]
        out[i, 1] <- sum(!is.na(vec)) #number present
        out[i, 2] <- lengthGain.na(vec) #number gained
        out[i, 3] <- lengthLoss.na(vec) #number lost
        out[i, 4] <- round(1 - prop.na(vec),2) #proportion present 
        out[i, 5] <- round(propGain.na(vec),2) #proportion gained
        out[i, 6] <- round(propLoss.na(vec),2) #proportion lost
	
        if (nr > 6) #if 2 or more classes
        {
            
            for (j in 1:nrow(colMatr))
            {
                
                vec <- dat[i, colMatr[ j, ] == 1]
                out[i, (6*j+1):(6*j+6)] <-
                    c(sum(!is.na(vec)), lengthGain.na(vec),
                      lengthLoss.na(vec), round(1 - prop.na(vec), 2),
                      round(propGain.na(vec), 2),
                      round(propLoss.na(vec), 2)
                      )
                
            }
            
        }
        
    }
    out

}

lengthGain.na <-
    function(x)
    sum(x == 1, na.rm = TRUE)

propGain.na <-
    function(x)
    mean(x == 1, na.rm = TRUE)

lengthLoss.na <-
    function(x)
    sum(x == -1, na.rm = TRUE)

propLoss.na <-
    function(x)
    mean(x == -1, na.rm = TRUE)

prop.na <-
    function(x)
    mean(is.na(x))

gainloss.func <-
    function (dat, cols, thres, quant = .5)
{

    if (length(thres) == 1)
        thres <- rep(thres, ncol(dat))
    if (length(thres) != ncol(dat))
        stop("Error: number of thresholds is not the same as number\
of samples")
    dt <- as.matrix(dat[ ,cols ])
    thr <- thres[cols]
    loss.med <- loss <- gain.med <- gain <- rep(0, nrow(dt))
    
    for (i in 1:nrow(dt))
        if (!all(is.na(dt[ i, ])))
        {
            
            x <- dt[ i, ]
            th <- thr[!is.na(x)]
            x <- x[!is.na(x)]
            tmp.gain <- x >= th
            gain[i] <- mean(tmp.gain)
            if (any(tmp.gain))
                gain.med[i] <- quantile(x[tmp.gain], 1 - quant)
            tmp.loss <- x <= -th
            loss[i] <- mean(tmp.loss)
            if (any(tmp.loss))
                loss.med[i] <- quantile(x[tmp.loss], quant)
            
        }
    
    list(gainP = gain,
         lossP = loss,
         gainMed = gain.med,
         lossMed = loss.med)
    
}

plotFreqStatColors <-
    function(aCGH.batch, resT, pheno, colored = TRUE, ...)
    plotfreq.stat(aCGH.batch, resT, pheno, colored = TRUE, ...)

plotFreqStatGrey <-
    function(aCGH.batch, resT, pheno, colored = FALSE, ...)
    plotfreq.stat(aCGH.batch, resT, pheno, colored = FALSE, ...)

plotFreqStat <-
    function(aCGH.obj, resT, pheno, chrominfo = human.chrom.info.Jul03,
             X = TRUE, Y = FALSE, threshold = TRUE, minChanged = 0, all = FALSE,
             rsp.uniq = unique(pheno), nlim = 1, cutplot = 0,
             titles = rsp.uniq, thres = .2, ylm = c(-1, 1),
             ngrid = 2, p.thres = c(.01, .05, .1), mincolors = .1,
             quant.col = .11, numaut = 22, onepage = TRUE, colored = TRUE,
             summarize.clones = TRUE)
{

    col.scheme <- 
        if (colored)
            list(pal =
                 c("red", "blue", "green",
                   "yellow")[1:length(p.thres)],
                 gain.low = "white",
                 gain.high = "green",
                 loss.low = "red",
                 loss.high = "white",
                 abline1 = "blue",
                 abline2 = "grey50",
                 mtext = "red",
                 kb.loc = "blue",
                 abline3 = "black",
                 abline4 = "grey50",
                 )
        else
            list(pal =
                 c("grey10", "grey40", "grey70",
                   "grey90")[1:length(p.thres)],
                 gain.low = "grey50",
                 gain.high = "grey0",
                 loss.low = "grey0",
                 loss.high = "grey50",
                 abline1 = "grey50",
                 abline2 = "grey50",
                 mtext = "black",
                 kb.loc = "black",
                 abline3 = "black",
                 abline4 = "grey50",
                 )
    data <- log2.ratios(aCGH.obj)
    datainfo <- clones.info(aCGH.obj)
    dataSign <- log2.ratios.imputed(aCGH.obj)
    rsp.uniq <- sort(rsp.uniq)
    
    ## creating response matrix colmatr

    if (!all)
	colmatr <-
            t(
              sapply(rsp.uniq,
                     function(rsp.uniq.level)
                     ifelse(pheno == rsp.uniq.level, 1, 0)
                     )
              )

    ## screening out clones that are gained or lost in < minChanged in
    ## classes under comparison
    ## indeces present:

    tmp <- apply(as.matrix(colmatr), 2, sum)
    indecesnow <- which(tmp == 1)
    data.thres <- threshold.func(data, thresAbs = thres)
    prop.ch <- changeProp.func(dat = data.thres, colMatr = colmatr)
    maxch <- changeProp.overall.func(dat = data.thres[ ,indecesnow ])
    clones.index <- which(maxch >= minChanged)

    ##removing clones to skip from the dataset

    data <- data[clones.index,]
    data.thres <- data.thres[clones.index,]
    dataSign <- dataSign[clones.index,]
    datainfo <- datainfo[clones.index,]

    ## start table:
    if (summarize.clones)
        bac.summary <-
            table.bac.func(dat = data.thres, colMatr = colmatr)

    ## creating color matrix for displaying intensity of gains and losses
    ## and for plotting p-values

    colors.gain <-
        maPalette(low = col.scheme$gain.low,
                  high = col.scheme$gain.high,
                  k = ngrid)
    colors.loss <-
        maPalette(low = col.scheme$loss.low,
                  high = col.scheme$loss.high,
                  k = ngrid)
###    sq.loss <- seq(-nlim, -mincolors, length = ngrid + 1)
###    sq.gain <- seq(mincolors, nlim, length = ngrid + 1)
###    matr.colors.loss <-
###        data.frame(sq.loss[ -length(sq.loss) ], sq.loss[-1],
###                   colors.loss)
###    matr.colors.gain <-
###        data.frame(sq.gain[ -length(sq.gain) ], sq.gain[-1],
###                   colors.gain)

    ## Now, start:
    ## if perform significance analysis on thresholded data only:

    if (threshold)
	dataSign <- threshold.func(dataSign, thres)
    nr <- nrow(colmatr)
    
    ## if 1 class only, no significance analysis:
###    if (nrow(colmatr) == 1)
###        sign <- F
###    if (sign)
###        nr <- nr + 1
    nr <- nr + 1
    
    tmp <- as.data.frame(matrix(0, ncol = 2, nrow = 1))
    colnames(tmp) <- c("gainP", "lossP")
    gainloss <-
        lapply(1:nrow(colmatr),
               function(j)
               gainloss.func(dat = data,
                             cols = which(colmatr[ j, ] == 1),
                             thres = thres,
                             quant = quant.col)
               )
    dt <- dataSign[ ,colmatr[1,] == 1, drop = FALSE ]
    rsp <- rep(1, ncol(dt))
    for (j in 2:nrow(colmatr))
    {
        
        dt <- cbind(dt, dataSign[ ,colmatr[ j, ] == 1 ])
        rsp <- c(rsp, rep(j, sum(colmatr[ j, ] == 1)))
        
    }
    rsp <- rsp - 1

    ## Process statistics
    ## for plotting test stats and p-values

    res <- resT[clones.index,]
    maxT <- res$adjp[order(res$index)]
    
    teststat <- abs(res$teststat[order(res$index)])
    st.now <-
        sapply(p.thres,
               function(threshold) {
                   
                   if (any(maxT <= threshold))
                       min(teststat[maxT <= threshold])
                   else
                       NA
                   
               }
               )
    pal.now <- col.scheme$pal

    ##append to bac summary file
    if (summarize.clones)
        bac.summary <-
            cbind(bac.summary, res$rawp[order(res$index)], maxT)
    
    ##Now preparing for plotting:

    numchr <- numaut
    if (X)
        numchr <- numchr + 1
    if (Y)
        numchr <- numchr + 1
    chrominfo <- chrominfo[ 1:numchr, ]

    ##compute cumulative kb locations
    start <- c(0, cumsum(chrominfo$length))
    kb.loc <- datainfo$kb
    for (i in 1:nrow(chrominfo))
        kb.loc[datainfo$Chrom == i] <-
            start[i] + datainfo$kb[datainfo$Chrom == i]

    ## preparation for graphs
    chrom.start <- c(0, cumsum(chrominfo$length))[1:numchr]
    chrom.centr <- chrom.start + chrominfo$centr
    chrom.mid <- chrom.start + chrominfo$length / 2

    ##now, plot
    par(mfrow = c((if (onepage) nr else 1), 1), lab = c(1, 8, 7),
        tcl = -.2,  xaxs = "i")

    for (g in 1:length(titles))
    {

        gl <- gainloss[[g]]
        tl <- as.character(titles[g])
        ylm[1] <- min(ylm, min(gl$lossP))
        ylm[2] <- max(ylm, max(gl$gainP))

###        col.nrow <-
###            sapply(gl$gainMed,
###                   function(cl) {
                       
###                       if (cl >= nlim)
###                           cl <- nlim - 10 ^ (-6)
###                       cnr <- 
###                           which(cl >= matr.colors.gain[ ,1 ] &
###                                 cl < matr.colors.gain[ ,2 ])
###                       if (length(cnr) > 0)
###                           cnr
###                       else
###                           1
                       
###                   }
###                   )
        ind <- which(gl$gainP >= cutplot)
        plot(kb.loc[ind], gl$gainP[ind],
             col = "green",
###             as.character(matr.colors.gain[ind, 3][col.nrow[ind]]),
             type = "h", xlab = "chromosome number",
             ylab = "Fraction gained or lost", pch = 18, main = tl,
             ylim = ylm,
             xlim = c(0, max(cumsum(chrominfo$length), kb.loc[ind],
             rm.na = TRUE))
             )
###        col.nrow <-
###            sapply(gl$lossMed,
###                   function(cl) {
                       
###                       if (cl <=- nlim)
###                           cl <- -nlim + 10 ^ (-6)
###                       cnr <-
###                           which(cl >= matr.colors.loss[ ,1 ] &
###                                 cl < matr.colors.loss[ ,2 ])
###                       if (length(cnr) > 0)
###                           cnr
###                       else
###                           ngrid
                       
###                   }
###                   )
        ind <- gl$lossP >= cutplot
        points(kb.loc[ind], -gl$lossP[ind],
               col = "red",
###               as.character(matr.colors.loss[ind, 3][col.nrow[ind]]),
               type = "h")

        abline(h = 0)
        abline(v = cumsum(chrominfo$length), col = col.scheme$abline1)
        abline(v = chrom.centr, lty = 2, col = col.scheme$abline2)

        for (i in seq(2, numaut, b = 2))
            mtext(paste("", i), side = 3, at = (chrom.mid[i]),
                  line = .3, col = col.scheme$mtext, cex.main = .5)
        for (i in seq(1, numaut, b = 2))
            mtext(paste("", i), side = 1, at = (chrom.mid[i]),
                  line = .3, col = col.scheme$mtext, cex.main = .5)
        if (X)
            if (i == numaut)
                mtext("X", side = 1, at = (chrom.mid[numaut + 1]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
            else
                mtext("X", side = 3, at = (chrom.mid[numaut + 1]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
        if (Y)
            if (i == numaut)
                mtext("Y", side = 3, at = (chrom.mid[numaut + 2]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
            else
                mtext("Y", side = 1, at = (chrom.mid[numaut + 2]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
        
    }
    plot(kb.loc, teststat, col = col.scheme$kb.loc,
         ylim = c(0, max(teststat)), type = "h",
         xlab = "chromosome number", ylab = "clone statistic",
         pch = 18, main = paste(titles, collapse = " vs "),
         xlim = c(0, max(cumsum(chrominfo$length), kb.loc, rm.na = TRUE))
         )
    if (length(st.now) > 0)
        abline(h = rev(st.now), col = rev(pal.now), lty = 2)
    abline(v = cumsum(chrominfo$length), col = col.scheme$abline3)
    abline(v = chrom.centr, lty = 2, col = col.scheme$abline4)

    for (i in seq(1, numaut, b = 2))
        mtext(paste("", i), side = 1, at = chrom.mid[i], line = .3,
              col = col.scheme$mtext, cex.main = .5)
    for (i in seq(2, numaut, b = 2))
        mtext(paste("", i), side = 3, at = chrom.mid[i], line = .3,
              col = col.scheme$mtext, cex.main = .5)
    if (X)
        if (i == numaut)
            mtext("X", side = 1, at = chrom.mid[numaut + 1],
                  line = .3, col = col.scheme$mtext, cex.main = .5)
        else
            mtext("X", side = 3, at = chrom.mid[numaut + 1],
                  line = .3, col = col.scheme$mtext, cex.main = .5)
    if (Y)
        if (i == numaut)
            mtext("Y", side = 3, at = chrom.mid[numaut + 2],
                  line = .3, col = col.scheme$mtext, cex.main = .5)
        else
            mtext("Y", side = 1, at = chrom.mid[numaut + 2],
                  line = .3, col = col.scheme$mtext, cex.main = .5)

    if (summarize.clones)
    {
        
        bac.summary <- as.data.frame(bac.summary)
        nms <-
            c("NumPresent", "NumGain", "NumLost", "PropPresent",
              "PropGain", "PropLost")
        cnames <- colnames(bac.summary)
        cnames[1:6] <- paste(nms, "All", sep = ".")
        if (nrow(colmatr) > 1)
            for (m in 1:length(rsp.uniq))
                cnames[ (6 * m + 1):(6 * (m + 1)) ] <-
                    paste(nms, rsp.uniq[m], sep = ".")
        cnames[(ncol(bac.summary)-1):ncol(bac.summary)] <-
            c("rawp", "adjp.maxT")
        colnames(bac.summary) <- cnames
        bac.summary <- cbind(datainfo, bac.summary)	
###    write.table(bac.summary, filetable, col.names = TRUE, row.names = FALSE,
###                sep = "\t", quote = FALSE)
        invisible(bac.summary)
        
    }
    
}

## This description is old!
##frequency plot for the whole genome using outside p-values and stats

##data
##rsp -- phenotype, NA are not allowed, have to be consequetive integers
##datainfo
##chrominfo
##titles -- the titles of the frequency plots -- has to have as many names as 
##levels in the response
##thres -- unique threshold or vector of tumor specific thresholds. In the latter
##case has to contain as many thershold as samples
##cutplot -- don't plots clones which gained/lost in fewer than <= fraction of cases
##sign -- to do significance comparison (T) or not (F). to do comparison uses
##multtest package
##nperm -- if sign =T, then how many permutations for maxT
##test -- name of the test
##ranks -- whether to work with ranked data If nor, "n"
##side -- two sisded (abs) test or 1-sided ("upper" or "lower")
##p.thres -- for which adjusted p-values show the cut-off
##filePS -- name of the ps/pdf file
##PS = "ps" or "pdf"
##X=T -- is X (23) chrom included?
##Y=F  -- is Y (24) chrom included?
##numaut=22 -- number of autosomes
##ngrid=50: density of colors
##nlim=1: upper limit for solors
##mincolors: minimum limit for colors
##by defult "white"corersponds to [-.2,.2] and red and green to [-1,-.2] and [.2,1[]
##respectively
## quant.col=.5: percentile for coloring gaind/lost clones -- <= .5, E.g. .25
##would correspond to taking 25th % for lost and 75% for gained samples

plotFreqGivenStat <-
    function(aCGH.obj, stat, statPerm, pheno, summarize.clones = FALSE,
             ...)
{

    maxstat <- apply(statPerm, 2, max, na.rm = TRUE)
    plotFreqStat(aCGH.obj,
                 resT =
                 list(teststat = stat,
                      adjp =
                      sapply(teststat,
                             function(t.i) sum(maxstat >= t.i)
                             ) / length(maxstat)
                      ),
                 pheno,
                 summarize.clones = FALSE,
                 ...)
    
}

plotfreqGivenStatFinalColors <-
    function(aCGH.obj, ...)
    plotfreq.givenstat.final.colors.func(data =
                                         log2.ratios(aCGH.obj),
                                         datainfo =
                                         clones.info(aCGH.obj),
                                         ...)
