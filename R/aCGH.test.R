create.resT <-
    function(resT.raw, p.adjust.method = "fdr")
{
    
    rawp <- resT.raw[ 2, ]
    adjp <- p.adjust(rawp, p.adjust.method)
    teststat <- resT.raw[ 1, ]
    
    data.frame(index = 1:ncol(resT.raw),
               teststat = teststat,
               rawp = rawp,
               adjp = adjp
               )[ order(adjp, rawp, teststat), ]

}

aCGH.test <- function(aCGH.obj, rsp, test = c("survdiff", "coxph", "linear.regression"),p.adjust.method = "fdr", subset = NULL, strt = NULL, ...)
{

    l2r <- as.matrix(log2.ratios.imputed(aCGH.obj))
    if (!is.null(subset))
        l2r <- l2r[ subset, ]
    test <- match.arg(test)
    pheno <- phenotype(aCGH.obj)
    resT <- 
        sapply(1:nrow(l2r),
               function(i) {
                   
###                   if (i %% 100 == 0)
###                       print(i)
                   clone <- l2r[ i, ]
                   fmla <-
                       if (!is.null(strt))
                           rsp ~ clone + strata(strt)
                       else
                           rsp ~ clone
                   switch(test,
                          survdiff = {
                              
                              survdiff.fit <- try(survdiff(fmla, ...))
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
                                               sum(etmp > 0) - 1))
                                  
                              }
                              
                          },
                          coxph = {

                              coxph.fit <- try(coxph(fmla, ...))
                              if (inherits(coxph.fit, "try-error"))
                                  c(0, 1)
                              else
                              {
                                  cf <-  coxph.fit$coef
				  cf.se <- sqrt(coxph.fit$var)
				  cf.std <- cf/cf.se
				  c(cf.std, 2*(1-pnorm(abs(cf.std))))
                                  #logtest <-
                                  #    -2 * (coxph.fit$loglik[1] -
                                  #          coxph.fit$loglik[2])
                                  #beta <- coxph.fit$coef
                                  #df <- length(beta[!is.na(beta)])
                                  #c(logtest, 1 - pchisq(logtest, df))
                                  
                              }
                              
                          },
                          linear.regression = {
                              
				reg <- lm(fmla, ...)
				cf <- (summary(reg))$coef
				c(cf[2,3], cf[2,4])
                        ##      fstat <-
                        ##       summary(lm(clone ~ rsp, ...))$fstatistic
                        ##   c(fstat[1],
                        ##    1 - pf(fstat[1], fstat[2], fstat[3]))
                              
                          }
###                          logistic.regression = {
###
###                              glm.fit <-
###                                  try(glm(fmla, family=binomial()))
###                              if (inherits(glm.fit, "try-error"))
###                                  c(0, 1)
###                              else
###                              {
###                                  
###                                  stat <-
###                                      2 * (glm.fit$null.deviance -
###                                           glm.fit$deviance)
###                                  c(stat,
###                                    1 - pchisq(stat,
###                                               glm.fit$df.null -
###                                               glm.fit$df.residual))
###
###                              }
###                              
###                          }
                          )
                   
               }
               )
    
    create.resT(resT, p.adjust.method)
    
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
               
               tmp <- rep(NA, nrow(dat))
               na.col <- is.na(dat[ ,i ])
               col <- dat[ ,i ][!na.col]
               tmp[!na.col] <-
                   ifelse(col >= thresAbs[i],
                          1,
                          ifelse(col <= -thresAbs[i], -1, 0)
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

gainLoss <-
    function (dat, cols, thres = 0.25)
{

    if (length(thres) == 1)
        thres <- rep(thres, ncol(dat))
    if (length(thres) != ncol(dat))
        stop("Error: number of thresholds is not the same as number\
of samples")
    dt <- as.matrix(dat[ ,cols ])
    thr <- thres[cols]
    loss <- gain <- rep(0, nrow(dt))
    
    for (i in 1:nrow(dt))
        if (!all(is.na(dt[ i, ])))
        {
            
            x <- dt[ i, ]
            th <- thr[!is.na(x)]
            x <- x[!is.na(x)]
            tmp.gain <- x >= th
            gain[i] <- mean(tmp.gain)
            #if (any(tmp.gain))
            #    gain.med[i] <- quantile(x[tmp.gain], 1 - quant)
            tmp.loss <- x <= -th
            loss[i] <- mean(tmp.loss)
            #if (any(tmp.loss))
            #    loss.med[i] <- quantile(x[tmp.loss], quant)
            
        }
    
    list(gainP = gain,lossP = loss)
       
    
}

plotFreqStatColors <-
    function(aCGH.batch, resT, pheno, colored = TRUE, ...)
    plotfreq.stat(aCGH.batch, resT, pheno, colored = TRUE, ...)

plotFreqStatGrey <-
    function(aCGH.batch, resT, pheno, colored = FALSE, ...)
    plotfreq.stat(aCGH.batch, resT, pheno, colored = FALSE, ...)

plotFreqStat <-
    function(aCGH.obj, resT = NULL, pheno = rep(1, ncol(aCGH.obj)),
             chrominfo = human.chrom.info.Jul03,
             X = TRUE, Y = FALSE, rsp.uniq = unique(pheno),
             all = length(rsp.uniq) == 1 && is.null(resT),
             titles = if (all) "All Samples" else rsp.uniq,
             cutplot = 0, thres = .25, factor=2.5, ylm = c(-1, 1), 
	     p.thres = c(.01, .05, .1), numaut = 22, onepage = TRUE,
             colored = TRUE
             ){

#check if sd.samples are non-empty:
    if (!is.null(sd.samples(aCGH.obj)))
	{
		thres <- factor*(sd.samples(aCGH.obj)$madGenome)
	}
    col.scheme <- 
        if (colored)
            list(pal =
                 c("red", "blue", "green", "orange")[
                                                     1:length(p.thres)
                                                     ],
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
                 c("grey10", "grey40", "grey70", "grey90")[
                                                           1:length(p.thres)
                                                           ],
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
    
    rsp.uniq <- sort(rsp.uniq)
    
    ## creating response matrix colmatr
    
    colmatr <-
        if (length(rsp.uniq) > 1)
            t(
              sapply(rsp.uniq,
                     function(rsp.uniq.level)
                     ifelse(pheno == rsp.uniq.level, 1, 0)
                     )
              )
        else
            matrix(rep(1, length(pheno)),
                   ncol = length(pheno),
                   nrow = 1
                   )
    
    ## screening out clones that are gained or lost in < minChanged in
    ## classes under comparison indeces present: DON't do this anymore
    ## since rest is compute seprately : NOT ANYMORE

    #tmp <- apply(as.matrix(colmatr), 2, sum)
    #indecesnow <- which(tmp == 1)
    #data.thres <- threshold.func(data, thresAbs = thres)
    #prop.ch <- changeProp.func(dat = data.thres, colMatr = colmatr)
    #maxch <- changeProp.overall.func(dat = data.thres[ ,indecesnow ])
    #clones.index <- which(maxch >= minChanged)
     

    ##removing clones to skip from the dataset NOT ANYMORE

    #data <- data[clones.index,]
    #data.thres <- data.thres[clones.index,]
    #datainfo <- datainfo[clones.index,]


    nr <- nrow(colmatr)
    
    ## if 1 class only or missing resT, no significance analysis. 
    ##Otherwise: extra figure
    if (!is.null(resT))
        nr <- nr + 1
    
    tmp <- as.data.frame(matrix(0, ncol = 2, nrow = 1))
    colnames(tmp) <- c("gainP", "lossP")
    gainloss <-
        lapply(1:nrow(colmatr),
               function(j)
               gainLoss(dat = data,
                             cols = which(colmatr[j, ] == 1),
                             thres = thres)
               )
    dt <- data[ ,colmatr[ 1, ] == 1, drop = FALSE ]
    rsp <- rep(1, ncol(dt))
    if (nrow(colmatr) > 1)
        for (j in 2:nrow(colmatr))
        {
            
            dt <- cbind(dt, data[ ,colmatr[ j, ] == 1])
            rsp <- c(rsp, rep(j, sum(colmatr[ j, ] == 1)))
            
        }
    rsp <- rsp - 1
    

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

        ind <- which(gl$gainP >= cutplot)
        plot(kb.loc[ind], gl$gainP[ind],
             col = "green",
             type = "h", xlab = "chromosome number",
             ylab = "Fraction gained or lost", pch = 18, main = tl,
             ylim = ylm,
             xlim = c(0, max(cumsum(chrominfo$length), kb.loc[ind],
             rm.na = TRUE)), xaxt="n")
             
	axis(side=1, at=kb.loc[ind][1], label="", tick=FALSE)
        ind <- gl$lossP >= cutplot
        points(kb.loc[ind], -gl$lossP[ind],
               col = "red",
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
            if (is.even(numaut))
                mtext("X", side = 1, at = (chrom.mid[numaut + 1]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
            else
                mtext("X", side = 3, at = (chrom.mid[numaut + 1]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
        if (Y)
            if (is.even(numaut))
                mtext("Y", side = 3, at = (chrom.mid[numaut + 2]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
            else
                mtext("Y", side = 1, at = (chrom.mid[numaut + 2]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
        
    }
    if (!is.null(resT))
    {
        
        ## Process statistics
        ## for plotting test stats and p-values
        
        res <- resT[order(resT$index) ,]
        #res <- res[res$index %in% clones.index ,]
        maxT <- res$adjp
        teststat <- abs(res$teststat)
        st <-
            sapply(p.thres,
                   function(threshold) {
                       
                       if (any(maxT <= threshold))
                           min(teststat[maxT <= threshold])
                       else
                           NA
                       
                   }
                   )
        plot(kb.loc, teststat, col = col.scheme$kb.loc,
             ylim = c(0, max(teststat)), type = "h",
             xlab = "chromosome number", ylab = "clone statistic",
             pch = 18, main = paste(titles, collapse = " vs "),
             xlim =
             c(0, max(cumsum(chrominfo$length), kb.loc, rm.na = TRUE)), xaxt="n"
             )
	axis(side=1, at=kb.loc[ind][1], label="", tick=FALSE)	
	st.now <- rev(st)
	pal.now <- rev(col.scheme$pal)
        if (length(st.now) > 0)
            abline(h = st.now, col = pal.now, lty = 2)

    
    abline(v = cumsum(chrominfo$length), col = col.scheme$abline3)
    abline(v = chrom.centr, lty = 2, col = col.scheme$abline4)

    for (i in seq(1, numaut, b = 2))
        mtext(paste("", i), side = 1, at = chrom.mid[i], line = .3,
              col = col.scheme$mtext, cex.main = .5)
    for (i in seq(2, numaut, b = 2))
        mtext(paste("", i), side = 3, at = chrom.mid[i], line = .3,
              col = col.scheme$mtext, cex.main = .5)
    if (X)
            if (is.even(numaut))
                mtext("X", side = 1, at = (chrom.mid[numaut + 1]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
            else
                mtext("X", side = 3, at = (chrom.mid[numaut + 1]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
        if (Y)
            if (is.even(numaut))
                mtext("Y", side = 3, at = (chrom.mid[numaut + 2]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
            else
                mtext("Y", side = 1, at = (chrom.mid[numaut + 2]),
                      line = .3, col = col.scheme$mtext, cex.main = .5)
        
	}
}

summarize.clones <-
    function(aCGH.obj, resT = NULL, pheno = rep(1, ncol(aCGH.obj)),
             rsp.uniq = unique(pheno), thres = .25,factor=2.5, 
             all = length(rsp.uniq) == 1 && is.null(resT), titles = if (all) "all" else rsp.uniq)
{

	#check if sd.samples are non-empty:
    if (!is.null(sd.samples(aCGH.obj)))
	{
		thres <- factor*(sd.samples(aCGH.obj)$madGenome)
	}
    
    data <- log2.ratios(aCGH.obj)
    datainfo <- clones.info(aCGH.obj)
    
    rsp.uniq <- sort(rsp.uniq)

    colmatr <-
        if (length(rsp.uniq) > 1)
            t(
              sapply(rsp.uniq,
                     function(rsp.uniq.level)
                     ifelse(pheno == rsp.uniq.level, 1, 0)
                     )
              )
        else
            matrix(rep(1, length(pheno)),
                   ncol = length(pheno),
                   nrow = 1
                   )

	data.thres <- threshold.func(data, thresAbs = thres)

####NOT ANYMORE
    #tmp <- apply(as.matrix(colmatr), 2, sum)
    #indecesnow <- which(tmp == 1)
    #prop.ch <- changeProp.func(dat = data.thres, colMatr = colmatr)
    #maxch <- changeProp.overall.func(dat = data.thres[ ,indecesnow ])
    #clones.index <- which(maxch >= minChanged)
    #data <- data[clones.index,]
    #data.thres <- data.thres[clones.index,]
    #dataSign <- dataSign[clones.index,]
    #datainfo <- datainfo[clones.index,]

    bac.summary <- table.bac.func(dat = data.thres, colMatr = colmatr)
    if (!is.null(resT))
    {
        
        ## Process statistics
        ## for plotting test stats and p-values
        
        res <- resT[order(resT$index) ,]
        #res <- res[res$index %in% clones.index ,]
        bac.summary <- cbind(bac.summary, res$teststat, res$rawp, res$adjp)

    }
    bac.summary <- as.data.frame(bac.summary)
    nms <-
        c("NumPresent", "NumGain", "NumLost", "PropPresent",
          "PropGain", "PropLost")
    cnames <- colnames(bac.summary)
    cnames[1:6] <- paste(nms, "All", sep = ".")
    if (nrow(colmatr) > 1)
        for (m in 1:length(rsp.uniq))
            cnames[ (6 * m + 1):(6 * (m + 1)) ] <-
                paste(nms, titles[m], sep = ".")
    if (!is.null(resT))
        cnames[(ncol(bac.summary) - 2):ncol(bac.summary)] <-
            c("stat", "rawp", "adjp")
    colnames(bac.summary) <- cnames
    bac.summary <- cbind(datainfo, bac.summary)
    invisible(bac.summary)
    
}

