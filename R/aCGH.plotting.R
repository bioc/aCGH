plotGenome <-
    function(aCGH.obj, samples = 1:num.samples(aCGH.obj), naut = 22,
             Y = TRUE, X = TRUE, data = log2.ratios(aCGH.obj),
             chrominfo = human.chrom.info.Jul03, yScale = c(-2, 2),
             samplenames = sample.names(aCGH.obj), ylb = "Log2Ratio")
{

    datainfo <- clones.info(aCGH.obj)
    ##total number of chromosomes to plot:
    
    nchr <- naut
    if (X)
	nchr <- nchr + 1
    if (Y)
        nchr <- nchr + 1
    
    nsamples <- length(samplenames)
    
    ##reordering according to genomic position
    ord <- order(datainfo$Chrom, datainfo$kb)
    chrom <- datainfo$Chrom[ord]
    kb <- datainfo$kb[ord]
    data <- data[ord,]
    
    ##screening out unampped clones
    ind.unmap <- which(is.na(chrom) | is.na(kb) | (chrom > (naut+2)))
    if (length(ind.unmap) > 0)
    {
        
        chrom <- chrom[-ind.unmap]
    	kb <- kb[-ind.unmap]
    	data <- data[-ind.unmap ,]
        
    }
    
    ##removing chromosome not to plot:
    data <- data[chrom <= nchr ,]
    kb <- kb[chrom <= nchr]
    chrom <- chrom[chrom <= nchr]

    chrominfo <- chrominfo[1:nchr,]
    chrom.start <- c(0, cumsum(chrominfo$length))[1:nchr]
    chrom.centr <- chrom.start + chrominfo$centr
    chrom.mid <- chrom.start + chrominfo$length / 2
    chrom.rat <- chrominfo$length / max(chrominfo$length)

    par(cex = .6, pch = 18, lab = c(1, 6, 7), cex.axis = 1.5,
        xaxs = "i")
    for (k in 1:length(samples))
    {

        vec <- data[ ,samples[k] ]
        name <- samplenames[samples[k]]
        clone.genomepos <- rep(0, length(kb))
        for (i in 1:nrow(chrominfo))
            clone.genomepos[chrom == i] <-
                kb[chrom == i] + chrom.start[i]

        ##Now, determine vertical scale for each chromosome:

        y.ranges <-
            sapply(1:nrow(chrominfo),
                   function(i)
                   range(vec[chrom == i], yScale, na.rm = TRUE)
                   )
        ylim <- range(y.ranges)
###        y.min <- y.ranges[1 ,]
###        y.max <- y.ranges[2 ,]
        
#########################
        
        plot(clone.genomepos / 1000, vec, ylim = ylim, xlab = "",
             ylab = "",
             xlim =
             c(min(clone.genomepos[clone.genomepos > 0], na.rm = TRUE) /
               1000,
               clone.genomepos[sum(clone.genomepos > 0)] / 1000),
             col = "black")
        
        
        ##title(main=paste(name, " ", sample[k], " - Whole Genome"),
        ##ylab=ylb, xlab="Chromosome", cex.lab=1.5,cex.main=2)
        title(main = paste(samples[k], " ", name), ylab = ylb,
              xlab = "", cex.lab = 1.5, cex.main = 2)
        
        for (i in seq(1, naut, b = 2))
            mtext(paste("", i), side = 1, at = chrom.mid[i] / 1000,
                  line = .3, col = "red")
        for (i in seq(2, naut, b = 2))
            mtext(paste("", i), side = 3, at = chrom.mid[i] / 1000,
                  line = .3, col = "red")
        
        if (X)
            mtext("X", side = 1, at = chrom.mid[naut + 1] / 1000,
                  line = .3, col = "red")
        if (Y)
            mtext("Y", side = 3, at = chrom.mid[naut + 2] / 1000,
                  line = .3, col = "red")
        
        abline(v = c(chrom.start / 1000,
               (chrom.start[nrow(chrominfo)] +
                chrominfo$length[nrow(chrominfo)]) / 1000), lty = 1)
        abline(h = seq(-1 , 1, b = .5), lty = 3)
        abline(v = (chrominfo$centromere + chrom.start) / 1000,
               lty = 3, col = "red")
        
    }
    invisible(list(x = clone.genomepos / 1000,
                   ylim = range(y.ranges)))
    
}

###############################

plotSummaryProfile <-
    function(aCGH.obj,
             response = as.factor(rep("All", ncol(aCGH.obj))),
             titles = unique(response[!is.na(response)]), X = TRUE,
             Y = FALSE, maxChrom = 23,
             chrominfo = human.chrom.info.Jul03,
             num.plots.per.page = length(titles))
{

    if (is.null(genomic.events(aCGH.obj)))
        stop("compute the genomic events first using\
find.genomic.events")
    ind.samp <- which(!is.na(response))
    resp.na <- response[ind.samp]
    response.uniq <- sort(unique(resp.na))
    ge <- genomic.events(aCGH.obj)

###    length.num.func <-
###        function(x, num)
###            sapply(num, function(nn) sum(x == nn & !is.na(x)))

    df.not.na <-
        data.frame(response = response,
                   numtrans =
                   apply(ge$num.transitions, 2, sum, na.rm = TRUE),
                   numtrans.binary =
                   apply(ge$num.transitions.binary, 2, sum, na.rm = TRUE),
                   numaber =
                   apply(ge$num.aberrations, 2, sum, na.rm = TRUE),
                   numaber.binary =
                   apply(ge$num.aberrations.binary, 2, sum, na.rm = TRUE),
                   numamplif =
                   apply(ge$num.amplifications[ 1:maxChrom, ], 2, sum,
                         na.rm = TRUE),
                   numamplif.binary =
                   apply(ge$num.amplifications.binary[ 1:maxChrom, ],
                         2, sum, na.rm = TRUE),
                   numoutlier =
                   apply(ge$num.outliers, 2, sum, na.rm = TRUE),
                   num.outliers.binary =
                   apply(ge$num.outliers.binary, 2, sum, na.rm = TRUE),
                   numchromgain =
                   apply(ge$whole.chrom.gain.loss[ 1:maxChrom, ], 2,
                         length.num.func, 1),
###                   apply(ge$whole.chrom.gain[ 1:maxChrom, ], 2, sum,
###                         na.rm = TRUE),
                   numchromloss =
                   apply(ge$whole.chrom.gain.loss[ 1:maxChrom, ], 2,
                         length.num.func, -1),
###                   apply(ge$whole.chrom.loss[ 1:maxChrom, ], 2, sum,
###                         na.rm = TRUE)
                   sizeamplicon =
                   apply(ge$size.amplicons[ 1:maxChrom, ], 2, sum,
                         na.rm = TRUE),
                   numamplicon =
                   apply(ge$num.amplicons[ 1:maxChrom, ], 2, sum,
                         na.rm = TRUE)
                   )[ which(!is.na(response)), ]
    attach(df.not.na)
    numchromchange <- numchromgain + numchromloss
    
    boxplot.this <-
        function(events, title, sig = 6)
        {

            p.value <-
                if (length(response.uniq) > 1)	
                    signif(kruskal.test(events ~ resp.na)$p.value,
                           sig)
                else
                    ""
            boxplot(events ~ resp.na, notch = TRUE, names = titles,
                    varwidth = TRUE, main = paste(title, p.value))
            
        }

#############################################
    ##Plot1:
    par(mfrow = c(2, 2))

    boxplot.this(numtrans, "Number of Transitions")
    boxplot.this(numtrans.binary,
                 "Number of Chrom containing Transitions")
    boxplot.this(numaber, "Number of Aberrations")
    boxplot.this(numchromchange, "Number of Whole Chrom Changes")

#############################################
    ##Plot2:
    
    boxplot.this(numamplif, "Number of Amplifications")
    boxplot.this(numamplif.binary,
                 "Number of Chrom containing Amplifications")
    boxplot.this(numamplicon, "Number of Amplicons")
    boxplot.this(sizeamplicon, "Amount of Genome Amplified")

#############################################

    plot.freq.this <-
        function(matr, i, ylb)
        {
            
            par(mfrow = c(num.plots.per.page, 1))

            out <-
                sapply(1:length(response.uniq),
                       function(j)
                       apply(matr[
                                  ,which(resp.na == response.uniq[j])
                                  ],
                             1,
                             prop.num.func,
                             i
                             )
                       )
            mx <- max(c(out), na.rm = TRUE)
            if (length(titles == 1))
                out <- cbind(out, out)
            
            plotGenome(aCGH.obj, samples = 1:length(titles),
                       yScale = c(0, mx), data = out, naut = 22,
                       X = X, Y = Y, ylb = ylb,
                       chrominfo = chrominfo[ 1:maxChrom, ],
                       samplenames = titles
                       )
            
        }
    
    ##Plot3: trans start
    plot.freq.this(ge$transitions$trans.matrix, 1,
                  "Proportion of Transition Starts")

    ##Plot4: trans end
    plot.freq.this(ge$transitions$trans.matrix, 2,
                  "Proportion of Transition Ends")

    ##Plot5: amplification
    plot.freq.this(ge$amplifications$amplif, 1,
                  "Proportion of Amplifications")
    mtext("Amplifications", side = 3, outer = TRUE)

    ##Plot6: aberration
    plot.freq.this(ge$aberrations$aber, 1,
                   "Proportion of Aberrations")
    mtext("Aberrations", side = 3, outer = TRUE)

    ##Plot7: whole chromosomal gain/loss:

    par(mfrow = c(num.plots.per.page, 2), lab = c(5,6,7))
    
    matr <- ge$whole.chrom.gain.loss[ 1:22, ]
    out.gain <- matrix(NA, nrow = nrow(matr), ncol = length(titles))
    out.loss <- matrix(NA, nrow = nrow(matr), ncol = length(titles))
    for (j in 1:length(response.uniq))
    {
        
        ind <- which(response == response.uniq[j])
        out.gain[ ,j ] <-
            apply(matr[ ,ind ], 1, length.num.func, 1) / ncol(matr)
        out.loss[ ,j ] <-
            apply(matr[ ,ind ], 1, length.num.func, -1) / ncol(matr)
        
    }
    mx.gain <- max(c(out.gain), na.rm = TRUE)
    mx.loss <- max(c(out.loss), na.rm = TRUE)
    mx <- max(mx.gain, mx.loss)
    for (j in 1:length(titles))
    {
        
        plot(1:22, out.gain[,j], pch = 20,
             main = as.character(titles[j]),
             xlab = "chromosome",
             ylab = "Proportion of whole chromosomes gains",
             ylim = c(0, mx), xlim = c(0,23))
        plot(1:22, out.loss[,j], pch = 20,
             main = as.character(titles[j]),
             xlab = "chromosome",
             ylab = "Proportion of whole chromosomes losses",
             ylim = c(0, mx), xlim = c(0,23))
        
    }
    detach(df.not.na)

}	

plotHmmStates <-
    function(aCGH.obj, sample.ind, chr = 1:num.chromosomes(aCGH.obj),
             statesres = hmm.merged(aCGH.obj), maxChrom = 23,
             chrominfo = human.chrom.info.Jul03, yScale = c(-2, 2),
             samplenames = sample.names(aCGH.obj)
             )
{

    if (is.null(statesres))
        stop("merge the hmm states first using merge.hmm.states\
function")
    if (length(sample.ind) > 1)
        stop("plotHmmStates currently prints only 1 sample at a\
time\n")
    if (is.null(genomic.events(aCGH.obj)))
        stop("compute the genomic.events of aCGH.obj first using\
find.genomic.events function")
    
###    hmm.res.merge <- merge.func(statesres = statesres, minDiff = .25)
###    states.bic <- hmm.res.merge$states.hmm
    ge <- genomic.events(aCGH.obj)
    aber <- ge$aberrations$aber
    amplif <- ge$amplifications$amplif
    trans <- ge$transitions$trans.matr
    outliers <- ge$outliers$outlier
    pred <- ge$outliers$pred.out

    chrom.rat <- chrominfo$length / max(chrominfo$length)
    chrom.start <- c(0, cumsum(chrominfo$length))[1:maxChrom]
    
    ##chrom.mid contains middle positions of the chromosomes relative to
    ##the whole genome (useful for plotting the whole genome)
    chrom.mid <- chrom.start + chrominfo$length[1:maxChrom] / 2
    chrom <- statesres[ ,1 ]
    par(lab = c(15, 6, 7), pch = 18, cex = 1, lwd = 1,
        mfrow = c(2, 1))
    
    sq.state <- seq(3, ncol(statesres), b = 6)
    sq.obs <- seq(8, ncol(statesres), b = 6)
    
    for (j in 1:length(chr))
    {

        ind.nonna <-
            which(!is.na(statesres[chrom == chr[j],
                                   sq.obs[sample.ind]]))

        kb <- statesres[chrom == chr[j], 2][ind.nonna] / 1000
        obs <- statesres[chrom == chr[j],
                         sq.obs[sample.ind]][ind.nonna]
        states <-
            statesres[chrom == chr[j],
                      sq.state[sample.ind]][ind.nonna]
        nstates <- length(unique(states)) 

        abernow <- aber[chrom == chr[j], sample.ind][ind.nonna]
        outliersnow <-
            outliers[chrom == chr[j], sample.ind][ind.nonna]
        amplifnow <- amplif[chrom == chr[j], sample.ind][ind.nonna]
        transnow <- trans[chrom == chr[j], sample.ind][ind.nonna]

        ## predicted values when non-aberration of outlier: otherwise
        ## observed
        
        prednow <- obs
        predicted <- pred[chrom == chr[j], sample.ind][ind.nonna]
        prednow[outliersnow == 0 & abernow == 0] <-
            predicted[outliersnow == 0 & abernow == 0]

        y.min <- min(yScale[1], min(obs))
        y.max <- max(yScale[2], max(obs))

        ##observed

        plot(kb, obs, xlab = "", ylab = "", ylim = c(y.min, y.max),
             type = "l", col = "blue",
             xlim = c(0, chrominfo$length[chr[j]] / 1000)
             )
        points(kb, obs, col = "black")
        title(main = paste("Sample", sample.ind,
              samplenames[sample.ind], "- Chr", chr[j],
              "Number of states", nstates),
              xlab = "kb (in 1000's)", ylab = "data (observed)"
              )
        
        abline(h = seq(-2, 2, b = .5), lty = 3)
        abline(v = chrominfo$centromere[chr[j]] / 1000, lty = 2,
               col = "red", lwd = 3)

        if (nstates > 1)
        {
            
            abline(v = kb[transnow == 1], col = "blue", lwd = 2)
            abline(v = kb[transnow == 2], col = "green", lty = 2,
                   lwd = .5)
            
        }

        ##amplif = red
        ##aber = orange
        ##outliers = yellow

        if (length(outliersnow[outliersnow == 1]) > 0)
            points(kb[outliersnow == 1], obs[outliersnow == 1],
                   col = "yellow")
        if (length(abernow[abernow == 1]) > 0)
            points(kb[abernow == 1], obs[abernow == 1],
                   col = "orange")
        if (length(amplifnow[amplifnow == 1]) > 0)
            points(kb[amplifnow == 1], obs[amplifnow == 1],
                   col="red")

        ##predicted states:
        
        plot(kb, prednow, xlab = "", ylab = "",
             ylim = c(y.min, y.max), type = "l", col = "blue",
             xlim =c(0, chrominfo$length[chr[j]] / 1000))
        
        points(kb, prednow, col = "black")
        title(xlab = "kb (in 1000's)", ylab = "data (smoothed)")
        abline(h = seq(-2, 2, b = .5), lty = 3)
        abline(v = chrominfo$centromere[chr[j]] / 1000, lty = 2,
               col = "red", lwd = 3)

        ##start (dotted blue) and end of states (green)
        if (nstates > 1)
        {
            
            abline(v = kb[transnow == 1], col = "blue", lwd = 2)
            abline(v = kb[transnow == 2], col = "green", lty = 2,
                   lwd = .5)
            
        }

        ##amplif = red
        ##aber = orange
        ##outliers = yellow

        if (length(outliersnow[outliersnow == 1]) > 0)
            points(kb[outliersnow == 1], obs[outliersnow == 1],
                   col = "yellow")
        if (length(abernow[abernow == 1]) > 0)
            points(kb[abernow == 1], obs[abernow == 1],
                   col = "orange")
        if (length(amplifnow[amplifnow == 1]) > 0)
            points(kb[amplifnow == 1], obs[amplifnow == 1],
                   col = "red")

    } 

}

plotHmmStatesNew <-
    function(aCGH.obj, sample.ind, chr = 1:num.chromosomes(aCGH.obj),
             statesres = hmm(aCGH.obj)$states.hmm[[1]],
             chrominfo = human.chrom.info.Jul03, yScale = c(-2, 2),
             samplenames = sample.names(aCGH.obj)
             )
{

    if (length(sample.ind) > 1)
        stop("plotHmmStatesNew currently prints only 1 sample at a\
time\n")
    if (is.null(genomic.events(aCGH.obj)))
        stop("compute the genomic.events of aCGH.obj first using\
find.genomic.events function")
    
    ge <- genomic.events(aCGH.obj)
    aber <- ge$aberrations$aber
    amplif <- ge$amplifications$amplif
    trans <- ge$transitions$trans.matr
    outliers <- ge$outliers$outlier
    pred <- ge$outliers$pred.out

    chrom.rat <- chrominfo$length / max(chrominfo$length)
    chrom.start <- c(0, cumsum(chrominfo$length))[chr]
    
    ##chrom.mid contains middle positions of the chromosomes relative to
    ##the whole genome (useful for plotting the whole genome)
    chrom.mid <- chrom.start + chrominfo$length[chr] / 2
    chrom <- statesres[ ,1 ]
    par(lab = c(15, 6, 7), pch = 18, cex = 1, lwd = 1,
        mfrow = c(2, 1))
    
    sq.state <- seq(3, ncol(statesres), b = 6)
    sq.obs <- seq(8, ncol(statesres), b = 6)
    
    for (j in 1:length(chr))
    {

        ind.nonna <-
            which(!is.na(statesres[chrom == chr[j],
                                   sq.obs[sample.ind]]))

        kb <- statesres[chrom == chr[j], 2][ind.nonna] / 1000
        obs <- statesres[chrom == chr[j],
                         sq.obs[sample.ind]][ind.nonna]
        states <-
            statesres[chrom == chr[j],
                      sq.state[sample.ind]][ind.nonna]
        nstates <- length(unique(states)) 

        abernow <- aber[chrom == chr[j], sample.ind][ind.nonna]
        outliersnow <-
            outliers[chrom == chr[j], sample.ind][ind.nonna]
        amplifnow <- amplif[chrom == chr[j], sample.ind][ind.nonna]
        transnow <- trans[chrom == chr[j], sample.ind][ind.nonna]

        ## predicted values when non-aberration of outlier: otherwise
        ## observed
        
        prednow <- obs
        predicted <- pred[chrom == chr[j], sample.ind][ind.nonna]
        prednow[outliersnow == 0 & abernow == 0] <-
            predicted[outliersnow == 0 & abernow == 0]

        y.min <- min(yScale[1], min(obs))
        y.max <- max(yScale[2], max(obs))

        ##observed

        plot(kb, obs, xlab = "", ylab = "", ylim = c(y.min, y.max),
             type = "l", col = "blue",
             xlim = c(0, chrominfo$length[chr[j]] / 1000)
             )
        points(kb, obs, col = "black")
        title(main = paste("Sample", sample.ind,
              samplenames[sample.ind], "- Chr", chr[j],
              "Number of states", nstates),
              xlab = "kb (in 1000's)", ylab = "data (observed)"
              )
        
        abline(h = seq(-2, 2, b = .5), lty = 3)
        abline(v = chrominfo$centromere[chr[j]] / 1000, lty = 2,
               col = "red", lwd = 3)

        if (nstates > 1)
        {
            
            abline(v = kb[transnow == 1], col = "blue", lwd = 2)
            abline(v = kb[transnow == 2], col = "green", lty = 2,
                   lwd = .5)
            
        }

        ##amplif = red
        ##aber = orange
        ##outliers = yellow

        if (length(outliersnow[outliersnow == 1]) > 0)
            points(kb[outliersnow == 1], obs[outliersnow == 1],
                   col = "yellow")
        if (length(abernow[abernow == 1]) > 0)
            points(kb[abernow == 1], obs[abernow == 1],
                   col = "orange")
        if (length(amplifnow[amplifnow == 1]) > 0)
            points(kb[amplifnow == 1], obs[amplifnow == 1],
                   col="red")

        ##predicted states:
        
        plot(kb, prednow, xlab = "", ylab = "",
             ylim = c(y.min, y.max), type = "l", col = "blue",
             xlim =c(0, chrominfo$length[chr[j]] / 1000))
        
        points(kb, prednow, col = "black")
        title(xlab = "kb (in 1000's)", ylab = "data (smoothed)")
        abline(h = seq(-2, 2, b = .5), lty = 3)
        abline(v = chrominfo$centromere[chr[j]] / 1000, lty = 2,
               col = "red", lwd = 3)

        ##start (dotted blue) and end of states (green)
        if (nstates > 1)
        {
            
            abline(v = kb[transnow == 1], col = "blue", lwd = 2)
            abline(v = kb[transnow == 2], col = "green", lty = 2,
                   lwd = .5)
            
        }

        ##amplif = red
        ##aber = orange
        ##outliers = yellow

        if (length(outliersnow[outliersnow == 1]) > 0)
            points(kb[outliersnow == 1], obs[outliersnow == 1],
                   col = "yellow")
        if (length(abernow[abernow == 1]) > 0)
            points(kb[abernow == 1], obs[abernow == 1],
                   col = "orange")
        if (length(amplifnow[amplifnow == 1]) > 0)
            points(kb[amplifnow == 1], obs[amplifnow == 1],
                   col = "red")

    } 

}

##plot.all.arrays.plot <-
##    function(aCGH.obj, chrominfo = human.chrom.info.Jul03,
##             sample.names = sample.names(aCGH.obj))
##{

##    if (is.null(hmm(aCGH.obj)))
##        stop("compute the hmm states first using find.hmm.states\
##function")
##    if (is.null(sd.samples(aCGH.obj)))
##        stop("compute first the standard deviations for samples using\
##computeSD.Samples function")

##    madGenome <- sd.samples(aCGH.obj)$madGenome
##    ord <- order(madGenome)
##    nm <-
##        apply(data.frame(sample.names, round(madGenome, 2),
##                         round((3 * madGenome), 2)),
##              1,
##              paste,
##              collapse = ":"
##              )
##    op <- par(mfrow = c(2, 1))
##    sapply(ord,
##           function(i) {
               
##               cat(i, "\n")
##               plotGenome(aCGH.obj, sample = i)
               
##           }
##           )
##    par(op)
    
##}

##plot.all.arrays <-
##    function(aCGH.obj, plot.file.name = "samples.noise.ps")
##{

##    postscript(plot.file.name, paper="letter")
##    plot.all.arrays.plot(aCGH.obj)
##    dev.off()
    
##}

##produce.pairs.plot <-
##    function(aCGH.obj)
##{

##    attach(aCGH.obj)
##    nm <- colnames(qual.rep)
##    for (i in 1:length(nm))
##    {
##        cat(i, "\n")
##        matr <-
##            t(as.matrix(log2.ratios)[
##                                     which(clones.info(aCGH.obj)$Clone ==
##                                           nm[i]),
##                                     ]
##              )
##        ylm <- c(min(c(matr), na.rm = TRUE), max(c(matr), na.rm = TRUE))
##        pairs(matr, ylim = c(ylm[1], ylm[2]),
##              xlim = c(ylm[1], ylm[2]), pch = 20)
##        title(paste(nm[i], " ;median maxdiff is ",
##                    qual.rep[2,i], " ;spearman corr is ",
##                    qual.rep[3,i])
##              )
##    }
##    detach(aCGH.obj)
    
##}

plotValGenome <-
    function(aCGH.obj, phen = rep(1, ncol(aCGH.obj)),
             data = log2.ratios(aCGH.obj),
             datainfo = clones.info(aCGH.obj),
             chrominfo = human.chrom.info.Jul03,
             cutoff = 1, ncolors = 50, byclass = TRUE, showaber = FALSE,
             amplif = 1, homdel = -1, vecchrom = 1:23,
             samplenames = sample.names(aCGH.obj),
             title = "Image Plots")
{
    
    resp0 <- phen
    resp <- resp0
    if (!byclass)
	resp <- rep(1, length(resp0))

    tbl.resp <- table(resp)
    ##label.col <- c("red", "green", "blue", "skyblue", "orange", "pink", "gray20")
    label.col <- rainbow(6)
    
    par(bg = "grey20")
    
    kb <- datainfo$kb
    data <- as.matrix(data)
    dt.cp <- data
    dt <- apply(data, 2,floor.func, cutoff)    
###    chromb <- rep(0,nrow(chrominfo))
    ##centrloc <- rep(0,nrow(chrominfo))
###    for (i in 1:nrow(chrominfo))
###    {
        
###	for (j in 1:i)
###            ##chromb[i] <- chromb[i]+length(datainfo$Chrom[datainfo$Chrom==vecchrom[j]])+.5
###            chromb[i] <- chromb[i]+sum(datainfo$Chrom == vecchrom[j])
###        ##centrloc[i] <- chromb[i]-length(datainfo$Chrom[datainfo$Chrom==vecchrom[i] & datainfo$kb >= chrominfo$centr[vcchrom[i]]])
	
###    }
    chromb <- c(0, cumsum(table(datainfo$Chrom)))
    ##chromb <- c(.5, chromb)
    ##chromb <- c(0, chromb)

    dt <- dt[ ,order(resp) ]
    dt.cp <- dt.cp[ ,order(resp) ]
    resp0 <- resp0[order(resp)]
    samplenames <- samplenames[order(resp)]
    resp <- resp[order(resp)]
    start <- 1
    
    ##mapping order
    ord <- rep(0, length(resp))
    for (i in 1:(length(tbl.resp)))
    {
	
	ind <- which(resp == i)
	cr <- as.dist(1 - corna(data[ ,ind ]))
	ord[start:sum(tbl.resp[1:i])] <-
            hclust(cr, method = "ave")$ord + start - 1
	start <- sum(tbl.resp[1:i]) + 1
	
    }
    dt <- dt[ ,ord ]
    dt.cp <- dt.cp[ ,ord ]
    resp <- resp[ord]
    resp0 <- resp0[ord]
    samplenames <- samplenames[ord]

    image(x = 1:length(kb), y = 1:length(resp), z = dt,
          col = maPalette(low = "red", high = "green", mid = "white",
          k = ncolors), axes = FALSE, xlab = "", ylab = "",
          zlim = c(-cutoff, cutoff))
    
    ##abline(h=seq(.5, 81.5, b=1), col="gray20", lwd=.2)
    if (showaber)
    {
        
        ##for (i in 1:nrow(dt))
        ##{
        for (j in 1:ncol(dt))
        {
            
            tmp <- dt.cp[,j]
            i <- which(tmp >= amplif & !is.na(tmp))
            if (length(i) > 0)
                ##if ((!is.na(dt.cp)) && (dt.cp[i,j] >= amplif))
                points(i, rep(j, length(i)), col = "yellow", pch = 20,
                       cex = .7)
            i <- which(tmp <= homdel & !is.na(tmp))
            if (length(i) > 0)
                ##if ((!is.na(dt.cp)) && (dt.cp[i,j] >= amplif))
                points(i, rep(j, length(i)), col = "skyblue",
                       pch = 20, cex = .7)
            
        }
        ##}
    }
    for (j in 1:ncol(dt))
    {

        col <- label.col[resp0[j] + 1]
	mtext(resp0[j], side = 2, at = j, line=.3, col = col,
              cex = .5, las = 2)
	mtext(paste((samplenames)[j], ""), side = 4, at = j,
              line = .3, col = col, cex = .3, las = 2)
	
    }
    ##title(main="Whole genome", xlab = "clone", ylab = "sample", col.lab="white", col.main="white")
    title(xlab = "clone", ylab = "sample", col.lab = "white",
          col.main = "white", main = title)
    ##abline(v=centrloc, col="white", lty=2, lwd=.5)
    abline(v = chromb, col = "black", lty = 1, lwd = .5)
    loc <- chromb[-1] - diff(chromb) / 2
    for (i in seq(2, nrow(chrominfo), b = 2))
        mtext(paste("", vecchrom[i]), side = 3, at = loc[i],
              line = .3,col = "white", cex.main = .5)
    for (i in seq(1, nrow(chrominfo), b = 2))
        mtext(paste("", vecchrom[i]), side = 1, at = loc[i],
              line = .3,col = "white", cex.main = .5)
    ##mtext("X", side = 1, at = (loc[nrow(chrominfo)]), line=.3,col="white", cex.main=.5)

}

plotValChrom <-
    function(aCGH.obj, phen = rep(1, ncol(aCGH.obj)),
             data = log2.ratios(aCGH.obj),
             datainfo = clones.info(aCGH.obj),
             chrominfo = human.chrom.info.Jul03, chrom = 1:23,
             cutoff = 1, ncolors = 50, amplif = 1, homdel = -1,
             byclass = TRUE, samplenames = sample.names(aCGH.obj),
             clonenames = datainfo$Clone, title = "Image Plot")
{
    
    ##label.col <- c("red", "green", "blue", "yellow", "skyblue", "orange", "pink", "gray20")
    label.col <- rainbow(6)
    par(bg = "grey20")    
    samplenames.cp <- samplenames
    for (chr in chrom)
    {
        
        resp0 <- phen
        resp <- resp0
        samplenames <- samplenames.cp
        if (!byclass)
            resp <- rep(1, length(resp0))
        tbl.resp <- table(resp)
        kb <- datainfo$kb[datainfo$Chrom == chr]
        dt <- as.matrix(data[ datainfo$Chrom == chr, ])
        clonenms <- clonenames[datainfo$Chrom == chr]
        
        dt.cp <- dt
        dt <- apply(dt.cp, 2, floor.func, cutoff)       
        if (chrominfo$centr[chr] >0)
###            centr <- length(kb[kb<=chrominfo$centr[chr]])
            centr <- sum(kb <= chrominfo$centr[chr])
        dt <- dt[,order(resp)]
        dt.cp <- dt.cp[ ,order(resp) ]
        resp0 <- resp0[order(resp)]
        samplenames <- samplenames[order(resp)]
        resp <- resp[order(resp)]
        start <- 1

        ##mapping order
        ord <- rep(0, length(resp))
        for (i in 1:length(tbl.resp))
        {
            
            ind <- which(resp == i)
            cr <- as.dist(1 - corna(dt.cp[ ,ind ]))
            ord[start:sum(tbl.resp[1:i])] <-
                hclust(cr, method = "ave")$ord + start - 1
            start <- sum(tbl.resp[1:i])+1
            
        }

        dt <- dt[ ,ord ]
        dt.cp <- dt.cp[ ,ord ]
        resp <- resp[ord]
        resp0 <- resp0[ord]
        samplenames <- samplenames[ord]

        image(x = 1:length(kb), y = 1:length(resp), z = dt,
              col = maPalette(low = "red", high = "green",
              mid = "white", k = ncolors), axes = FALSE, xlab = "",
              ylab = "", zlim = c(-cutoff, cutoff))

        if (chrominfo$centr[chr] > 0)
            abline(v = centr, col = "black")
        for (i in 1:nrow(dt))
        {

            ##if ((i %% 2) == 0)
            ##{
            ##	
            ##	mtext(paste(clonenms[i], ""), side = 1, at = i, line=.3, col="white", cex=.5, las=2)
            ##}
            ##else
            ##{	
            ##	mtext(paste(clonenms[i], ""), side = 3, at = i, line=.3, col="white", cex=.5, las=2)
            ##}
            mtext(paste(clonenms[i], ""), side = 1, at = i, line = .3,
                  col = "white", cex = .25, las = 2)
            for (j in 1:ncol(dt.cp))
            {
                
		if (i == 1)
		{

                    col <- label.col[resp0[j] + 1]
                    mtext(resp0[j], side = 2, at = j, line = .3,
                          col = col, cex = .5, las = 2)
                    mtext(paste(samplenames[j], ""), side = 4, at = j,
                          line = .3, col = col, las = 2, cex = .5)
                    
		}
		if (!is.na(dt.cp[i, j]) && dt.cp[i, j] >= amplif)
                    points(i, j, col = "yellow", pch = 20, cex = .7)
		if (!is.na(dt.cp[i, j]) && dt.cp[i, j] <= homdel)
                    points(i, j, col = "skyblue", pch = 20, cex = .7)
		
            }
            
        }
        title(main = paste(title, " Chromosome ", chr),
              col.main = "white")

    }

}

plotChrom <-
    function(aCGH.obj, sample = 1:ncol(aCGH.obj), chr = 1,
             yScale = c(-1, 1), data = log2.ratios(aCGH.obj),
             datainfo = clones.info(aCGH.obj),
             chrominfo = human.chrom.info.Jul03,
             samplenames = sample.names(aCGH.obj))
{
    
    nsamples <- length(sample)
    ord <- order(datainfo$Chrom, datainfo$kb)
    chrom <- datainfo$Chrom[ord]
    kb <- datainfo$kb[ord]
    data <- data[ ord, ]
    
    par(mfrow = c(nsamples, 1))
    par(cex = .6, pch = 18, lab = c(1,6,7), cex.axis = 1.5)
    kb <- kb[chrom == chr]
    centr.loc <- chrominfo$centromere[chr]
    for (k in 1:nsamples)
    {
        
        vec <- data[chrom == chr, sample[k]]
        name <- samplenames[sample[k]]
        
        ##Now, determine vertical scale for each chromosome:
        
        y.min <- min(yScale[1], minna(vec))
        y.max <- max(yScale[2], maxna(vec))

        plot(kb / 1000, vec, ylim = c(y.min, y.max), xlab = "",
             ylab = "",
             xlim = c(min(kb[kb > 0], na.rm = TRUE), kb[sum(kb > 0)]) /
             1000,
             col = "black", cex = 1.5)
        lines(kb / 1000, vec, col = "blue", lwd = .5)
        abline(v = centr.loc / 1000, col = "red", lty = 2)
        abline(h = 0, col = "black", lty = 2)
        abline(h = seq(-.6, .6, b = .2), lty = 3)
        title(main = paste(name, " Chr ", chr), ylab = "Log2Ratio",
              xlab = "Chromosome", cex.lab = 1.5, cex.main = 2)
        
    }

}

plotGene <-
    function(aCGH.obj, phen = rep(1, ncol(aCGH.obj)),
             data = log2.ratios(aCGH.obj), cutoff = 1, ncolors = 50,
             byclass = TRUE, method = "ave", showaber = FALSE, amplif = 1,
             homdel = -1, samplenames = sample.names(aCGH.obj),
             title = "Image Plots")
    
{
    
    resp0 <- phen
    resp <- resp0
    if (!byclass)
	resp <- rep(1, length(resp0))

    tbl.resp <- table(resp)
    label.col <-
        c("red", "blue", "green", "skyblue", "orange", "pink", "gray20")
    ##label.col <- rainbow(6)
    par(bg = "grey20")

    data <- as.matrix(data)
    dt.cp <- data
    dt <- apply(data, 2, floor.func, cutoff)
    dt <- dt[ ,order(resp) ]
    resp0 <- resp0[ order(resp) ]
    samplenames <- samplenames[ order(resp) ]
    resp <- resp[ order(resp) ]

    start <- 1
    ##mapping order
    ord <- rep(0, length(resp))
    for (i in 1:length(tbl.resp))
    {
	
	ind <- which(resp == i)
	cr <- as.dist(1 - corna(data[ ,ind ]))
	ord[start:sum(tbl.resp[1:i])] <-
            hclust(cr, method = method)$ord + start - 1
	start <- sum(tbl.resp[1:i]) + 1
	
    }
    dt <- dt[ ,ord ]
    resp <- resp[ord]
    resp0 <- resp0[ord]
    samplenames <- samplenames[ord]

    image(x = 1:nrow(dt), y = 1:length(resp), z = dt,
          col = maPalette(low = "red", high = "green", mid = "white",
          k = ncolors), axes = FALSE, xlab = "", ylab = "",
          zlim = c(-cutoff, cutoff))
    ##abline(h=seq(.5, 81.5, b=1), col="gray20", lwd=.2)

    if (showaber)
    {
        ##for (i in 1:nrow(dt))
        ##{
        for (j in 1:ncol(dt))
        {
            
            tmp <- dt.cp[ ,j ]
            i <- which(tmp >= amplif & !is.na(tmp))
            if (length(i) > 0)
                ##if ((!is.na(dt.cp)) && (dt.cp[i,j] >= amplif))
                points(i, rep(j, length(i)), col = "yellow", pch = 20,
                       cex = .7)
            i <- which(tmp <= homdel & !is.na(tmp))
            if (length(i) > 0)
                ##if ((!is.na(dt.cp)) && (dt.cp[i,j] >= amplif))
                points(i, rep(j, length(i)), col = "skyblue",
                       pch = 20, cex = .7)
            
        }
        ##}
    }
    for (j in 1:ncol(dt))
    {

        col <- label.col[resp0[j] + 1]
	mtext(resp0[j], side = 2, at = j, line = .3, col = col,
              cex = .5, las = 2)
	mtext(paste(samplenames[j], ""), side = 4, at = j, line = .3,
              col = col, cex = .25, las = 2)
	
    }
    ##title(main="Whole genome", xlab = "clone", ylab = "sample", col.lab="white", col.main="white")
    title(xlab = "clone", ylab = "sample", col.lab = "white",
          col.main = "white", main = title)

}

plotGeneSign <-
    function(aCGH.obj, phen = rep(1, ncol(aCGH.obj)),
             data = log2.ratios(aCGH.obj), cutoff = 1, ncolors = 50,
             byclass = TRUE, method = "ave", showaber = FALSE, amplif = 1,
             homdel = -1, samplenames = sample.names(aCGH.obj),
             title = "Image Plots", sign = FALSE, dataSign = data,
             nperm = 1000, test = "f", ranks = "y", side = "abs",
             p.thres = c(.01, .05, .1, .2),
             clusterindex = rep(1, nrow(data)))
{

    resp0 <- phen
    resp <- resp0
    if (!(byclass))
	resp <- rep(1, length(resp0))

    tbl.resp <- table(resp)
    label.col <-
        c("red", "blue", "green", "skyblue", "orange", "pink",
          "gray20")
    ##label.col <- rainbow(6)
    ##par(bg="grey20")
    if (sign)
	par(mfrow = c(2, 1))

    data <- as.matrix(data)
    dt.cp <- data
    dt <- apply(data, 2,floor.func, cutoff)    

    dt <- dt[,order(resp)]
    resp0 <- resp0[order(resp)]
    samplenames <- samplenames[order(resp)]
    resp <- resp[order(resp)]

    ##to order within class:

    start <- 1
    ##mapping order
    ord <- rep(0, length(resp))
    for (i in 1:(length(tbl.resp)))
    {
	
	ind <- which(resp == i)
        ##cr <- as.dist(1-cor.na(data[,ind]))
	cr <- dist(t(data[,ind]))
	ord[start:sum(tbl.resp[1:i])] <- hclust(cr, method=method)$ord+start-1
	start <- sum(tbl.resp[1:i])+1
	
	
	
    }
    dt <- dt[,ord]
    resp <- resp[ord]
    resp0 <- resp0[ord]
    samplenames <- samplenames[ord]


    image(x=(1:nrow(dt)), y=1:length(resp), z=dt, col = maPalette(low = "red", high = "green", mid = "white", k =ncolors), axes = FALSE, xlab = "", ylab = "", zlim=c(-cutoff,cutoff))
    ##abline(h=seq(.5, 81.5, b=1), col="gray20", lwd=.2)

    if (showaber)
    {
        ##for (i in 1:nrow(dt))
        ##{
        for (j in 1:ncol(dt))
        {
            
            tmp <- dt.cp[,j]
            i <- (1:length(tmp))[tmp >= amplif & !is.na(tmp)]
            if (length(i) > 0)
                ##if ((!is.na(dt.cp)) && (dt.cp[i,j] >= amplif))
            {
                points(i, rep(j, length(i)), col="yellow", pch=20, cex=.7)
            }
            i <- (1:length(tmp))[tmp <= homdel & !is.na(tmp)]
            if (length(i) > 0)
                ##if ((!is.na(dt.cp)) && (dt.cp[i,j] >= amplif))
            {
                points(i, rep(j, length(i)), col="skyblue", pch=20, cex=.7)
            }
            
        }
        ##}
    }
    for (j in 1:ncol(dt))
    {
	mtext((resp0)[j], side = 2, at = j, line=.3, col=label.col[((resp0)[j]+1)], cex=.5, las=2)
	mtext(paste((samplenames)[j], ""), side = 4, at = j, line=.3, col=label.col[((resp0)[j]+1)], cex=.25, las=2)
	
    }

    if (length(unique(clusterindex)) > 1)
    {
	clusterloc <- cumsum(table(clusterindex))+.5
	clusterloc <- clusterloc[-length(clusterloc)]
	abline(v=clusterloc, col="blue", lwd=.5)
    }
    

    title(xlab = "gene", ylab = "sample", col.lab="black", col.main="black", main=title)

##################now, significance:
    if (sign)
    {
	pal <- c("red", "green", "yellow", "blue")
	pal <- pal[1:length(p.thres)]
	
	res <-  mt.maxT(X=dataSign, classlabel=phen,test=test,side=side,fixed.seed.sampling="y",B=nperm, na=.mt.naNUM, nonpara=ranks)
	maxT <- res$adjp[order(res$index)]	
	
        ##rawp <- res$rawp[order(res$index)]
	teststat <- abs(res$teststat[order(res$index)])
	st <- rep(NA, length(p.thres))
	for (s in 1:length(p.thres))
	{
            if (length(maxT[maxT<=p.thres[s]]) > 0)
            {
                st[s] <- min(teststat[maxT<=p.thres[s]])
            }
	}
        
	st.now <- st
	pal.now <- pal
	par(xaxs="i")
	plot(1:length(teststat),teststat, col="blue", ylim=c(0,max(teststat)), type="h", xlab="gene", ylab="gene statistic", pch=18, col.lab="black", col.axis="black")
	
	if (length(st.now) > 0)
	{
            abline(h=rev(st.now), col=rev(pal.now), lty=2)
	}
	
    }
    if (length(unique(clusterindex)) > 1)
    {
        abline(v=clusterloc, col="red", lwd=.5)
    }

}
