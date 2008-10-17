
maPalette <-
    function (low = "white", high = c("green", "red"), mid = NULL,
              k = 50)
{
    
    low <- col2rgb(low) / 255
    high <- col2rgb(high) / 255
    if (is.null(mid))
    {
        
        r <- seq(low[1], high[1], len = k)
        g <- seq(low[2], high[2], len = k)
        b <- seq(low[3], high[3], len = k)
        
    }
    if (!is.null(mid))
    {
        
        k2 <- round(k / 2)
        mid <- col2rgb(mid) / 255
        r <-
            c(seq(low[1], mid[1], len = k2),
              seq(mid[1], high[1], len = k2))
        g <-
            c(seq(low[2], mid[2], len = k2),
              seq(mid[2], high[2], len = k2))
        b <-
            c(seq(low[3], mid[3], len = k2),
              seq(mid[3], high[3], len = k2))
        
    }
    rgb(r, g, b)
    
}


##########


##maColorBar <-
##    function (x, horizontal = TRUE, col = heat.colors(50),
##              scale = 1:length(x), k = 10, ...)
##{
    
##    if (is.numeric(x))
##        colmap <- col
##    else
##    {
        
##        colmap <- x
##        low <- range(scale)[1]
##        high <- range(scale)[2]
##        x <- seq(low, high, length = length(x))
        
##    }
##    if (length(x) > k)
##        x.small <- seq(x[1], x[length(x)], length = k)
##    else x.small <- x
##    if (horizontal) {
##        image(x, 1, matrix(x, length(x), 1), axes = FALSE, xlab = "",
##              ylab = "", col = colmap, ...)
##        axis(1, at = rev(x.small), labels = signif(rev(x.small),
##                                   2), srt = 270, col="white")
##    }
##    if (!horizontal) {
##        image(1, x, matrix(x, 1, length(x)), axes = FALSE, xlab = "",
##              ylab = "", col = colmap, ...)
##        par(las = 1)
##        axis(4, at = rev(x.small), labels = signif(rev(x.small),
##                                   2), col="white")
##        par(las = 0)
##    }
##    box()
##}

##########################################################################
############################OLD::::
##########################################################################

##rgcolors.my.func <-
##    function (n = 50)
##{
    
##    k <- round(n / 2)
##    g <- c(rep(0, k), seq(0, 1, length = k))
##    r <- c(rev(seq(0, 1, length = k)), rep(0, k))
##    rgb(r, g, rep(0, 2 * k))

##}

#################################################################################
##################################################################################


##SHOULD it BE (0, lim) rather tha (0,1)

##auxilliary function for image

##rgcolors.my.black.func <-
##    function (n = 50, lm = 1)
##{
    
##    k <- round(n/2)
##    g <- c(rep(0, k), seq(0, lm, length = k))
##    r <- c(rev(seq(0, lm, length = k)), rep(0, k))
##    b <- rep(0, 2 * k)
##    rgb(r, g, b)

##}



##################################################################################
##################################################################################

##auxilliary function for image

##colorbar.my.black.func <- function(n=50, lim=1, y=1, coloraxis="black")
##{
##    mx <- lim
##    mn <- -mx
##    x<-seq(mn,mx,length=n)
##    image(1:length(x),y,matrix(x,length(x),1),col=rgcolors.my.black.func(n,lim),axes=FALSE,xlab="",ylab="")
##    ##axis(2,at=length(x):1,labels=rev(x),srt=270)
##    ##axis(1,at=length(x):1,labels=rev(x),srt=270, col.axis=coloraxis)
##    box()
##}

##if we want to create image:

##postscript("colorrange.ps", paper="letter")
##par(bg="grey20")    
##colorbar.my.black.func(lim=1,n=20, coloraxis="white")    
##box()
##dev.off()


##postscript("colorrange.ps", paper="letter")
##par(bg="grey20")    
##lm <- 1
##k <- 50
##maColorBar(x=seq(-lm,lm,len=k), col = maPalette(low="red", high="green", mid="white", k=k), h=T)  
##dev.off()
##


##old color-bar function:

##colorbar.my.black.func <- function(lim=1, y=1, coloraxis="black")
##{
##  max <- lim
##  min <- -max
##  n <- (1/(max-min))*100
##  step <- 10*max/(2*50)
##  x<-seq(min,max,by=step)
##  image(1:length(x),y,matrix(x,length(x),1),col=
##rgcolors.my.black.func(n),axes=FALSE,xlab="",ylab="")
##  ##axis(2,at=length(x):1,labels=rev(x),srt=270)
##  axis(1,at=length(x):1,labels=rev(x),srt=270, col.axis=coloraxis)
##  box()
##}

########################OLD ENDED########################################
##################################################################################
##################################################################################

##flooring function. auxilliary function for image

##floorFunc <- function(x, floor)
##{
##    x[x > floor & !is.na(x)] <- floor
##    x[x < -floor & !is.na(x)] <- -floor
##    x
##}


##################################################################################
##################################################################################

##create image plots per chromosomes

##data- p by n matrix of lograt. missing values are ok (codede as NA)
##phen- phenotype appearing on the left side of the image. no missing values are allowed
##datainfo -- p by >= 3 matrix. Has to contain columns Clone, Chrom and lb
##chrominfo -- file describing lengths of chromosomes and centromere location
##chrom -- which chromosomes to show
##cutoff- where to cut-off the values too. Lower values make image colors brighter. 
##amplif -- anything >= is marked with yellow dots
##homde -- anything <= is marked with light blue dots
##bycllass-- whether to order samples and cluster them inside the class (T) or cluster
##all samples (F). Clustering is done using agglom. hierachical clustering with average
##linkage
##samplenames-- shown on the right. anything can be there
##clonenames -- names of the clones. Shown on the bottom
##ttl -- title of the image
##filePS -- name of the PS file

##Note that clustering is done based on individual chromosomes and order of ssamples
##varies from chromosome to chromosome

plotvalChrom.func <-
    function(data, phen=rep(1, ncol(data)), datainfo=clones.info,
             chrominfo=human.chrom.info.Jul03, chrom=1:20, cutoff=1, ncolors=50,
             amplif=1, homdel=-1, byclass=TRUE,
             samplenames=dimnames(data)[[2]],
             clonenames=datainfo$Clone, ttl="Image Plot",
             filePS="plotvalschrom.ps", ps=TRUE)
{
    ##label.col <- c("red", "green", "blue", "yellow", "skyblue", "orange", "pink", "gray20")
    label.col <- rainbow(6)
    if (ps)
    {
	postscript(filePS, paper="letter")
    }
    else
    {
	pdf(filePS, width=11, height=8.5)
    }
    par(bg="grey20")    
    samplenames.cp <- samplenames
    for (chr in chrom)
    {   
        resp0 <- phen
        resp <- resp0
        samplenames <- samplenames.cp
        if (!(byclass))
        {
            resp <- rep(1, length(resp0))
        }

        tbl.resp <- table(resp)
        kb <- datainfo$kb[datainfo$Chrom==chr]
        dt <- as.matrix(data[datainfo$Chrom==chr,])
        clonenms <- clonenames[datainfo$Chrom==chr]


        dt.cp <- dt
        dt <- apply(dt.cp, 2,floorFunc, cutoff)       
        if (chrominfo$centr[chr] >0)
        {
            centr <- length(kb[kb<=chrominfo$centr[chr]])
        }
        dt <- dt[,order(resp)]
        dt.cp <- dt.cp[,order(resp)]
        resp0 <- resp0[order(resp)]
        samplenames <- samplenames[order(resp)]
        resp <- resp[order(resp)]


        start <- 1

        ##mapping order
        ord <- rep(0, length(resp))
        for (i in 1:(length(tbl.resp)))
        {
            
            ind <- which(resp == i)
            cr <- as.dist(1-cor.na(dt.cp[,ind]))
            ord[start:sum(tbl.resp[1:i])] <- hclust(cr, method="ave")$ord+start-1
            start <- sum(tbl.resp[1:i])+1
        }

        dt <- dt[,ord]
        dt.cp <- dt.cp[,ord]
        resp <- resp[ord]
        resp0 <- resp0[ord]
        samplenames <- samplenames[ord]

        image(x=(1:length(kb)), y=1:length(resp), z=dt, col = maPalette(low = "red", high = "green", mid = "white" ,k=ncolors), axes = FALSE, xlab = "", ylab = "", zlim=c(-cutoff, cutoff))

        if (chrominfo$centr[chr] >0)
        {
            abline(v=centr, col="black")
        }
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
            mtext(paste(clonenms[i], ""), side = 1, at = i, line=.3, col="white", cex=.25, las=2)
            for (j in 1:ncol(dt.cp))
            {
		if (i ==1)
		{
                    mtext((resp0)[j], side = 2, at = j, line=.3, col=label.col[((resp0)[j]+1)], cex=.5, las=2)
                    mtext(paste((samplenames)[j], ""), side = 4, at = j, line=.3, col=label.col[((resp0)[j]+1)], las=2, cex=.5)
                    
		}
		if (!is.na(dt.cp[i,j]) && dt.cp[i,j]>=amplif)
		{	
                    
                    points(i, j, col="yellow", pch=20, cex=.7)
		}
		if (!is.na(dt.cp[i,j]) && dt.cp[i,j] <= homdel)
		{	
                    
                    points(i, j, col="skyblue", pch=20, cex=.7)
		}
		
            }
        }
        title(main=paste(ttl, " Chromosome ", chr), col.main="white")


    }
    dev.off()
}

##################################################################################
##################################################################################

#############################################################################

##creates the color scale graph
##k=50
##lm=1
##maColorBar(x=seq(-lm,lm,len=k), col = maPalette(low="red", high="green", mid="white", k=k), h=TRUE) 
##postscript("colorrange.ps", paper="letter")
##par(bg="grey20")    

##dev.off()
##
##
#############################################################################

##################################################################################
################################################################################
##################################################################################
##################################################################################

##plots images of correlation matrices

##X -- is p by n data matrix or p by p correlation matrix 
##new =TRUE if correlation matrix and F if data matrix and correlation needs to be computed
##

##plot.my.cor <-
##    function (X, new = TRUE, nrgcols = 50, labels = FALSE, labcols = 1,
##              title = "", ...)
##{
    
##    n <- ncol(X)
##    corr <- X
##    if (new)
##        corr <- cor.na(X)
##    image(1:n, 1:n, corr[, n:1], col = maPalette(low = "red", high = "green", mid = "black", k=nrgcols),
##          axes = FALSE, xlab = "", ylab = "",...)
##    if (length(labcols) == 1) {
##        axis(2, at = n:1, labels = labels, las = 2, cex.axis = 0.6,
##             col.axis = labcols)
##        axis(3, at = 1:n, labels = labels, las = 2, cex.axis = 0.6,
##             col.axis = labcols)
##    }
##    if (length(labcols) == n) {
##        cols <- unique(labcols)
##        for (i in 1:length(cols)) {
##            which <- (1:n)[labcols == cols[i]]
##            axis(2, at = (n:1)[which], labels = labels[which],
##                 las = 2, cex.axis = 0.6, col.axis = cols[i])
##            axis(3, at = which, labels = labels[which], las = 2,
##                 cex.axis = 0.6, col.axis = cols[i])
##        }
##    }
##    mtext(title, side = 3, line = 3)
##    box()
    
##}
################################

##################################################################################
##################################################################################

##averages correlation in a given bin. 

##cor.all -- p by p matrix of cross-correlations
##bin -- size of the bin (in kb)
##datainfo -- as before
##chrominfo -- as before
##lm -- the floor

##note that each chromosome is subdivided into the equasized number of bins
##so that their size is as close to specified bin as possible. The resulting
##bin sizes will slightly vary across chromosomes


##cor.bin.func <-
##    function(cor.all, bin, datainfo, chrominfo, lm=1)
##{
    
##    bin.num <- round( chrominfo$length/bin)
##    for (i in 1:length(bin.num))
##    {
##        if (bin.num[i] == 0)	
##        {
##            bin.num[i] <- 1
##        }
##    }
##    bin.size <- chrominfo$length/bin.num
    
##    cor.bin <- matrix(0, nrow=sum(bin.num), ncol=sum(bin.num))
    
##    chrom.bound <-  c(0,cumsum(bin.num))
##    chrom.bound.mid <- chrom.bound[-length(chrom.bound)]+diff(chrom.bound)/2
##    chrom.bound <- chrom.bound[2:(length(chrom.bound)-1)]
    

##    kb <- datainfo$kb
##    kbCum <- datainfo$kbCum
##    chr <- datainfo$Chrom
    
##    chruniq <- unique(chr)	
    

##    mx <- 1
    
##    for (ix in 1:length(chruniq))
##    {
        
##        for (jx in 1:bin.num[ix])
##        {
##            startx <- bin.size[ix]*(jx-1)
##            endx <- bin.size[ix]*jx
##            indx <- (1:length(kb))[kb >=startx & kb < endx & chr ==ix]
##            my <- 1
##            for (iy in 1:length(chruniq))
##            {
                
##                for (jy in 1:bin.num[iy])
##                {
##                    starty <- bin.size[iy]*(jy-1)
##                    endy <- bin.size[iy]*jy
##                    indy <- (1:length(kb))[kb >=starty & kb < endy & chr ==iy]
##                    if ((length(indx) > 0) && (length(indy) > 0))
##                    { 
##                        cor.mean <- mean(cor.all[indx, indy], na.rm=TRUE)
##                        if (cor.mean > lm)
##                        {
##                            cor.mean <- lm
##                        }
##                        if (cor.mean < -lm)
##                        {
##                            cor.mean <- -lm
##                        }
                        
##                    }
##                    else
##                    {
##                        cor.mean <- NA
##                    }
##                    cor.bin[mx,my] <- cor.mean
                    
##                    my <- my+1
##                }
##            }
##            mx <- mx+1
##        }
##    }
##    list(cor.bin=cor.bin, chrom.bound=chrom.bound, chrom.bound.mid=chrom.bound.mid)
##}


##################################################################################
##################################################################################

##plots genome

##plotGenome0.func <- function(sample, yScale = c(-2,2), namePSfile = "try.ps", data=Log2Rat, datainfo=clones.info, chrominfo=chrom.info, samplenames=sampleNames, naut = 22, X=TRUE, Y=TRUE, ps=TRUE, plotend=TRUE)
##{

##    nsamples <- length(sample)
##    ord <- order(datainfo$Chrom, datainfo$kb)
##    chrom <- datainfo$Chrom[ord]
##    kb <- datainfo$kb[ord]
##    data <- data[ord,]

##    chrom.rat <- chrominfo$length/max(chrominfo$length)
##    chrom.start <- rep(0, nrow(chrominfo))
##    for (i in 2:length(chrom.start))
##    {
##        chrom.start[i] <- sum(chrominfo$length[1:(i-1)])
##    }
##    ##
##    ##
##    ##chrom.mid contains middle positions of the chromosomes relative to
##    ##the whole genome (useful for plotting the whole genome)
##    ##
##    chrom.mid <- rep(0, nrow(chrominfo))
##    for (i in 1:length(chrom.start))
##    {
##        chrom.mid[i] <- chrom.start[i]+chrominfo$length[i]/2
##    }
##    if (plotend)
##    {
##        ##start a postscript file
##	if (ps)
##	{
##            postscript(namePSfile, paper="letter")
##	}
##	else
##	{
##            pdf(namePSfile, height=8.5, width=11)
##	}

##	par(mfrow=c(nsamples,1))
##    }
##    par(cex=.6, pch=18, lab=c(1,6,7), cex.axis=1.5, xaxs="i")
##    for (k in 1:nsamples)
##    {

##        vec <- data[,sample[k]]
##        name <- samplenames[sample[k]]

##        clone.genomepos <- rep(0, length(kb))
##        for (i in 1:nrow(chrominfo))
##        {
##            clone.genomepos[chrom==i] <- kb[chrom==i]+chrom.start[i]
##        }
######################

######################
##        ##Now, determine vertical scale for each chromosome:

##        y.min <- rep(yScale[1], nrow(chrominfo))
##        y.max <- rep(yScale[2], nrow(chrominfo))
        
##        for (i  in 1:nrow(chrominfo))
##        {
##            if (min.na(vec[(chrom==i)]) < y.min[i])
##            {
##                y.min[i] <- min.na(vec[(chrom==i)])
##            }
##            if (max.na(vec[(chrom==i)]) > y.max[i])
##            {
##                y.max[i] <- max.na(vec[(chrom==i)])
##            }
            
##        }
        
##        ##set genome scale to the min and mx values of the rest of the chromosomes:
        
##        ygenome.min <- min.na(y.min)
##        ygenome.max <- max.na(y.max)
        
        
        
###########################
        
##        plot(clone.genomepos/1000, vec, ylim=c(ygenome.min,ygenome.max), xlab="", ylab="", xlim=c(min(clone.genomepos[clone.genomepos>0], na.rm=TRUE)/1000, clone.genomepos[length(clone.genomepos[clone.genomepos>0])]/1000), col="black")
        
##        ##title(main=paste(name, " ", sample[k], " - Whole Genome"), ylab="Log2Ratio", xlab="Chromosome", cex.lab=1.5,cex.main=2)
##        title(main=paste(sample[k], " ", name), ylab="Log2Ratio", xlab="", cex.lab=1.5,cex.main=2)
        
##        for (i in seq(1,naut,b=2))
##        {
##            mtext(paste("", i), side = 1, at = (chrom.mid[i]/1000), line=.3, col="red")
##        }
##        for (i in seq(2,naut,b=2))
##        {
##            mtext(paste("", i), side = 3, at = (chrom.mid[i]/1000), line=.3, col="red")
##        }

##        if (X)
##        {
##            mtext("X", side = 1, at = (chrom.mid[naut+1]/1000), line=.3, col="red")
##        }
##        if (Y)
##        {
##            mtext("Y", side = 3, at = (chrom.mid[naut+2]/1000), line=.3, col="red")
##        }

##        abline(v=c(chrom.start/1000, (chrom.start[nrow(chrominfo)]+chrominfo$length[nrow(chrominfo)])/1000), lty=1)
##        ##abline(h=seq(ygenome.min,ygenome.max, b=.2), lty=3)
##        abline(h=seq(-1,1, b=.5), lty=3)

##        abline(v=(chrominfo$centromere+ chrom.start)/1000, lty=3, col="red")
##    }
##    if (plotend)
##    {
##	dev.off()
##    }
##}


##################################################################################
##################################################################################

##plots chromosomes

##plotChrom0.func <-
##    function(sample, chr, yScale = c(-1,1), namePSfile = "try.ps",
##             data, datainfo, chrominfo=human.chrom.info.Jul03,
##             samplenames=sampleNames)
##{
    
##    nsamples <- length(sample)
##    ord <- order(datainfo$Chrom, datainfo$kb)
##    chrom <- datainfo$Chrom[ord]
##    kb <- datainfo$kb[ord]
##    data <- data[ord,]
    
##    ##start a postscript file
    
##    postscript(namePSfile, paper="letter")
##    par(mfrow=c(nsamples,1))
##    par(cex=.6, pch=18, lab=c(1,6,7), cex.axis=1.5)
##    kb <- kb[chrom==chr]
##    centr.loc <- chrominfo$centromere[chr]
##    for (k in 1:nsamples)
##    {
        
##        vec <- data[chrom==chr,sample[k]]
##        name <- samplenames[sample[k]]
        
######################
##        ##Now, determine vertical scale for each chromosome:
        
##        y.min <- min(yScale[1],min.na(vec))
##        y.max <- max(yScale[2],max.na(vec))
        
        
###########################
        
##        plot(kb/1000, vec, ylim=c(y.min,y.max), xlab="", ylab="", xlim=c(min(kb[kb>0], na.rm=TRUE)/1000, kb[length(kb[kb>0])]/1000), col="black", cex=1.5)
##        lines(kb/1000, vec, col="blue",lwd=.5)
##        abline(v=centr.loc/1000, col="red", lty=2)
##        abline(h=0, col="black", lty=2)
##        abline(h=seq(-.6,.6, b=.2), lty=3)
##        title(main=paste(name, " Chr ",chr), ylab="Log2Ratio", xlab="Chromosome", cex.lab=1.5,cex.main=2)
        
##    }
##    dev.off()
##}

##################################################################################
##################################################################################

##creating chromFull.info file

##plot entire genome or chrom. or all chrom's:

##chrom.rat <- chrom.info$length/max(chrom.info$length)  
##i.e. for each chromosome it repreesents the fraction of length of the
##longest chromosome
##
##chrom.start contains starting positions of the chromosomes relative to the
##whole genome (0 for the first)
##chrom.start <- rep(0, 23)
##for (i in 2:length(chrom.start))
##{
##	chrom.start[i] <- sum(chrom.info$length[1:(i-1)])
##}
##
##chrom.mid contains middle positions of the chromosomes relative to
##the whole genome (useful for plotting the whole genome)
##chrom.mid <- rep(0, 23)
##for (i in 1:length(chrom.start))
##{
##	chrom.mid[i] <- chrom.start[i]+chrom.info$length[i]/2
##}
##
##chromFull.info <- as.data.frame(cbind(chrom.info, chrom.start, chrom.mid, chrom.rat))
##dimnames(chromFull.info)[[2]] <- c("chr", "length", "centromere", "start", "mid", "rat")
##################################################################################
##################################################################################

##plots
##either all chromsoomes and genome/whole genome/or individual chromosome
##data -- as before
##map -- clones.info
##samplename -- index or samplename of the smaple to plot
##sampNm -- vector of sample names
##whatToPlot: "All" -- all chromosomes and Gneome/ "G" -- Genome or chromosomal number

plotCGH.func <-
    function (data=data.cgh, map=map.cgh, chrominfo=human.chrom.info.Jul03,
              samplename, sampNm=sampleNames, whatToPlot= "All",
              yScale = c(-2,2), namePSfile = "plotCGH.ps")
{
################General Comments############################################
    ##Jane Fridlyand, 08/13/2001
    ##
    ##Function to plot CGH data so that the horizontail box sizes are proportinal
    ##to the physical lengths of chromosomes###
    ##
    ##Lots of comments here (required files etc)
    ##
##############################################################################
###############################NOTE##################################
    ##
    ##This function does no trouble shooting -- when one of the necessary files
    ##is in the wrong format and can not be read in or just not available,
    ##it fails.
##################End of NOTE#######################################
######################################################################
#################ARGUMENTS############################################
    ##
    ##Parametesr that probably should not be controlled by an average user:
    ##
    ##1. data: name of the data file containing intensity ratio data
    ## default is data.cgh
    ##samples are in columns and clones are in rows
    ##
    ##2. map: name of the file containing 2 columns: chr and kb
    ##default is map.cgh
    ##
    ##3. chrominfo -- name of the file containing 6 columns: chr, length, 
    ##   centromere, start, mid, rat (ratio w.r.t to the longest chrom.
    ##
#########################################
    ##
    ##Required:
    ##
    ##1. samplename or index of a sample to plot 
    ##
    ##
##########Parameters that should be controlled by an average user:
    ##
    ##1. whatToPlot: to plot
    ##    a) 23 chromosomes and whole genome (All) -- default
    ##    b) 1 chromosome only (number of a chromosome) 
    ##    c) whole genome only (G)
    ##
    ##
    ##2. yScale: 
    ##   number to which scale is fixed. default is is c(-2,2)
    ##   scales are adjusted for chromsomes which are outside that range
    ##
    ##3. namePSfile -- gets used only if PS file is to be created: -- name 
    ##   default is plotCGH.ps 
    ##
#######################################################################
#########################################################################
#########################EXAMPLES of USAGE #############################
    ##
    ##Suppose a CGH  file "/path/filename" needs to be analyzed.
    ##
###
    ##Vanilla Example:
    ##
    ##to plot 23 chromosomes and genome with vertical scale of (-2,2) and produce
    ## output ps file "plotCGH.ps" for a 5-th sample
    ##
    ##plotCGH.func(samplename=5)
    ##
###
    ##to plot 23 chromosomes and genome with vertical scale of (-1,1) 
    ##and produce output ps file "path/new.ps" for a samplename "mysample.txt"
    ##
    ##plotCGH.func(samplename="mysample.txt",yScale=c(-1,1),namePSfile ="/path/new.ps")
    ##
####
    ##plot chromosome 1 only with vertical scale (-2,2) for a smaplename 6:
    ##
    ##plotCGH.func(samplename=6,whatToPlot = 1)
    ##
###
    ##plot whole Genome only:
    ##
    ##plotCGH.func(samplename=6,whatToPlot = "G")
    ##
#########################END of EXAMPLES of USAGE ############################
#######################################################################
#########################################################################

#####################START#################################################

#########creating chromFull.info file

    chrom.rat <- chrominfo$length/max(chrominfo$length)  
    ##i.e. for each chromosome it repreesents the fraction of length of the
    ##longest chromosome
    ##
    ##chrom.start contains starting positions of the chromosomes relative to the
    ##whole genome (0 for the first)
    chrom.start <- rep(0, 23)
    for (i in 2:length(chrom.start))
    {
	chrom.start[i] <- sum(chrominfo$length[1:(i-1)])
    }
    ##
    ##chrom.mid contains middle positions of the chromosomes relative to
    ##the whole genome (useful for plotting the whole genome)
    chrom.mid <- rep(0, 23)
    for (i in 1:length(chrom.start))
    {
	chrom.mid[i] <- chrom.start[i]+chrominfo$length[i]/2
    }

    chromFull.info <- as.data.frame(cbind(chrominfo, chrom.start, chrom.mid, chrom.rat))
    dimnames(chromFull.info)[[2]] <- c("chr", "length", "centromere", "start", "mid", "rat")

########################################


    ##computing positions in genome for each clone:

    clone.genomepos <- rep(0, length(map$kb))
    for (i in 1:23)
    {
	clone.genomepos[map$Chrom==i] <- map$kb[map$Chrom==i]+chromFull.info$start[i]
    }

##########
    ##Now, determine vertical scale for each chromosome:

    y.min <- rep(yScale[1], 23)
    y.max <- rep(yScale[2], 23)

##############
    ##figure out the sample
    ##
    smpnames <- sampNm
    if ((samplename >= 1) && (samplename <= nrow(data)))
        ##samplename was the index
    {
	smp <- samplename
	samplename <- smpnames[smp]
        ##so now samplename is a name
    }
    else #samplename 
    {
	
	smp <- (1:length(smpnames))[smpnames==samplename]
    }
##############
    ##values to plot:

    vals <- data[,smp]

#############
    ##adjust scales of chromosomes that have values outside a fixed scale

    for (i in 1:23)
    {
	y.min[i] <- min(c(vals[map$Chrom==i],yScale[1]), na.rm=TRUE)
	y.max[i] <- max(c(vals[map$Chrom==i],yScale[2]), na.rm=TRUE)
    }

    ##set genome scale to the max and min values across chrom's

    ygenome.min <- min(y.min, na.rm=TRUE)
    ygenome.max <- max(y.max, na.rm=TRUE)	

#########################
    ##start a postscript file

##########start plotting:

    ##plot one chromosome only:

    if ((whatToPlot >= 1) && (whatToPlot <= 23)) 
    {
	postscript(namePSfile, paper="letter")
	par(lab=c(15,6,7), pch=18, cex=1, lwd=1)
	plot(map$kb[map$Chrom==whatToPlot]/1000, vals[map$Chrom==whatToPlot], ylim=c(y.min[whatToPlot],y.max[whatToPlot]), xlab="", ylab="", col="black", xlim=c(0, chromFull.info$length[whatToPlot]/1000))
	lines(map$kb[map$Chrom==whatToPlot]/1000, vals[map$Chrom==whatToPlot],col="blue")
	title(main=paste("Sample ", samplename, " ", smp, " Chr ",whatToPlot), xlab="kb (in 1000's)", ylab="Log2Ratio")
	abline(h=seq(-.6,.6, b=.2), lty=3)
	abline(h=0, col="black", lty=2)
	abline(v=chromFull.info$centromere[whatToPlot]/1000, lty=2, col="red")
	dev.off()
    }

    ##if plot whole genome only:

    else if (whatToPlot == "G")
    {
	postscript(namePSfile, paper="letter")
	par(cex=.6, pch=18, lab=c(1,6,7), cex.axis=1.5, xaxs="i")
	plot(clone.genomepos[map$Chrom<=23]/1000, vals[map$Chrom<=23], ylim=c(ygenome.min,ygenome.max), xlab="", ylab="", xlim=c(min(clone.genomepos[clone.genomepos>0], na.rm=TRUE)/1000, clone.genomepos[length(clone.genomepos[clone.genomepos>0])]/1000), col="black")
        title(main=paste("Sample ", samplename, " ", smp, "Whole Genome"), ylab="Log2Ratio", xlab="Chromosome", cex.lab=1.5,cex.main=2)
	for (i in seq(1,21,b=2))
	{	
            mtext(paste("", i), side = 1, at = (chromFull.info$mid[i]/1000), line=.3, col="red")
	}
	mtext("X", side = 1, at = (chromFull.info$mid[length(chromFull.info$mid)]/1000), line=.3, col="red")
	abline(v=c(chromFull.info$start/1000, (chromFull.info$start[23]+chromFull.info$length[nrow(chrominfo)])/1000), lty=1)
	abline(h=seq(ygenome.min,ygenome.max, b=.2), lty=3)
	abline(v=(chromFull.info$centromere+ chromFull.info$start)/1000, lty=3, col="red")
	dev.off()
    }
    else #if all chromosomes and genome are plotted:
    {
	postscript(namePSfile, paper="letter", horizontal=FALSE)
        ##just a safety line
	close.screen(all=TRUE)
        ##"inch" factor for to determine size of the plot in inches (for "pin" parameter)
	fact <- 3.9
        ##split the screen

	split.screen(c(9,1))
	screen(1)
	split.screen(c(1,2))
	screen(2)
	split.screen(c(1,2))
	screen(3)
	split.screen(c(1,2))
	screen(4)
	split.screen(c(1,2))
	screen(5)
	split.screen(c(1,3))
	screen(6)
	split.screen(c(1,3))
	screen(7)
	split.screen(c(1,4))
	screen(8)
	split.screen(c(1,3))
	screen(28)
	split.screen(c(1,2))
	screen(29)
	split.screen(c(1,2))

        ##plot chromosomes

	scr.seq <- c(10:27, 31:34, 30)  
	j.seq <- 1:23
	for (j in j.seq)
	{
            
            screen(scr.seq[j])
            par(cex=.5, pch=20, lab=c(15,4,7), tcl=-.2, las=1, oma=c(0,0,0,0), cex.axis=1.3, cex.main=1.3, mgp=c(0,.15,0), lwd=.5)
            par(pin=c(chromFull.info$rat[j]*fact, .65))
            
            plot(map$kb[map$Chrom==j]/1000, vals[map$Chrom==j], ylim=c(y.min[j],y.max[j]), xlab="", ylab="", type="l", col="blue", xlim=c(0, chromFull.info$length[j]/1000))
            points(map$kb[map$Chrom==j]/1000, vals[map$Chrom==j], col="black")
            
            if (j < 23)
            {
		title(main=paste("Chr",j), line=.1)
            }
            else
            {
                title(main="Chr. X", line=.1)
            }
            abline(h=seq(y.min[j],y.max[j], b=.5), lty=3)
            abline(v=0, lty=2)		
            abline(v=chromFull.info$centromere[j]/1000, lty=2, col="red")
            
	}
	
        ##plot genome:
	screen(9 )

	par(cex=.5, pch=20, lab=c(1,4,7), tcl=-.2, las=1, cex.axis=1.3, mgp=c(0,.15,0), cex.main=1.3, xaxs="i")
	par(pin=c(7.8, .55))
	plot(clone.genomepos[map$Chrom<=23]/1000, vals[map$Chrom<=23], ylim=c(ygenome.min,ygenome.max), xlab="", ylab="", xlim=c(min(clone.genomepos[clone.genomepos>0], na.rm=TRUE)/1000, clone.genomepos[length(clone.genomepos[clone.genomepos>0])]/1000), col="black", type="l", lwd=1)
	title(main="Whole Genome (not to horizontal scale)",line=.1)
	for (i in seq(1,21,b=2))
	{	
            mtext(paste("", i), side = 1, at = (chromFull.info$mid[i]/1000), line=.3, col="red", cex.main=.5)
	}
	mtext("X", side = 1, at = (chromFull.info$mid[nrow(chromFull.info)]/1000), line=.3, col="red",cex.main=.5)
	abline(v=c(chromFull.info$start/1000, (chromFull.info$start[23]+chromFull.info$length[nrow(chromFull.info)])/1000), lty=1)
	abline(h=seq(ygenome.min,ygenome.max, b=.5), lty=3)
	abline(v=(chromFull.info$centromere+chromFull.info$start)/1000, lty=3, col="red")

        mtext(paste("Sample ", samplename, " ", smp, "Log2Ratio of Intensities vs Position in 1000's kb"), outer=TRUE, line=-1.2, cex=.8)
	dev.off()
    }
####################
    

}

#################################################################################
#################################################################################
##auxilliary function for frequency plot

##gainLoss <- function (dat, cols,thres, quant=.5)
##{

##if (length(thres) == 1)
##{

##	thres <- rep(thres, ncol(dat))
##}
##if (length(thres) != ncol(dat))
##{
##	print("Error: number of thresholds is not the same as number of tumors")
##	exit()
##}

##dt <- as.matrix(dat[,cols])
##thr <- thres[cols]
##gain <- rep(0, nrow(dt))
##gain.med <- gain
##loss <- rep(0, nrow(dt))
##loss.med <- loss

##for (i in 1:nrow(dt))
##{
##	x <- dt[i,]
##	th <- thr[!is.na(x)]
##	x <- x[!is.na(x)]
##	tmp.gain <- x-th
##	tmp.loss <- x+th
##	gain[i] <- length(tmp.gain[tmp.gain>=0])/length(x)
##	if (gain[i] > 0)
##	{
##		#gain.med[i] <- median(x[tmp.gain>=0])
##		gain.med[i] <- quantile(x[tmp.gain>=0], 1-quant)
##	}
##	loss[i] <- length(tmp.loss[tmp.loss<=0])/length(x)
##	if (loss[i] > 0)
##	{
##		#loss.med[i] <- median(x[tmp.loss<=0])
##		loss.med[i] <- quantile(x[tmp.loss<=0], quant)
##	}
##}

##list(gainP=gain, lossP=loss, gainMed=gain.med, lossMed=loss.med)

##}
#################################################################################
#################################################################################



##################################################################################
##################################################################################

##frequency plot for the whole genome

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
##nperm -- if sign =TRUE, then how many permutations for maxT
##test -- name of the test
##ranks -- whether to work with ranked data If nor, "n"
##side -- two sisded (abs) test or 1-sided ("upper" or "lower")
##p.thres -- for which adjusted p-values show the cut-off
##filePS -- name of the ps/pdf file
##PS = "ps" or "pdf"
##X=TRUE -- is X (23) chrom included?
##Y=FALSE  -- is Y (24) chrom included?
##numaut=22 -- number of autosomes

plotfreq.stat.final.func <-
    function(data, rsp,datainfo, chrominfo, titles, thres=.2,cutplot =
             0, ylm=c(-1,1), sign=FALSE, dataSign=data, nperm=1000,
             test="f", ranks="y", side="abs",
             p.thres=c(.01,.05,.1,.2), X=TRUE, Y=FALSE, numaut=22,
             filePS="gainslosses.ps", PS="ps", onepage=TRUE)
{

    rsp.uniq <- unique(rsp)
    colmatr <- matrix(0, nrow=length(rsp.uniq), ncol=length(rsp))
    for (i in 1:nrow(colmatr))
    {
	colmatr[i,rsp==rsp.uniq[i]] <- 1
    }



    pal <- c("red", "blue", "green", "yellow")


    if (nrow(colmatr) == 1)
    {
	sign <- F
    }

    nr <- nrow(colmatr)
    if (sign)
    {	
	nr <- nr+1
	
    }

    tmp <- matrix(0, ncol=2,nrow=1)   
    tmp <- as.data.frame(tmp) 
    dimnames(tmp)[[2]] <- c("gainP", "lossP")   
    gainloss <- rep(list(tmp),nrow(colmatr))            

    for (j in 1:nrow(colmatr))
    {
	
	cols <- (1:ncol(colmatr))[colmatr[j,]==1]
	gainloss[[j]] <- gainLoss(dat=data, cols=cols,thres=thres)
	
    }

    if (sign)
    {
	dt <- dataSign[,colmatr[1,]==1]
	rsp <- rep(1, ncol(dt))
	for (j in 2:nrow(colmatr))
	{
            dt <- cbind(dt, dataSign[,colmatr[j,]==1])
            rsp <- c(rsp, rep(j, ncol(dataSign[,colmatr[j,]==1])))
	}
	rsp <- rsp-1
	res <-  mt.maxT(X=dt,classlabel=rsp,test=test,side=side,fixed.seed.sampling="y",B=nperm, na=.mt.naNUM, nonpara=ranks)
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
        
        ##st.now <- st[!is.na(st)]
        ##pal.now <- pal[!is.na(st)]

	st.now <- st
	pal.now <- pal
	
    }

    numchr <- numaut
    if (X)
    {
	numchr <- numchr+1
    }
    if (Y)
    {
	numchr <- numchr+1
    }

    chrominfo <- chrominfo[1:numchr,]

    ##compute cumulative kb locations
    start <- c(0, cumsum(chrominfo$length))
    kb.loc <- datainfo$kb
    for (i in 1:numchr)
    {
        tmp <- start[i]+datainfo$kb[datainfo$Chrom==i]
        kb.loc[datainfo$Chrom==i] <- tmp
    }
    ##preparation for graphs
    chrom.start <- rep(0,numchr)
    for (i in 2:length(chrom.start))
    {
        chrom.start[i] <- sum(chrominfo$length[1:(i-1)])

    }
    chrom.centr <- rep(0,numchr)
    for (i in 1:length(chrom.centr))
    {
        chrom.centr[i] <- chrom.start[i]+chrominfo$centr[i]

    }

    chrom.mid <- rep(0, numchr)
    for (i in 1:length(chrom.start))
    {
        chrom.mid[i] <- chrom.start[i]+chrominfo$length[i]/2
    }

    ##now, plot
    ##nc <- max(length(titles)/2,1)
    nc <- 1
    
    if (PS == "ps")
    {
        postscript(filePS,paper="letter")
    }
    else if (PS == "pdf")
    {
        pdf(filePS, width = 8.5, height =11)
    }
    if (onepage)
    {
        par(mfrow=c(nr,nc), lab=c(1,8,7), tcl=-.2,  xaxs="i")
    }
    else
    {
        par(mfrow=c(1,nc), lab=c(1,8,7), tcl=-.2,  xaxs="i")
    }
    for (g in 1:length(titles))
    {
        gl <- gainloss[[g]]
        tl <- titles[g]
        ylm[1] <- min(ylm, min(gl$lossP))
        ylm[2] <- max(ylm, max(gl$gainP))
        
        plot(kb.loc[gl$gainP>=cutplot],gl$gainP[gl$gainP>=cutplot], col="green", type="h", xlab="chromosome number", ylab="Fraction gained or lost", pch=18, main=tl, ylim=ylm, xlim=c(0, max(cumsum(chrominfo$length))))
        
        points(kb.loc[gl$lossP>=cutplot],-gl$lossP[gl$lossP>=cutplot], col="red", type="h")
        
        abline(h=0)
        abline(h=seq(-.8,.8,b=.2), lty=2,lwd=.5)
        abline(v=cumsum(chrominfo$length), col="blue")
        abline(v=chrom.centr, lty=2, col="grey50")
        for (i in seq(1,(numaut),b=2))
        {
            mtext(paste("", i), side = 1, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
        }
        for (i in seq(2,(numaut),b=2))
        {
            mtext(paste("", i), side = 3, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
        }
        
        if(X)
        {
            if (i == numaut)
            {
                mtext("X", side = 1, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
            }
            else
            {
                mtext("X", side = 3, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
            }
        }
        if (Y)
        {
            if (i == numaut)
            {
                mtext("Y", side = 3, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
            }
            else
            {
                mtext("Y", side = 1, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
            }
            
        }
        
        
    }
    if (sign)
    {
        plot(kb.loc,teststat, col="blue", ylim=c(0,max(teststat)), type="h", xlab="chromosome number", ylab="clone statistic", pch=18, main=paste(titles, collapse=" vs "))
        if (length(st.now) > 0)
        {
            abline(h=rev(st.now), col=rev(pal.now), lty=2)
        }
        abline(v=cumsum(chrominfo$length), col="black")
        abline(v=chrom.centr, lty=2, col="grey50")
        for (i in seq(1,(numaut),b=2))
        {
            mtext(paste("", i), side = 1, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
        }
        for (i in seq(2,(numaut),b=2))
        {
            mtext(paste("", i), side = 3, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
        }
        
        if(X)
        {
            if (i == numaut)
            {
                mtext("X", side = 1, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
            }
            else
            {
                mtext("X", side = 3, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
            }
        }
        if (Y)
        {
            if (i == numaut)
            {
                mtext("Y", side = 3, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
            }
            else
            {
                mtext("Y", side = 1, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
            }
            
        }
        
        
    }

    dev.off()

}


#############################################################################
#############################################################################
##frequency plot for individual chromosomes

##data
##rsp -- phenotype, NA are not allowed, have to be consequetive integers
##chrom -- vector of chromosomes to show. currently bugged; so has to be 1:smthg
##datainfo
##chrominfo
##titles -- the titles of the frequency plots -- has to have as many names as 
##levels in the response
##thres -- unique threshold or vector of tumor specific thresholds. In the latter
##case has to contain as many thershold as samples
##sign -- to do significance comparison (T) or not (F). to do comparison uses
##multtest package
##nperm -- if sign =TRUE, then how many permutations for maxT
##test -- name of the test
##ranks -- whether to work with ranked data If nor, "n"
##side -- two sisded (abs) test or 1-sided ("upper" or "lower")
##p.thres -- for which adjusted p-values show the cut-off
##filePS -- name of the ps/pdf file
##PS = "ps" or "pdf"

plotfreq.stat.chrom.final.func <-
    function(data, rsp, chrom, datainfo, chrominfo, titles, thres=.2,
             ylm=c(-1,1), sign=TRUE, nperm=1000, test="f", ranks="y",
             side="abs", p.thres=c(.01,.05,.1,.2), dataSign=data,
             filePS="gainslosses.ps", PS="ps", onepage=TRUE)
{
##########compute gainloss functions:
    rsp.uniq <- unique(rsp)
    colmatr <- matrix(0, nrow=length(rsp.uniq), ncol=length(rsp))
    for (i in 1:nrow(colmatr))
    {
	colmatr[i,rsp==rsp.uniq[i]] <- 1
    }



    pal <- c("red", "blue", "green", "yellow")

    if (nrow(colmatr) == 1)
    {
	sign <- F
    }

    nr <- nrow(colmatr)
    if (sign)
    {	
	nr <- nr+1
    }

    tmp <- matrix(0, ncol=2,nrow=1)   
    tmp <- as.data.frame(tmp) 
    dimnames(tmp)[[2]] <- c("gainP", "lossP")   
    gainloss <- rep(list(tmp),nrow(colmatr))            

    for (j in 1:nrow(colmatr))
    {
	
	cols <- (1:ncol(colmatr))[colmatr[j,]==1]
	gainloss[[j]] <- gainLoss(dat=data, cols=cols,thres=thres)
	
    }

    if (sign)
    {
	dt <- dataSign[,colmatr[1,]==1]
	rsp <- rep(1, ncol(dt))
	for (j in 2:nrow(colmatr))
	{
            dt <- cbind(dt, dataSign[,colmatr[j,]==1])
            rsp <- c(rsp, rep(j, ncol(dataSign[,colmatr[j,]==1])))
	}
	rsp <- rsp-1
	res <-  mt.maxT(X=dt,classlabel=rsp,test=test,side=side,fixed.seed.sampling="y",B=nperm, na=.mt.naNUM, nonpara=ranks)
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
        
        ##st.now <- st[!is.na(st)]
        ##pal.now <- pal[!is.na(st)]

	st.now <- st
	pal.now <- pal
	
	

    }
    
    if (PS == "ps")
    {
        postscript(filePS,paper="letter")
    }
    else if (PS == "pdf")
    {
        pdf(filePS, width = 8.5, height = 11)
    }
    else
    {
        print("no legimate file format is specified")
        exit()
    }
    nc <- 1
    ##	for (ch in 1:length(chrom))
    for (ch in chrom)
    {
        ##chrom.ind <- (1:nrow(datainfo))[datainfo$Chrom== chrom[ch]]
        chrom.ind <- (1:nrow(datainfo))[datainfo$Chrom== ch]
        
        kb <- datainfo$kb[chrom.ind]
        chrom.centr <- chrominfo$centr[ch]
        if (onepage)
        {
            par(mfrow=c(nr,nc), xaxs="i")
        }
        else
        {
            par(mfrow=c(1,nc), xaxs="i")
        }
        
        for (g in 1:length(titles))
        {
            gl <- gainloss[[g]]
            tl <- titles[g]
            plot(kb,(gl$gainP[chrom.ind]), col="green", ylim=ylm, type="h", xlab="chromosome number", ylab="Fraction gained or lost", pch=18, main=paste(tl, " chrom ", ch))
            
            points(kb,(-gl$lossP[chrom.ind]), col="red", type="h")
            
            abline(h=0)
            abline(h=seq(-.8,.8,b=.2), lty=2,lwd=.5)
            abline(v=chrom.centr, lty=2, col="grey50")
            
            
        }
	if (sign)
	{
            plot(kb,teststat[chrom.ind], col="blue", ylim=c(0,max(teststat[chrom.ind])), type="h", xlab="kb", ylab="clone statistic", pch=18, main=paste(paste(titles, collapse=" vs "), " -- chrom ", ch))
            if (length(st.now) > 0)
            {
                
                abline(h=rev(st.now), col=rev(pal.now), lty=2)
            }
            
            abline(v=chrom.centr, lty=2, col="grey50")
            
            
	}
        
    }
    dev.off()

}


#############################################################################
#############################################################################
##################################################################################
##################################################################################

##frequency plot for the whole genome

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
##nperm -- if sign =TRUE, then how many permutations for maxT
##test -- name of the test
##ranks -- whether to work with ranked data If nor, "n"
##side -- two sisded (abs) test or 1-sided ("upper" or "lower")
##p.thres -- for which adjusted p-values show the cut-off
##filePS -- name of the ps/pdf file
##PS = "ps" or "pdf"
##X=TRUE -- is X (23) chrom included?
##Y=FALSE  -- is Y (24) chrom included?
##numaut=22 -- number of autosomes
##ngrid=50: density of colors
##nlim=1: upper limit for solors
##mincolors: minimum limit for colors
##by defult "white"corersponds to [-.2,.2] and red and green to [-1,-.2] and [.2,1[]
##respectively
## quant.col=.5: percentile for coloring gaind/lost clones -- <= .5, E.g. .25
##would correspond to taking 25th % for lost and 75% for gained samples


##plotfreq.stat.final.colors.func <- function(data, rsp,datainfo, chrominfo, titles, thres=.2,cutplot = 0, ylm=c(-1,1), sign=FALSE, dataSign=data, nperm=1000, test="f", ranks="y", side="abs", p.thres=c(.01,.05,.1,.2), filePS="gainslosses.ps", ngrid=50, nlim=1, mincolors=.2, quant.col=.5, X=TRUE, Y=FALSE, numaut=22, PS="ps", onepage=TRUE)
##{

##rsp.uniq <- unique(rsp)
##colmatr <- matrix(0, nrow=length(rsp.uniq), ncol=length(rsp))
##for (i in 1:nrow(colmatr))
##{
##	colmatr[i,rsp==rsp.uniq[i]] <- 1
##}



##pal <- c("red", "blue", "green", "yellow")
##pal <- pal[1:length(p.thres)]

##colors.gain <- maPalette(low = "white", high = "green", k=ngrid)
##colors.loss <- maPalette(low = "red", high = "white", k=ngrid)

##sq.loss <- seq(-nlim, -mincolors, length=(ngrid+1))
##sq.gain <- seq(mincolors, nlim, length=(ngrid+1))

######################################################

##matr.colors.loss <- data.frame(sq.loss[-length(sq.loss)], sq.loss[-1], colors.loss)
##matr.colors.gain <- data.frame(sq.gain[-length(sq.gain)], sq.gain[-1], colors.gain)

##if (nrow(colmatr) == 1)
##{
##	sign <- F
##}

##nr <- nrow(colmatr)
##if (sign)
##{	
##	nr <- nr+1

##}

##tmp <- matrix(0, ncol=2,nrow=1)   
##tmp <- as.data.frame(tmp) 
##dimnames(tmp)[[2]] <- c("gainP", "lossP")   
##gainloss <- rep(list(tmp),nrow(colmatr))            

##for (j in 1:nrow(colmatr))
##{

##	cols <- (1:ncol(colmatr))[colmatr[j,]==1]
##	gainloss[[j]] <- gainLoss(dat=data, cols=cols,thres=thres, quant=quant.col)

##}



##if (sign)
##{
##	dt <- dataSign[,colmatr[1,]==1]
##	rsp <- rep(1, ncol(dt))
##	for (j in 2:nrow(colmatr))
##	{
##		dt <- cbind(dt, dataSign[,colmatr[j,]==1])
##		rsp <- c(rsp, rep(j, ncol(dataSign[,colmatr[j,]==1])))
##	}
##	rsp <- rsp-1
##	res <-  mt.maxT(X=dt,classlabel=rsp,test=test,side=side,fixed.seed.sampling="y",B=nperm, na=.mt.naNUM, nonpara=ranks)
##	maxT <- res$adjp[order(res$index)]
##	#rawp <- res$rawp[order(res$index)]
##	teststat <- abs(res$teststat[order(res$index)])
##	st <- rep(NA, length(p.thres))
##	for (s in 1:length(p.thres))
##	{
##		if (length(maxT[maxT<=p.thres[s]]) > 0)
##		{
##			st[s] <- min(teststat[maxT<=p.thres[s]])
##		}
##	}

##	#st.now <- st[!is.na(st)]
##	#pal.now <- pal[!is.na(st)]

##	st.now <- st
##	pal.now <- pal

##}



##numchr <- numaut
##if (X)
##{
##	numchr <- numchr+1
##}
##if (Y)
##{
##	numchr <- numchr+1
##}

##chrominfo <- chrominfo[1:numchr,]


###compute cumulative kb locations
##        start <- c(0, cumsum(chrominfo$length))
##        kb.loc <- datainfo$kb
##        for (i in 1:nrow(chrominfo))
##        {
##                tmp <- start[i]+datainfo$kb[datainfo$Chrom==i]
##                kb.loc[datainfo$Chrom==i] <- tmp
##        }
###preparation for graphs
##        chrom.start <- rep(0,nrow(chrominfo))
##        for (i in 2:length(chrom.start))
##        {
##                chrom.start[i] <- sum(chrominfo$length[1:(i-1)])

##        }
##        chrom.centr <- rep(0,nrow(chrominfo))
##        for (i in 1:length(chrom.centr))
##        {
##                chrom.centr[i] <- chrom.start[i]+chrominfo$centr[i]

##        }

##        chrom.mid <- rep(0, nrow(chrominfo))
##        for (i in 1:length(chrom.start))
##        {
##              chrom.mid[i] <- chrom.start[i]+chrominfo$length[i]/2
##        }

###now, plot
##        #nc <- max(length(titles)/2,1)
##	nc <- 1

##	if (PS == "ps")
##	{
##		postscript(filePS,paper="letter")
##	}
##	else if (PS == "pdf")
##	{
##		pdf(filePS, width = 8.5, height =11)
##      	}
##	if (onepage)
##	{
##		par(mfrow=c(nr,nc), lab=c(1,8,7), tcl=-.2,  xaxs="i")
##	}
##	else
##	{
##		par(mfrow=c(1,nc), lab=c(1,8,7), tcl=-.2,  xaxs="i")
##	}
##        for (g in 1:length(titles))
##        {
##                gl <- gainloss[[g]]
##		tl <- titles[g]
##		ylm[1] <- min(ylm, min(gl$lossP))
##		ylm[2] <- max(ylm, max(gl$gainP))

##		cl <- gl$gainMed	

##		col.nrow <- rep(0, length(cl))
##		for (i in 1:length(cl))
##		{
##			if (cl[i]>=nlim)
##			{
##				cl[i] <- nlim-10^(-6)
##			}
##			if (length((1:nrow(matr.colors.gain))[cl[i]>=matr.colors.gain[,1] & cl[i]<matr.colors.gain[,2]]) > 0)
##			{
##				col.nrow[i] <- (1:nrow(matr.colors.gain))[cl[i]>=matr.colors.gain[,1] & cl[i]<matr.colors.gain[,2]]
##			}
##			else
##			{
##				col.nrow[i] <- 1
##			}
##		}		


##		plot(kb.loc[gl$gainP>=cutplot],gl$gainP[gl$gainP>=cutplot], col=as.character(matr.colors.gain[gl$gainP>=cutplot,3][col.nrow[gl$gainP>=cutplot]]), type="h", xlab="chromosome number", ylab="Fraction gained or lost", pch=18, main=tl, ylim=ylm, xlim=c(0, max(cumsum(chrominfo$length))))

##		cl <- gl$lossMed

##		col.nrow <- rep(0, length(cl))
##		for (i in 1:length(cl))
##		{
##			if (cl[i]<=-nlim)
##			{
##				cl[i] <- -nlim+10^(-6)
##			}
##			if (length((1:nrow(matr.colors.loss))[cl[i]>=matr.colors.loss[,1] & cl[i]<matr.colors.loss[,2]]) > 0)
##			{
##				col.nrow[i] <- (1:nrow(matr.colors.loss))[cl[i]>=matr.colors.loss[,1] & cl[i]<matr.colors.loss[,2]]
##			}
##			else
##			{
##				col.nrow[i] <- ngrid
##			}
##		}		


##                points(kb.loc[gl$lossP>=cutplot],-gl$lossP[gl$lossP>=cutplot], col=as.character(matr.colors.loss[gl$lossP>=cutplot,3][col.nrow[gl$lossP>=cutplot]]), type="h")

##                abline(h=0)
##                abline(h=seq(-.8,.8,b=.2), lty=2,lwd=.5)
##                abline(v=cumsum(chrominfo$length), col="blue")
##                abline(v=chrom.centr, lty=2, col="grey50")

##		for (i in seq(2,(numaut),b=2))
##                {
##                     mtext(paste("", i), side = 3, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
##		}
##		for (i in seq(1,(numaut),b=2))
##                {
##                     mtext(paste("", i), side = 1, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
##		}

##		if(X)
##		{
##			if (i == numaut)
##			{
##				mtext("X", side = 1, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
##			}
##			else
##			{
##				mtext("X", side = 3, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
##			}
##		}
##		if (Y)
##		{
##			if (i == numaut)
##			{
##				mtext("Y", side = 3, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
##			}
##			else
##			{
##				mtext("Y", side = 1, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
##			}

##		}

##        }
##	if (sign)
##	{
##		plot(kb.loc,teststat, col="blue", ylim=c(0,max(teststat)), type="h", xlab="chromosome number", ylab="clone statistic", pch=18, main=paste(titles, collapse=" vs "), xlim=c(0, max(cumsum(chrominfo$length))))
##		if (length(st.now) > 0)
##		{
##			abline(h=rev(st.now), col=rev(pal.now), lty=2)
##		}
##		abline(v=cumsum(chrominfo$length), col="black")
##                abline(v=chrom.centr, lty=2, col="grey50")

##		for (i in seq(1,(numaut),b=2))
##                {
##                     mtext(paste("", i), side = 1, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
##		}
##		for (i in seq(2,(numaut),b=2))
##                {
##                     mtext(paste("", i), side = 3, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
##		}

##		if(X)
##		{
##			if (i == numaut)
##			{
##				mtext("X", side = 1, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
##			}
##			else
##			{
##				mtext("X", side = 3, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
##			}
##		}
##		if (Y)
##		{
##			if (i == numaut)
##			{
##				mtext("Y", side = 3, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
##			}
##			else
##			{
##				mtext("Y", side = 1, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
##			}

##		}


##	}

##	dev.off()

##}


#############################################################################
#############################################################################
#############################################################################
#############################################################################
##################################################################################
##################################################################################

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
##nperm -- if sign =TRUE, then how many permutations for maxT
##test -- name of the test
##ranks -- whether to work with ranked data If nor, "n"
##side -- two sisded (abs) test or 1-sided ("upper" or "lower")
##p.thres -- for which adjusted p-values show the cut-off
##filePS -- name of the ps/pdf file
##PS = "ps" or "pdf"
##X=TRUE -- is X (23) chrom included?
##Y=FALSE  -- is Y (24) chrom included?
##numaut=22 -- number of autosomes
##ngrid=50: density of colors
##nlim=1: upper limit for solors
##mincolors: minimum limit for colors
##by defult "white"corersponds to [-.2,.2] and red and green to [-1,-.2] and [.2,1[]
##respectively
## quant.col=.5: percentile for coloring gaind/lost clones -- <= .5, E.g. .25
##would correspond to taking 25th % for lost and 75% for gained samples

plotfreq.givenstat.final.colors.func <-
    function(data, rsp,datainfo, chrominfo, titles, thres=.2,cutplot =
             0, ylm=c(-1,1), sign=FALSE, stat=stats, statPerm=statsPerm,
             p.thres=c(.01,.05,.1,.2), filePS="gainslosses.ps",
             ngrid=50, nlim=1, mincolors=.2, quant.col=.5, X=TRUE, Y=FALSE,
             numaut=22, PS="ps", onepage=TRUE)
{

    rsp.uniq <- unique(rsp)
    colmatr <- matrix(0, nrow=length(rsp.uniq), ncol=length(rsp))
    for (i in 1:nrow(colmatr))
    {
	colmatr[i,rsp==rsp.uniq[i]] <- 1
    }



    pal <- c("red", "blue", "green", "yellow")
    pal <- pal[1:length(p.thres)]

    colors.gain <- maPalette(low = "white", high = "green", k=ngrid)
    colors.loss <- maPalette(low = "red", high = "white", k=ngrid)

    sq.loss <- seq(-nlim, -mincolors, length=(ngrid+1))
    sq.gain <- seq(mincolors, nlim, length=(ngrid+1))
    
####################################################

    matr.colors.loss <- data.frame(sq.loss[-length(sq.loss)], sq.loss[-1], colors.loss)
    matr.colors.gain <- data.frame(sq.gain[-length(sq.gain)], sq.gain[-1], colors.gain)

    if (nrow(colmatr) == 1)
    {
	sign <- F
    }

    nr <- nrow(colmatr)
    if (sign)
    {	
	nr <- nr+1
	
    }

    tmp <- matrix(0, ncol=2,nrow=1)   
    tmp <- as.data.frame(tmp) 
    dimnames(tmp)[[2]] <- c("gainP", "lossP")   
    gainloss <- rep(list(tmp),nrow(colmatr))            

    for (j in 1:nrow(colmatr))
    {
	
	cols <- (1:ncol(colmatr))[colmatr[j,]==1]
	gainloss[[j]] <- gainLoss(dat=data, cols=cols,thres=thres, quant=quant.col)
	
    }



    if (sign)
    {
	
	teststat <- stat
	maxstat <- apply(statPerm,2,max, na.rm=TRUE)
	maxT <- rep(NA, length(maxstat))
	for (i in 1:length(teststat)) 
	{
            maxT[i] <- length(maxstat[maxstat>=teststat[i]])
	}
	maxT <- maxT/length(maxstat)
	
	st <- quantile(maxstat, (1-p.thres))
	st.now <- st
	pal.now <- pal
	
    }



    numchr <- numaut
    if (X)
    {
	numchr <- numchr+1
    }
    if (Y)
    {
	numchr <- numchr+1
    }

    chrominfo <- chrominfo[1:numchr,]


    ##compute cumulative kb locations
    start <- c(0, cumsum(chrominfo$length))
    kb.loc <- datainfo$kb
    for (i in 1:nrow(chrominfo))
    {
        tmp <- start[i]+datainfo$kb[datainfo$Chrom==i]
        kb.loc[datainfo$Chrom==i] <- tmp
    }
    ##preparation for graphs
    chrom.start <- rep(0,nrow(chrominfo))
    for (i in 2:length(chrom.start))
    {
        chrom.start[i] <- sum(chrominfo$length[1:(i-1)])

    }
    chrom.centr <- rep(0,nrow(chrominfo))
    for (i in 1:length(chrom.centr))
    {
        chrom.centr[i] <- chrom.start[i]+chrominfo$centr[i]

    }

    chrom.mid <- rep(0, nrow(chrominfo))
    for (i in 1:length(chrom.start))
    {
        chrom.mid[i] <- chrom.start[i]+chrominfo$length[i]/2
    }

    ##now, plot
    ##nc <- max(length(titles)/2,1)
    nc <- 1
    
    if (PS == "ps")
    {
        postscript(filePS,paper="letter")
    }
    else if (PS == "pdf")
    {
        pdf(filePS, width = 8.5, height =11)
    }
    
    if (onepage)
    {
        par(mfrow=c(nr,nc), lab=c(1,8,7), tcl=-.2,  xaxs="i")
    }
    else
    {
        par(mfrow=c(1,nc), lab=c(1,8,7), tcl=-.2,  xaxs="i")
    }
    for (g in 1:length(titles))
    {
        gl <- gainloss[[g]]
        tl <- titles[g]
        ylm[1] <- min(ylm, min(gl$lossP))
        ylm[2] <- max(ylm, max(gl$gainP))
        
        cl <- gl$gainMed	
        
        col.nrow <- rep(0, length(cl))
        for (i in 1:length(cl))
        {
            if (cl[i]>=nlim)
            {
                cl[i] <- nlim-10^(-6)
            }
            if (length((1:nrow(matr.colors.gain))[cl[i]>=matr.colors.gain[,1] & cl[i]<matr.colors.gain[,2]]) > 0)
            {
                col.nrow[i] <- (1:nrow(matr.colors.gain))[cl[i]>=matr.colors.gain[,1] & cl[i]<matr.colors.gain[,2]]
            }
            else
            {
                col.nrow[i] <- 1
            }
        }		

        
        plot(kb.loc[gl$gainP>=cutplot],gl$gainP[gl$gainP>=cutplot],
             col=as.character(matr.colors.gain[gl$gainP>=cutplot,3][col.nrow[gl$gainP>=cutplot]]),
             type="h", xlab="chromosome number",
             ylab="Fraction gained or lost", pch=18, main=tl,
             ylim=ylm, xlim=c(0, max(cumsum(chrominfo$length))))
        
        
        cl <- gl$lossMed
        
        col.nrow <- rep(0, length(cl))
        for (i in 1:length(cl))
        {
            if (cl[i]<=-nlim)
            {
                cl[i] <- -nlim+10^(-6)
            }
            if (length((1:nrow(matr.colors.loss))[cl[i]>=matr.colors.loss[,1] & cl[i]<matr.colors.loss[,2]]) > 0)
            {
                col.nrow[i] <-
                    which(cl[i]>=matr.colors.loss[,1] &
                          cl[i]<matr.colors.loss[,2])
            }
            else
            {
                col.nrow[i] <- ngrid
            }
        }		

        
        points(kb.loc[gl$lossP>=cutplot],-gl$lossP[gl$lossP>=cutplot],
               col=as.character(matr.colors.loss[gl$lossP>=cutplot,3][col.nrow[gl$lossP>=cutplot]]),
               type="h")
        
        abline(h=0)
        abline(h=seq(-.8,.8,b=.2), lty=2,lwd=.5)
        abline(v=cumsum(chrominfo$length), col="blue")
        abline(v=chrom.centr, lty=2, col="grey50")

        for (i in seq(2,(numaut),b=2))
        {
            mtext(paste("", i), side = 3, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
        }
        for (i in seq(1,(numaut),b=2))
        {
            mtext(paste("", i), side = 1, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
        }
        
        if(X)
        {
            if (i == numaut)
            {
                mtext("X", side = 1, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
            }
            else
            {
                mtext("X", side = 3, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
            }
        }
        if (Y)
        {
            if (i == numaut)
            {
                mtext("Y", side = 3, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
            }
            else
            {
                mtext("Y", side = 1, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
            }
            
        }
        
    }
    if (sign)
    {
        plot(kb.loc,teststat, col="blue", ylim=c(0,max(teststat)), type="h", xlab="chromosome number", ylab="clone statistic", pch=18, main=paste(titles, collapse=" vs "), xlim=c(0, max(cumsum(chrominfo$length))))
        if (length(st.now) > 0)
        {
            abline(h=rev(st.now), col=rev(pal.now), lty=2)
        }
        abline(v=cumsum(chrominfo$length), col="black")
        abline(v=chrom.centr, lty=2, col="grey50")

        for (i in seq(1,(numaut),b=2))
        {
            mtext(paste("", i), side = 1, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
        }
        for (i in seq(2,(numaut),b=2))
        {
            mtext(paste("", i), side = 3, at = (chrom.mid[i]), line=.3, col="red", cex.main=.5)
        }
        
        if(X)
        {
            if (i == numaut)
            {
                mtext("X", side = 1, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
            }
            else
            {
                mtext("X", side = 3, at = (chrom.mid[numaut+1]), line=.3, col="red", cex.main=.5)
            }
        }
        if (Y)
        {
            if (i == numaut)
            {
                mtext("Y", side = 3, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
            }
            else
            {
                mtext("Y", side = 1, at = (chrom.mid[numaut+2]), line=.3, col="red", cex.main=.5)
            }
            
        }

        
    }

    dev.off()

}
