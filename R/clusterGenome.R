clusterGenome <-
    function(aCGH.obj, response = as.factor(rep("All", ncol(aCGH.obj))),
             chrominfo = human.chrom.info.Jul03, cutoff=1,
             lowCol = "red", highCol = "green", midCol = "black",
             ncolors = 50, byclass = FALSE, showaber = FALSE, amplif = 1,
             homdel = -0.75, samplenames = sample.names(aCGH.obj),
             vecchrom = 1:23, titles = "Image Plot", methodS = "ward",
             dendPlot = TRUE, imp = TRUE, categoricalPheno = TRUE)
{
	aCGH.obj <- aCGH.obj[,!is.na(response)]
	samplenames <- samplenames[!is.na(response)]
	response <- response[!is.na(response)]
    if (categoricalPheno)
    {
        resp0 <- response
        resp0.num <- as.numeric(as.factor(resp0))
        resp <- as.numeric(as.factor(resp0))
        if (!(byclass))
        {
            resp <- rep(1, length(resp0))
        }

        tbl.resp <- table(resp)
        label.col <- rainbow(length(unique(resp)))
    } 
    else
    {
	byclass <- FALSE
	resp0 <- response
    	resp0.num <- resp0 
    	
	resp <- rep(1, length(resp0))
    

    	tbl.resp <- table(resp)
    ##label.col <- c("red", "green", "blue", "skyblue", "orange", "pink", "gray20")
    	label.col <- rainbow(length(unique(resp)))
    }

    #par(bg="grey20")
    
    datainfo <- clones.info(aCGH.obj)
    if (imp)
    {
    	data <- log2.ratios.imputed(aCGH.obj)
    }
    else
    {
	data <- log2.ratios(aCGH.obj)
    }
    indUse <- NA
    chromb <- 0
    for (i in 1:length(vecchrom))
    {
	indUse <- c(indUse, which(datainfo$Chrom == vecchrom[i]))
        chromb <- c(chromb, length(which(datainfo$Chrom == vecchrom[i])))
        
    }
    indUse <- indUse[-1]
    chromb <- cumsum(chromb)   

    datainfo <- datainfo[indUse,]
    data <- data[indUse,]	
    kb <- datainfo$kb
    data <- as.matrix(data)

    if (dendPlot)
    {
	#ind.pres <- which(!is.na(response))
	cr <- dist(t( data))
        hcl <- hclust(cr, method=methodS)
    }

    dt.cp <- data
    dt <- apply(data, 2,floor.func, cutoff)    
    

    dt <- dt[,order(resp)]
    dt.cp <- dt.cp[,order(resp)]

    resp0 <- resp0[order(resp)]
    resp0.num <- resp0.num[order(resp)]

    samplenames <- samplenames[order(resp)]
    resp <- resp[order(resp)]


    start <- 1
    ##mapping order
    ord <- rep(0, length(resp))
    for (i in 1:(length(tbl.resp)))
    {
	
	ind <- which(resp == i)
	#cr <- as.dist(1-cor.na(data[,ind]))
        cr <- dist(t(data[,ind]))
	ord[start:sum(tbl.resp[1:i])] <- hclust(cr, method=methodS)$ord+start-1
	start <- sum(tbl.resp[1:i])+1
	
	
	
    }
    dt <- dt[,ord]
    dt.cp <- dt.cp[,ord]

    resp <- resp[ord]
    resp0 <- resp0[ord]
    resp0.num <- resp0.num[ord]
    samplenames <- samplenames[ord]
    image(x=(1:length(kb)), y=1:length(resp), z=dt, col =
          maPalette(low = lowCol, high = highCol, mid = midCol
                    ,k=ncolors), axes = FALSE, xlab = "", ylab = "",
          zlim=c(-cutoff,cutoff))
    
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
	mtext((resp0)[j], side = 2, at = j, line=.3, col=label.col[((resp0.num)[j])], cex=.5, las=2)
	mtext(paste((samplenames)[j], ""), side = 4, at = j, line=.3, col=label.col[((resp0.num)[j])], cex=.3, las=2)
	
    }
    ##title(main="Whole genome", xlab = "clone", ylab = "sample", col.lab="white", col.main="white")
    title(xlab = "clone", ylab = "sample", main=titles)
    ##abline(v=centrloc, col="white", lty=2, lwd=.5)
    
    abline(v=chromb, col="white", lty=1, lwd=1)
    loc <- chromb[-1]-diff(chromb)/2
    if (length(vecchrom) > 1)
    {
    for (i in seq(2,length(vecchrom),b=2))
    {
        
        mtext(paste("", vecchrom[i]), side = 3, at = (loc[i]), line=.3,cex.main=.5)
    }
    }
    for (i in seq(1,length(vecchrom),b=2))
    {
        
        mtext(paste("", vecchrom[i]), side = 1, at = (loc[i]), line=.3, cex.main=.5)
    }

    ##mtext("X", side = 1, at = (loc[nrow(chrominfo)]), line=.3,col="white", cex.main=.5)


  if (dendPlot)
  {
	if (length(unique(resp0)) > 1)
	{
	   plot(hcl, labels=response, main="Dendogram")
	}
	else
	{
	  plot(hcl, labels=(sample.names(aCGH.obj)), main="Dendogram")
	}		
  }
}
