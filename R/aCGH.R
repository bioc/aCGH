
create.aCGH <-
    function(log2.ratios, clones.info, phenotype = NULL)
{
    
    if (nrow(log2.ratios) != nrow(clones.info))
        stop("Number of rows of log2.ratios and clones.info differ!")
    if (!is.null(phenotype) && ncol(log2.ratios) != nrow(phenotype))
        stop("Number of columns of log2.ratios and number of rows in\
phenotype differ!")
##    if (!all(rownames(log2.ratios) == clones.info$Clone))
##        rownames(log2.ratios) <- clones.info$Clone
    value <-
        list(log2.ratios = log2.ratios,
             clones.info = clones.info,
             phenotype = phenotype)
    class(value) <- "aCGH"
    attr(value, "call") <- sys.call()
    value
    
}

log2.ratios <- function(aCGH.obj) aCGH.obj$log2.ratios
##"log2.ratios<-" <-
##    function(aCGH.obj, value)
##{

##    if (!is.aCGH(aCGH.obj))
##	stop("object is not of class aCGH")
##    if (any(dim(value) != dim(aCGH.obj$log2.ratios)))
##        stop("invalid replacement dimensions")
##    aCGH.obj$log2.ratios <- value
##    aCGH.obj

##}

clones.info <- function(aCGH.obj) aCGH.obj$clones.info
##"clones.info<-" <-
##    function(aCGH.obj, value)
##{

##    if (!is.aCGH(aCGH.obj))
##	stop("object is not of class aCGH")
##    if (any(dim(value) != dim(aCGH.obj$clones.info)))
##        stop("invalid replacement dimensions")
##    aCGH.obj$clones.info <- value
##    aCGH.obj

##}

is.aCGH <- function(aCGH.obj) inherits(aCGH.obj, "aCGH")

dim.aCGH <- function(aCGH.obj) dim(aCGH.obj$log2.ratios)

num.clones <- nrow.aCGH <-
    function(aCGH.obj) nrow(aCGH.obj$log2.ratios)

num.samples <- ncol.aCGH <-
    function(aCGH.obj) ncol(aCGH.obj$log2.ratios)

num.chromosomes <- function(aCGH.obj) length(unique(aCGH.obj$clones.info$Chrom))

clone.names <- row.names.aCGH <- rownames.aCGH <-
    function(x) x$clones.info$Clone
"clone.names<-" <- "row.names<-.aCGH" <- "rownames<-.aCGH" <-
    function(x, value)
{
    
    if (!is.aCGH(x))
	stop("object is not of class aCGH")
    if (length(value) != length(x$clones.info$Clone))
        stop("invalid replacement dimensions")
    row.names(x$clones.info$Clone) <- as.factor(value)
    x
    
}

colnames.aCGH <- col.names.aCGH <- sample.names <-
    function(aCGH.obj) colnames(aCGH.obj$log2.ratios)
"colnames<-.aCGH" <- "col.names<-.aCGH" <- "sample.names<-" <-
    function(aCGH.obj, value)
{
    
    if (!is.aCGH(aCGH.obj))
	stop("object is not of class aCGH")
    if (length(value) != ncol(aCGH.obj$log2.ratios))
        stop("invalid replacement dimensions")
    colnames(aCGH.obj$log2.ratios) <- value
    aCGH.obj
    
}

impute.lowess <-
    function(aCGH.obj, chrominfo = human.chrom.info.Jul03,
             maxChrom = 23, smooth = 0.1)
{

    data.imp <- log2.ratios <- log2.ratios(aCGH.obj)
    clones.info <- clones.info(aCGH.obj)
    uniq.chrom <- unique(clones.info$Chrom)
    for (j in uniq.chrom[uniq.chrom <= maxChrom])
    {
        
        cat("Processing chromosome ", j, "\n")
        centr <- chrominfo$centromere[j]
        indl <-
            which(clones.info$Chrom == j & clones.info$kb < centr)
        indr <-
            which(clones.info$Chrom == j & clones.info$kb > centr)
        kbl <- clones.info$kb[indl]
        kbr <- clones.info$kb[indr]
	
        for (i in 1:ncol(log2.ratios))
        {
            
            ##print(i)
            if (length(indl) > 0)
            {
                
                vecl <- log2.ratios[indl, i]
                ind <- which(!is.na(vecl))
                if (length(ind) > 0)
                    data.imp[indl, i][-ind] <-
                        approx(lowess(kbl[ind], vecl[ind], f = smooth),
                               xout = kbl[-ind])$y
                
            }
            if (length(indr) > 0)
            {
                
                vecr <- log2.ratios[indr, i]
                ind <- which(!is.na(vecr))
                if (length(ind) > 0)
                    data.imp[indr, i][-ind] <-
                        approx(lowess(kbr[ind], vecr[ind], f = smooth),
                               xout = kbr[-ind])$y
                
            }
            
        }
        
    }

#################
###now, if any missing values are left 
    
    prop.miss <- apply(data.imp, 2, prop.na)
    ## if any samples contain missing values
    if (max(prop.miss, na.rm = TRUE) > 0)
    {
        
        for (i in 1:ncol(data.imp))
        {
            
            vec <- data.imp[ ,i ]
            ind <- which(is.na(vec))
            if (length(ind) > 0)
            {
                
                vec[ind] <-
                    sapply(ind,
                           function(i) {

                               chr <- clones.info$Chrom[i]
                               kb <- clones.info$kb[i]
                               if (kb >= chrominfo$centromere[chr])
                                   median(vec[clones.info$Chrom == chr
                                              & clones.info$kb >=
                                              chrominfo$centromere[chr]],
                                          na.rm = TRUE)
                               else
                                   median(vec[clones.info$Chrom == chr
                                              & clones.info$kb <
                                              chrominfo$centromere[chr]],
                                          na.rm = TRUE)
                               
                           }
                           )

                ##if all values on chrom are missing
                vec[is.na(vec)] <- 0
                data.imp[,i] <- vec
            
            }
            
        }
        
    }
    prop.miss <- apply(data.imp, 2, prop.na)
    if (max(prop.miss) > 0)
        print(paste("Missing values still remain in samples ",
                    which(prop.miss > 0)))
    
    data.imp
    
}

log2.ratios.imputed <-
    function(aCGH.obj)
    aCGH.obj$log2.ratios.imputed

"log2.ratios.imputed<-" <-
    function(aCGH.obj, value)
{

    if (!is.aCGH(aCGH.obj))
	stop("object is not of class aCGH")
    if (!is.null(aCGH.obj$log2.ratios.imputed) &&
        any(dim(value) != dim(aCGH.obj$log2.ratios.imputed)))
        stop("invalid replacement dimensions")
    aCGH.obj$log2.ratios.imputed <- value
    aCGH.obj

}

find.hmm.states <-
    function(aCGH.obj, ...)
    hmm.run.func(aCGH.obj$log2.ratios, aCGH.obj$clones.info, ...)

hmm <- function(aCGH.obj) aCGH.obj$hmm
"hmm<-" <-
    function(aCGH.obj, value)
{

    if (!is.aCGH(aCGH.obj))
	stop("object is not of class aCGH")
    if (!is.null(aCGH.obj$hmm))
    {
        
        nstates.ok <-
            length(value$nstates.hmm) ==
                length(aCGH.obj$hmm$nstates.hmm) &&
        all(sapply(1:length(aCGH.obj$hmm$nstates.hmm),
                   function(i)
                   all(dim(aCGH.obj$hmm$nstates.hmm[[i]]) ==
                       dim(value$nstates.hmm[[i]]))
                   ))
        states.ok <-
            length(value$states.hmm) ==
                length(aCGH.obj$hmm$states.hmm) &&
        all(sapply(1:length(aCGH.obj$hmm$states.hmm),
                   function(i)
                   all(dim(aCGH.obj$hmm$states.hmm[[i]]) ==
                       dim(value$states.hmm[[i]]))
                   ))
        if (!nstates.ok || !states.ok)
            stop("invalid replacement dimensions")
        
    }
    aCGH.obj$hmm <- value
    aCGH.obj

}

mergeHmmStates <-
    function(aCGH.obj, model.use = 1, minDiff = .25)
{

    if (is.null(hmm(aCGH.obj)))
	stop("compute the hmm states first using find.hmm.states\
function")

    hmm <- hmm(aCGH.obj)
    mergeFunc(statesres = hmm$states.hmm[[model.use]],
              minDiff = minDiff)$states.hmm

}

hmm.merged <- function(aCGH.obj) aCGH.obj$hmm.merged
"hmm.merged<-" <-
    function(aCGH.obj, value)
{

    if (!is.aCGH(aCGH.obj))
	stop("object is not of class aCGH")
    aCGH.obj$hmm.merged <- value
    aCGH.obj

}

computeSD.Samples <-
    function(aCGH.obj, maxChrom = 22, maxmadUse = .3, maxmedUse = .5,
             maxState = 3, minClone = 20)
{

    if (is.null(hmm.merged(aCGH.obj)))
        stop("merge the hmm states first using merge.hmm.states\
function")

    ##computing SD of the tumor and sd on individual chromosomes

    computeSD.func(statesres = hmm.merged(aCGH.obj),
                   maxmadUse = maxmadUse, maxmedUse = maxmedUse,
                   maxState = maxState, minClone = minClone,
                   maxChrom = maxChrom)
    
}

sd.samples <- function(aCGH.obj) aCGH.obj$sd.samples
"sd.samples<-" <-
    function(aCGH.obj, value)
{

    if (!is.aCGH(aCGH.obj))
	stop("object is not of class aCGH")
    if (!is.null(aCGH.obj$sd.samples))
    {
        
        sd.samples.ok <-
            length(value) == length(aCGH.obj$sd.samples) &&
        all(
            sapply(1:length(aCGH.obj$sd.samples),
                   function(i)
                   all(dim(aCGH.obj$sd.samples[[i]]) == dim(value[[i]]))
                   )
                )
        if (!sd.samples.ok)
            stop("invalid replacement dimensions")

    }
    aCGH.obj$sd.samples <- value
    aCGH.obj

}

find.genomic.events <-
    function(aCGH.obj, maxChrom = 23, factor = 5, maxClones = 1,
             maxLen = 1000, absValSingle = 1, absValRegion = 1,
             diffVal1 = 1, diffVal2 = .5, maxSize = 10000,
             pChrom.min = .9, medChrom.min = .1)
{

    if (is.null(hmm.merged(aCGH.obj)))
        stop("merge the hmm states first using merge.hmm.states\
function")
    if (is.null(sd.samples(aCGH.obj)))
        stop("compute the std. deviations of aCGH samples using\
computeSD.Samples function")
    l2r <- log2.ratios(aCGH.obj)
    clones.info <- clones.info(aCGH.obj)
    sd.samples <- sd.samples(aCGH.obj)
    ncols <- ncol(l2r)
    statesMatrix <- hmm.merged(aCGH.obj)
    
    ##identifies outliers (factor times SD from the the median of the state)
    cat("Finding outliers\n")
    outliers <-
        findOutliers.func(thres = sd.samples$madGenome,
                          factor = factor, statesres = statesMatrix)
    ##identifies focal low level aberrations
    cat("Finding focal low level aberrations\n")
    aberrations <-
        findAber.func(maxClones = maxClones, maxLen = maxLen,
                      statesres = statesMatrix)
    ##identifies transitions: start and end of the states  
    cat("Finding transitions\n")
    transitions <-
        findTrans.func(outliers = outliers$outlier,
                       aber = aberrations$aber,
                       statesres = statesMatrix)
    ##identifies focal amplifications
    cat("Finding focal amplifications\n")
    amplifications <-
        findAmplif.func(absValSingle = absValSingle,
                        absValRegion = absValRegion,
                        diffVal1 = diffVal1, diffVal2 = diffVal2,
                        maxSize = maxSize,
                        translen.matr = transitions$translen.matrix,
                        trans.matr = transitions$trans.matr,
                        aber = aberrations$aber,
                        outliers = outliers$outlier,
                        pred = outliers$pred.out,
                        pred.obs = outliers$pred.obs.out,
                        statesres = statesMatrix)
    
    ##number of transitions per chromosome
    numTrans <- matrix(0, nrow = maxChrom, ncol = ncols)
    ##number of amplifications per chromosome
    numAmplif <- matrix(0, nrow = maxChrom, ncol = ncols)
    ##number of aberrations per chromosome
    numAber <- matrix(0, nrow = maxChrom, ncol = ncols)
    ##number of outliers per chromosome
    numOutlier <- matrix(0, nrow = maxChrom, ncol = ncols)
    ##number of chromosomes containing at least one transition
    numTrans.binary <- matrix(0, nrow = maxChrom, ncol = ncols)
    ##number of chromosomes containing at least one focal amplification
    numAmplif.binary <- matrix(0, nrow = maxChrom, ncol = ncols)
    ##number of chromosomes containing at least one focal aberration
    numAber.binary <- matrix(0, nrow = maxChrom, ncol = ncols)
    ##number of chromosomes containing at least one outlier
    numOutlier.binary <- matrix(0, nrow = maxChrom, ncol = ncols)
    ## whole chromosome gain or loss
    wholeChromGainLoss <- matrix(0, nrow = maxChrom, ncol = ncols)
    sizeAmplicon <- numAmplicon <-
        matrix(0, nrow = maxChrom, ncol = ncols)
    propClones <-  matrix(0, nrow = maxChrom, ncol = ncols)
    pvClones <-  matrix(0, nrow = maxChrom, ncol = ncols)
    medClones <- matrix(0, nrow = maxChrom, ncol = ncols)
    p.min <- pChrom.min
    pv.max <- .0001
    med.min <- medChrom.min
    chr <- clones.info(aCGH.obj)$Chrom
    kb <- clones.info(aCGH.obj)$kb

    for (j in 1:maxChrom)
    {

        cat("Processing chromosome ", j, "\n")
        ind <- chr == j
        trans <- transitions$trans.matrix[ ind, ,drop = FALSE]
        amplif <- amplifications$amplif[ ind, ,drop = FALSE]
        aber <- aberrations$aber[ ind, ,drop = FALSE]
        outlier <- outliers$outlier[ ind, ,drop = FALSE]
        for (i in 1:ncols)
        {

            numTrans[j, i] <- sum(trans[ ,i] == 1, na.rm=TRUE)
            if (numTrans[j, i] > 0)
                numTrans.binary[j, i] <- 1
            else # if no transitions
            {

                ##observed values
                obs <- l2r[ind, i]
                ##exclude aberrations and outliers
                obs <- obs[aber[ ,i ] == 0 & outlier[ ,i ] == 0]
                ##exclude missing values
                obs <- na.omit(obs)
                p <-
                    if (median(obs) >= 0)
                        mean(obs >= 0)
                    else
                        mean(obs < 0)
                propClones[j, i] <- p
                pv <-
                    1 - pnorm((p - .5) / sqrt((.5 ^ 2) / length(obs)))
                pvClones[j, i] <- pv
                medClones[j, i] <- median(obs)
                if ((p >= p.min) && (abs(median(obs)) >= med.min))
                {
                    
                    if (pv <= pv.max)
                        wholeChromGainLoss[j, i] <- 
                            if (median(obs) >= 0)
                                1
                            else
                                -1
                    
                }
                else
                    wholeChromGainLoss[j, i] <- 0
                
            }
            numAmplif[j,i] <- sum(amplif[ ,i ] == 1, na.rm = TRUE)
            if (numAmplif[j,i] > 0)
                numAmplif.binary[j,i] <- 1
            numAber[j,i] <- sum(aber[ ,i ] == 1, na.rm = TRUE)
            if (numAber[j,i] > 0)
                numAber.binary[j,i] <- 1
            numOutlier[j,i] <- sum(outlier[ ,i ] == 1, na.rm = TRUE)
            if (numOutlier[j,i] > 0)
                numOutlier.binary[j,i] <- 1
	    amp <- amplif[,i]
	    ind.na <- which(is.na(amp))
	    amp <- amp[-ind.na]
            try1 <- diff(amp)
            tmps <- which(try1 == 1) + 1
            tmpe <- which(try1 == -1)
            if (length(tmps) > length(tmpe))
                ##last clone
                tmpe <- c(tmpe, length(ind))
            if (length(tmps) < length(tmpe))
                ##first clone
                tmps <- c(1, tmps)
            if (length(tmpe) == length(tmps))
            {

                kb.ind <- kb[ind][-ind.na]
                tmplen <-
                    (kb.ind[tmpe] - kb.ind[tmps]) +
                        rep(1000, length(tmpe))
                sizeAmplicon[j, i] <- sum(tmplen)
                numAmplicon[j, i] <- length(tmpe)
                
            }
            
	}
        
    }

    list(num.transitions = numTrans,
         num.amplifications = numAmplif,
         num.aberrations = numAber,
         num.outliers = numOutlier,
         num.transitions.binary = numTrans.binary,
         num.amplifications.binary = numAmplif.binary,
         num.aberrations.binary = numAber.binary,
         num.outliers.binary = numOutlier.binary,
         whole.chrom.gain.loss = wholeChromGainLoss,
         size.amplicons = sizeAmplicon,
         num.amplicons = numAmplicon,
         outliers = outliers,
         aberrations = aberrations,
         transitions = transitions,
         amplifications = amplifications
         )
    
}

genomic.events <- function(aCGH.obj) aCGH.obj$genomic.events
"genomic.events<-" <-
    function(aCGH.obj, value)
{

    if (!is.aCGH(aCGH.obj))
	stop("object is not of class aCGH")
    if (!is.null(aCGH.obj$genomic.events))
    {
        
        events.ok <-
            length(value) == length(aCGH.obj$genomic.events) &&
        all(
            sapply(1:length(aCGH.obj$genomic.events),
                   function(i)
                   all(dim(aCGH.obj$genomic.events[[i]]) ==
                       dim(value[[i]]))
                   )
            )
        if (!events.ok)
            stop("invalid replacement dimensions")

    }
    aCGH.obj$genomic.events <- value
    aCGH.obj

}

phenotype <- function(aCGH.obj) aCGH.obj$phenotype
"phenotype<-" <-
    function(aCGH.obj, value)
{

    if (!is.aCGH(aCGH.obj))
	stop("object is not of class aCGH")
    if (nrow(value) != ncol(aCGH.obj$log2.ratios))
        stop("number of observations differs between the old and new\
phenotypes")
    aCGH.obj$phenotype <- value
    aCGH.obj

}

subset.hmm <-
    function(x, ...)
{

    if (is.null(x))
        return(NULL)
    ll <- list(...)
    i <- 
        if (is.null(ll$i))
            1:nrow(x$states.hmm[[1]])
        else
            ll$i
    j <- 
        if (is.null(ll$j))
            1:ncol(x$nstates.hmm[[1]])
        else
            ll$j
    chroms <- 
        if (is.null(ll$chroms))
            1:nrow(x$nstates.hmm[[1]])
        else
            ll$chroms
    with(x,
         list(nstates.hmm =
                   lapply(nstates.hmm,
                          function(nstates) nstates[chroms ,j]
                          ),
                   states.hmm =
                   lapply(states.hmm,
                          function(states)
                          states[i,
                                 c(1:2,
                                   as.vector(
                                             sapply(j,
                                                    function(k)
                                                    (3 + (k - 1) * 6):(2 + k * 6)
                                                    )
                                             )
					)
                                   
                                 ]
                          )
                   )
              )
         

}

subset.hmm.merged <-
    function(x, ...)
{

    if (is.null(x))
        return(NULL)

    ll <- list(...)
    i <- 
        if (is.null(ll$i))
            1:nrow(x$states.hmm[[1]])
        else
            ll$i
    j <- 
        if (is.null(ll$j))
            1:ncol(x$nstates.hmm[[1]])
        else
            ll$j
    
    x[i,
      c(1:2,
        as.vector(
                  sapply(j,
                         function(k)
                         (3 + (k - 1) * 6):(2 + k * 6)
                         )
                  )
        )
      ]
    
}

"[.aCGH" <-
    function(aCGH.obj, i, j, keep = FALSE)
{
    
    drop.i <- missing(i)
    drop.j <- missing(j)
    if (drop.i && drop.j)
        aCGH.obj
    else
    {

        if (drop.i)
            i <- 1:nrow(aCGH.obj)
        else
            if (mode(i) == "logical")
                i <- which(i)
        if (drop.j)
            j <- 1:ncol(aCGH.obj)
        else
            if (mode(j) == "logical")
                j <- which(j)
        res <-
            if (keep)
                list(log2.ratios = aCGH.obj$log2.ratios[i, j],
                     clones.info = aCGH.obj$clones.info[ i, ],
                     qual.rep = NULL,
                     bad.quality.index = NULL,
                     log2.ratios.imputed =
                     if (is.null(aCGH.obj$log2.ratios.imputed)) NULL
                     else aCGH.obj$log2.ratios.imputed[i, j],
                     sd.samples =
                     if (is.null(aCGH.obj$sd.samples)) NULL
                     else 
                     with(aCGH.obj$sd.samples,
                          list(madChrom = madChrom[ ,j ],
                               madGenome = madGenome[j]
                               )
                          ),
                     genomic.events =
                     if (is.null(aCGH.obj$genomic.events)) NULL
                     else 
                     lapply(aCGH.obj$genomic.events,
                            function(el)
                            if (is.matrix(el)) el[ ,j ]
                            else
                            lapply(el, function(el1) el1[i, j])
                            ),
                     hmm =
                     subset.hmm(hmm(aCGH.obj), i = i, j = j,
                                chroms =
                                which(table(clones.info(aCGH.obj)$Chrom[i]) > 0)
                                ),
                     hmm.merged =
                     subset.hmm.merged(hmm.merged(aCGH.obj), i = i,
                                       j = j),
                     phenotype =
                     if (is.null(aCGH.obj$phenotype)) NULL
                     else aCGH.obj$phenotype[j, , drop = FALSE]
                     )
            else
            {
                
                #warning("For now just subsetting the log2.ratios\
#and phenotype. Please rerun the find.hmm.states function!")
		warning("subsetting the log2.ratios only")
                list(log2.ratios =
                     aCGH.obj$log2.ratios[i, j, drop = FALSE],
                     clones.info =
                     aCGH.obj$clones.info[i, , drop = FALSE],
                     qual.rep = NULL,
                     bad.quality.index = NULL,
                     log2.ratios.imputed =
                     if (is.null(aCGH.obj$log2.ratios.imputed)) NULL
                     else aCGH.obj$log2.ratios.imputed[i, j, drop = FALSE],
                     sd.samples = NULL,
                     genomic.events = NULL,
                     hmm = NULL,
                     phenotype =
                     if (is.null(aCGH.obj$phenotype)) NULL
                     else aCGH.obj$phenotype[j, , drop = FALSE]
                     )

            }
        attr(res, "call") <- sys.call()
        class(res) <- "aCGH"
        res

    }
    
}

print.aCGH <-
    function(x, ...)
{

    cat("aCGH object\nCall: ")
    print(attr(x, "call"), ...)
    cat("\nNumber of Arrays", ncol(x),
        "\nNumber of Clones", nrow(x), "\n")

}

summary.aCGH <-
    function(object, ...)
{
    
    print.aCGH(object, ...)
    if (!is.null(log2.ratios.imputed(object)))
        cat("Imputed data exist\n")
    else
        cat("Imputed data does not exist\n")
    if (!is.null(hmm(object)))
        cat("HMM states assigned\n")
    else
        cat("HMM states are not assigned\n")
    if (!is.null(sd.samples(object)))
        cat("samples standard deviations are computed\n")
    else
        cat("samples standard deviations are not computed\n")
    if (!is.null(genomic.events(object)))
        cat("genomic events are assigned\n")
    else
        cat("genomic events are not assigned\n")
    if (!is.null(phenotype(object)))
        cat("phenotype exists\n")
    else
        cat("phenotype does not exists\n")

}

plot.aCGH <-
    function(x, ...)
{

    ll <- list(...)
    #dat <- 
    #    if (!is.null(ll$imp) && ll$imp)
    #        as.matrix(log2.ratios.imputed(x))
    #    else
    #        as.matrix(log2.ratios(x))
    dat <- as.matrix(log2.ratios(x))
    Colv <- 
        if (!is.null(ll$Colv) && ll$Colv)
            ll$Colv
        else
            NULL
###    heatmap(dat, Rowv = NA, Colv = Colv, main = "Heatmap",
###            labCol = sample.names(x))
###    if (!is.null(genomic.events(x)))
###        plotSummaryProfile(x)
    plotFreqStat(x)
    
}

minna <-
    function(x)
    min(x, na.rm = TRUE)

maxna <-
    function(x)
    max(x, na.rm = TRUE)

corna <-
    function(x)
    cor(x, use = "pairwise.complete.obs")

floor.func <-
    function(x, floor, x.na = x[!is.na(x)])
{
    x[!is.na(x)] <-
        ifelse(x.na > floor,
               floor,
               ifelse(x.na < -floor, -floor, x.na)
               )
    x
    
}

length.num.func <-
    function(x, num)
    sapply(num, function(i) sum(x == i, na.rm = TRUE))

prop.num.func <-
    function(x, num)
    sapply(num, function(i) mean(x == i, na.rm = TRUE))

as.eSet <-
    function(aCGH.obj)
{

    ver <- R.Version()
    if (as.numeric(ver$major) < 1 && as.numeric(ver$minor) < 9)
        stop("Using as.eSet() requires R version >= 1.9!")
    
    new.l2r <-
        new("exprList",
            .Data =
            list(exprs = log2.ratios(aCGH.obj),
                 log2.ratios.imputed = log2.ratios.imputed(aCGH.obj),
                 clones.info = clones.info(aCGH.obj),
                 hmm = hmm(aCGH.obj),
                 hmm.merged = hmm.merged(aCGH.obj),
                 sd.samples = sd.samples(aCGH.obj),
                 genomic.events = genomic.events(aCGH.obj)),
            eMetadata =
            data.frame(name =
                       c("log2 ratios", "log2 ratios imputed",
                         "clones info", "hmm states",
                         "merged hmm states",
                         "samples noise std. deviation",
                         "genomic events"),
                       etype =
                       c("random numbers", "imputed numbers",
                         "clone information", "integers",
                         "integers", "random numbers",
                         "integers")))
    pheno <- phenotype(aCGH.obj)
    pheno.names <- colnames(pheno)
    varLabels <-
        lapply(1:ncol(pheno),
               function(i) {

                   attr.class <- class(pheno[, i])
                   descr <-
                       if (attr.class == "factor")
                           paste(nlevels(pheno[, i]), "levels")
                       else
                           attr.class
                   
                   paste(pheno.names[i], "; ", descr, sep = "")
                   
               })
    names(varLabels) <- paste("cov", 1:ncol(pheno), sep = "")
    phenoData <-
        new('phenoData',
            pData = pheno,
            varLabels = varLabels,
            varMetadata = data.frame(varNames = pheno.names))
    
    new("eSet", eList = new.l2r, phenoData = phenoData)
    
}
