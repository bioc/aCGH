dotify.names <-
    function(nms)
    gsub("_", ".", nms)

read.Sproc.files <-
    function(fnames, maxsd = .2, minreplic = 2,
             cols = c("Log2Rat", "Log2StdDev", "NReplic", "Bad.P"))
    sapply(fnames,
           function(fname) {
               
               cat("Trying to read ", fname, "\n")
               dt.tmp <-
                   read.table(fname, h = TRUE, sep = "\t", quote = "",
                              comment.char = "", fill = TRUE,
                              blank.lines.skip = FALSE)
               colnames(dt.tmp) <- dotify.names(colnames(dt.tmp))
               dat <- dt.tmp[-nrow(dt.tmp), cols]
               log2rat <- dat[ ,1 ]
               log2stddev <- dat[ ,2 ]
               nreplic <- dat[ ,3 ]
               flag <- dat[ ,4 ]
               tmp1 <-
                   flag == 1 &
               ((log2stddev > maxsd) | (nreplic < minreplic))
               log2rat[tmp1] <- NA
               log2rat
               
           }
           )

##flag.func <-
##    function(dat, maxsd = .2, minreplic = 2, colvals = 1, colsd = 2,
##             colrep = 3, colbad = 4)
##{

##    seq.val <- seq(colvals, ncol(dat), by = 4)
##    seq.sd <- seq(colsd, ncol(dat), by = 4)
##    seq.rep <- seq(colrep, ncol(dat), by = 4)
##    seq.bad <- seq(colbad, ncol(dat), by = 4)

##    log2rat <- as.matrix(dat[ ,seq.val ])
##    log2stddev <- as.matrix(dat[ ,seq.sd ])
##    nreplic <- dat[ ,seq.rep ]
##    flag <- dat[ ,seq.bad ]

##    for (i in 1:ncol(flag))
##    {
        
#####        tmp <- flag[,i]
#####        tmp[((log2stddev[,i] > maxsd) | (nreplic[,i] < minreplic))] <- 1
#####        flag[ ,i ] <- tmp
#####        log2rat[tmp == 1, i] <- NA
##        tmp1 <- flag[,i] == 1 & ((log2stddev[,i] > maxsd) | (nreplic[,i] < minreplic))
##        log2rat[tmp1, i] <- NA
#####        log2stddev[tmp == 1, i] <- NA
        
##    }
    
##    log2rat
    
##}

extract.clones.info <-
    function(dt.tmp)
{
    
    dt.tmp <- dt.tmp[ -nrow(dt.tmp), ]
    colnames(dt.tmp) <- dotify.names(colnames(dt.tmp))
    clones.info <-
        dt.tmp[ , c("Clone", "Target", "Chromosome", "KB.POSITION")]
    
}

aCGH.read.Sprocs <-
    function(fnames, latest.mapping.file = NULL, maxsd = .2,
             minreplic = 2, chrom.remove.threshold = 24,
             prop.missing = .25, sample.names = fnames,
             sample.quality.threshold = .4,
             cols = c("Log2Rat", "Log2StdDev", "NReplic", "Bad.P"))
{
    
    maxdiff.func <-
        function(x)
            abs(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

    mincorr.func <-
        function(A)
        {
            
            crmn <- 2
            for (i in 1:(ncol(A) - 1))
                for (j in (i+1):ncol(A))
                {
                    
                    vec1 <- A[ ,i ]
                    vec2 <- A[ ,j ]
                    ind <- which(!is.na(vec1) & !is.na(vec2))
                    vec1 <- rank(vec1[ind])
                    vec2 <- rank(vec2[ind])
                    cr <- cor(vec1, vec2)
                    if (cr < crmn)
                        crmn <- cr
                    
                }
            
            crmn
            
        }

###    if (is.null(sample.names))
###        sample.names <- 
###            sapply(strsplit(fnames, "/"),
###                   function(vv) {

###                       name.split <-
###                           strsplit(vv[[length(vv)]], "\\.")[[1]]
###                       num.splits <- length(name.split)
###                       if (num.splits > 2)
###                           name.split <-
###                               paste(name.split[-num.splits], ".")
###                       name.split

###                   }
###                   )

    ## screening out clones with < 2 replicates or SD > .2 or
    ## SPROC indicator of 1

    log2.ratios <-
        read.Sproc.files(fnames, maxsd = maxsd, minreplic = minreplic,
                         cols = cols)
    colnames(log2.ratios) <- sample.names

    ## Extract the clones information from the first file in the list
    
    clones.info <- 
        extract.clones.info(read.table(fnames[1], h = TRUE, sep = "\t",
                                       quote = "", comment.char = "",
                                       fill = TRUE,
                                       blank.lines.skip = FALSE)
                            )
    
    ## if clones have newer mapping associated with them
    if (!is.null(latest.mapping.file))
    {

        latest.mapping <-
            read.table(latest.mapping.file, sep = "\t", h = TRUE,
                       quote = "", comment.char = "")[ ,1:4 ]
        colnames(latest.mapping) <-
            dotify.names(colnames(latest.mapping))
        ind.match <-
            match(as.character(clones.info$Clone),
                  as.character(latest.mapping$USER.CLONE.ID)
                  )
        ind <- ind.match[!is.na(ind.match)]
        clones.info <- latest.mapping[ ind, ]
        log2.ratios <- log2.ratios[ !is.na(ind.match), ]

    }
    colnames(clones.info) <- c("Clone", "Target", "Chrom", "kb")
    rownames(log2.ratios) <- clones.info$Clone
    
    ## remove unmapped clones
    ind.unmap <-
        which(clones.info$Chrom > chrom.remove.threshold |
              is.na(clones.info$Chrom) | is.na(clones.info$kb))
    clones.info <- clones.info[ -ind.unmap, ]
    log2.ratios <- log2.ratios[ -ind.unmap, ]

    ## reorder by chromosome and chromosomal position
    ord <- order(clones.info$Chrom, clones.info$kb)
    clones.info <- clones.info[ ord, ]
    log2.ratios <- log2.ratios[ ord, ]

    ## mark those samples that have bad quality
    bad.quality.index <-
        which(apply(log2.ratios,
                    2,
                    function(col)
                    mean(is.na(col)) > sample.quality.threshold
                    )
              )

    ## screen out clones missing in > prop.missing% of the samples:

    prop.miss <- apply(log2.ratios, 1, prop.na)
    clones.info.screen <-
        clones.info[ prop.miss <= prop.missing, ]
    log2.ratios <- log2.ratios[ prop.miss <= prop.missing, ]

    ## determine duplicates and average/remove them

    tbl <- table(clones.info.screen$Clone)
    qual.rep <- NULL
    if (any(tbl > 1))
    {
        
        tbl <- tbl[tbl > 1]
        qual.rep <- matrix(0, ncol = length(tbl), nrow = 2)
        nms <- names(tbl)
        cat("\nAveraging duplicated clones\n")
        for (i in 1:length(tbl))
        {

            ind1 <- which(clones.info.screen$Clone == nms[i])
            cat(as.character(clones.info.screen$Clone[ind1[1]]),
                "\t", ind1, "\n")
            vec <- apply(log2.ratios[ ind1, ], 2, mean, na.rm = TRUE)
            md <- apply(log2.ratios[ ind1, ], 2, maxdiff.func)
            md <- md[md > 0]
            qual.rep[1, i] <- round(median(md, na.rm = TRUE), 2)
            qual.rep[2, i] <-
                round(mincorr.func( t(log2.ratios[ ind1, ]) ), 2)
            for (j in 1:length(ind1))
                log2.ratios[ ind1[j], ] <- vec
            
        }
        qual.rep <- rbind(tbl, qual.rep)

    }
    ## contains median of abs max difference among all replicates and
    ## min correlations between replicates

    dupl  <- duplicated(clones.info.screen$Clone)
    clones.info.screen <- clones.info.screen[ !dupl, ]
    log2.ratios <- log2.ratios[ !dupl, ]

    if (!is.null(sample.names))
        colnames(log2.ratios) <- sample.names

    clones.info.screen$Clone <- factor(clones.info.screen$Clone)
    clones.info.screen$Target <- factor(clones.info.screen$Target)
    value <- create.aCGH(log2.ratios, clones.info.screen)
    attr(value, "call") <- sys.call()
    value$qual.rep <- qual.rep
    value$bad.quality.index <- bad.quality.index
    value

}
