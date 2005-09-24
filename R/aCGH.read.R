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
                              blank.lines.skip = FALSE, dec = ".")
               colnames(dt.tmp) <- dotify.names(colnames(dt.tmp))
               dat <- dt.tmp[-nrow(dt.tmp), cols]
               log2rat <- as.numeric(dat[, 1])
               log2stddev <- as.numeric(dat[, 2])
               nreplic <- as.numeric(dat[ ,3 ])
               flag <- as.numeric(dat[ ,4 ])
               tmp1 <-
                   flag == 1 | (log2stddev > maxsd) | (nreplic < minreplic)
               log2rat[tmp1] <- NA
               log2rat
               
           }
           )

extract.clones.info <-
    function(dt.tmp)
{
    
    dt.tmp <- dt.tmp[ -nrow(dt.tmp), ]
    colnames(dt.tmp) <- dotify.names(colnames(dt.tmp))
    clones.info <-
        dt.tmp[ , c("Clone", "Target", "Chromosome", "KB.POSITION")]
    
}

maxdiff.func <-
    function(x)
    abs(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
mincorr.func <-
    function(A)
    min(cor(A, use = "pair", method = "spearman"))

aCGH.read.Sprocs <-
    function(fnames, latest.mapping.file = NULL, maxsd = .2,
             minreplic = 2, chrom.remove.threshold = 24,
             prop.missing = .25, sample.names = fnames,
             sample.quality.threshold = .4,
             cols = c("Log2Rat", "Log2StdDev", "NReplic", "Bad.P"),
             unmapScreen = TRUE, dupRemove = TRUE)
{
    
    ## screening out clones with < 2 replicates or SD > .2 or
    ## SPROC indicator of 1

    log2.ratios <-
        read.Sproc.files(fnames, maxsd = maxsd, minreplic = minreplic,
                         cols = cols)
###    colnames(log2.ratios) <- sample.names
    	
    ## Extract the clones information from the first file in the list
    
    clones.info <- 
        extract.clones.info(read.table(fnames[1], h = TRUE, sep = "\t",
                                       quote = "", comment.char = "",
                                       fill = TRUE,
                                       blank.lines.skip = FALSE))
    
    ## if clones have newer mapping associated with them

    if (!is.null(latest.mapping.file))
    {

        latest.mapping <-
            read.table(latest.mapping.file, sep = "\t", h = TRUE,
                       quote = "", comment.char = "")[, 1:4]
        colnames(latest.mapping) <-
            dotify.names(colnames(latest.mapping))
        ind.match <-
            match(as.character(clones.info$Clone),
                  as.character(latest.mapping$USER.CLONE.ID))
        ind <- ind.match[!is.na(ind.match)]
        clones.info <- latest.mapping[ind ,]
        log2.ratios <- log2.ratios[!is.na(ind.match), , drop = FALSE]

    }
    colnames(clones.info) <- c("Clone", "Target", "Chrom", "kb")
    rownames(log2.ratios) <- clones.info$Target
    
    if (unmapScreen)
    {
        
        ## remove unmapped clones

    	ind.unmap <-
            which(clones.info$Chrom > chrom.remove.threshold |
                  is.na(clones.info$Chrom) | is.na(clones.info$kb))
	if (length(ind.unmap) > 0)
	{
            
            clones.info <- clones.info[-ind.unmap ,]
            log2.ratios <- log2.ratios[-ind.unmap, , drop = FALSE]
            
	}
        
    }

    ## reorder by chromosome and chromosomal position

    ord <- order(clones.info$Chrom, clones.info$kb)
    clones.info <- clones.info[ord ,]
    log2.ratios <- log2.ratios[ord , , drop = FALSE]

    ## mark samples that have bad quality

    bad.quality.index <-
        which(apply(log2.ratios,
                    2,
                    function(col)
                    mean(is.na(col)) > sample.quality.threshold))

    ## screen out clones missing in > prop.missing% of the samples:

    prop.miss <-
        apply(log2.ratios, 1, prop.na)
    clones.info <-
        clones.info[prop.miss <= prop.missing ,]
    log2.ratios <-
        log2.ratios[prop.miss <= prop.missing, , drop = FALSE]

    ## determine duplicates and average/remove them

    tbl <- table(clones.info$Clone)
    qual.rep <- NULL
    if (any(tbl > 1))
    {
        
        tbl <- tbl[tbl > 1]
        qual.rep <- matrix(0, ncol = length(tbl), nrow = 2)
        nms <- names(tbl)
        cat("\nAveraging duplicated clones\n")
        for (i in 1:length(tbl))
        {

            ind1 <- which(clones.info$Clone == nms[i])
            cat(as.character(clones.info$Clone[ind1[1]]),
                "\t", ind1, "\n")
            lr <- log2.ratios[ind1, , drop = FALSE]
            vec <- apply(lr, 2, mean, na.rm = TRUE)
            md <- apply(lr, 2, maxdiff.func)
            md <- md[md > 0]
            qual.rep[1, i] <- round(median(md, na.rm = TRUE), 2)
            qual.rep[2, i] <- round(mincorr.func(t(lr)), 2)
	    if (dupRemove)
            	for (j in 1:length(ind1))
                    log2.ratios[ind1[j] ,] <- vec
            
        }
        qual.rep <- rbind(tbl, qual.rep)

    }
    ## contains median of abs max difference among all replicates and
    ## min correlations between replicates

    if (dupRemove)
    {
        
	dupl  <- duplicated(clones.info$Clone)
    	clones.info <- clones.info[ !dupl, ]
    	log2.ratios <- log2.ratios[ !dupl, , drop = FALSE]
        
    }
    if (!is.null(sample.names))
        colnames(log2.ratios) <- sample.names

    clones.info$Clone <- factor(clones.info$Clone)
    clones.info$Target <- factor(clones.info$Target)
    value <- create.aCGH(log2.ratios, clones.info)
    attr(value, "call") <- sys.call()
    value$qual.rep <- qual.rep
    value$bad.quality.index <- bad.quality.index
    value

}

###########################

aCGH.process <-
    function(aCGH.obj, chrom.remove.threshold = 24,
             prop.missing = .25, sample.quality.threshold = .4,
             unmapScreen=TRUE, dupRemove = TRUE)
{
    
    log2.ratios <- log2.ratios(aCGH.obj)
    
    ## Extract the clones information from the first file in the list
    
    clones.info <- clones.info(aCGH.obj)
            
    if (unmapScreen)
    {
    ## remove unmapped clones
    	ind.unmap <-
        which(clones.info$Chrom > chrom.remove.threshold |
              is.na(clones.info$Chrom) | is.na(clones.info$kb))
	if (length(ind.unmap) > 0)
	{
    		clones.info <- clones.info[ -ind.unmap, ]
    		log2.ratios <- log2.ratios[ -ind.unmap, , drop = FALSE]
	}
     }
    ## reorder by chromosome and chromosomal position
    ord <- order(clones.info$Chrom, clones.info$kb)
    clones.info <- clones.info[ ord, ]
    log2.ratios <- log2.ratios[ ord, , drop = FALSE]

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
    clones.info <-
        clones.info[ prop.miss <= prop.missing, ]
    log2.ratios <-
        log2.ratios[ prop.miss <= prop.missing, , drop = FALSE]

    ## determine duplicates and average/remove them

    tbl <- table(clones.info$Clone)
    qual.rep <- NULL
    if (any(tbl > 1))
    {
        
        tbl <- tbl[tbl > 1]
        qual.rep <- matrix(0, ncol = length(tbl), nrow = 2)
        nms <- names(tbl)
	if (dupRemove)
            cat("\nAveraging duplicated clones\n")
	for (i in 1:length(tbl))
        {

            ind1 <- which(clones.info$Clone == nms[i])
            cat(as.character(clones.info$Clone[ind1[1]]),
                "\t", ind1, "\n")
            vec <-
                apply(log2.ratios[ ind1, , drop = FALSE],
                      2,
                      mean,
                      na.rm = TRUE)
            md <-
                apply(log2.ratios[ ind1, , drop = FALSE],
                      2,
                      maxdiff.func)
            md <- md[md > 0]
            qual.rep[1, i] <- round(median(md, na.rm = TRUE), 2)
            qual.rep[2, i] <-
                round(mincorr.func( t(log2.ratios[ ind1, , drop = FALSE]) ), 2)
	    if (dupRemove)
            	for (j in 1:length(ind1))
                    log2.ratios[ ind1[j], , drop = FALSE] <- vec
            
        }
        qual.rep <- rbind(tbl, qual.rep)

    }
    ## contains median of abs max difference among all replicates and
    ## min correlations between replicates

    if (dupRemove)
    {
        
	dupl  <- duplicated(clones.info$Clone)
    	clones.info <- clones.info[ !dupl, ]
    	log2.ratios <- log2.ratios[ !dupl, , drop = FALSE]
        
    }

    clones.info$Clone <- factor(clones.info$Clone)
    clones.info$Target <- factor(clones.info$Target)
    value <- create.aCGH(log2.ratios, clones.info, phenotype(aCGH.obj))
    attr(value, "call") <- sys.call()
    value$qual.rep <- qual.rep
    value$bad.quality.index <- bad.quality.index
    value

}

###########################

