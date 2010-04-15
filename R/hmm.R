##performing HMM:

#####################################################################
#####################################################################

hmm.run.func <-
    function(dat, datainfo = clones.info, vr = .01, maxiter = 100,
             aic = TRUE, bic = TRUE, delta = NA, eps = 0.01)
{
    chrom.uniq <- unique(datainfo$Chrom)
    states <- matrix(NA, nrow=nrow(dat), ncol=(2+6*ncol(dat)))
    states[,1:2] <- cbind(datainfo$Chrom, datainfo$kb)
    nstates <- matrix(NA, nrow=length(chrom.uniq), ncol=ncol(dat))

    states.list <- list(states)
    nstates.list <- list(nstates)

    ##list consists of experiments starting with aic then, if bic, scroll
    ##over deltas. delta= 1 by default.

    nlists <- 0
    if (aic)
    {
        nlists <- 1
    }
    if (bic)
    {
        if (is.na(delta))
        {
            delta <- c(1)
        }
        else
        {
            delta <- c(1, delta)
        }
        for (j in 1:length(delta))
        {
            nlists <- nlists+1
        }
    }
    
    if (nlists > 1)
    {
        for (j in 2:nlists)
        {
            states.list[[j]] <- states.list[[1]]
            nstates.list[[j]] <- nstates.list[[1]]
        }
    }
    
    for (i in 1:ncol(dat))
    {
        
        cat("sample is ", i, "\tChromosomes: ")
        colstart <- 2+(i-1)*6+1
        colend <- 2+i*6
        for (j in 1:length(chrom.uniq))
        {
            
            cat(j, " ")
            
            res <-
                try(states.hmm.func(sample=i, chrom=j, dat=dat,
                                    datainfo=datainfo,vr=vr,
                                    maxiter=maxiter, aic=aic, bic=bic,
                                    delta=delta, nlists=nlists, eps = eps))
            if (!(inherits(res, "try-error")))
                for (m in 1:nlists)
                {

                    states.list[[m]][((1:nrow(states))[states[,1]==j]),colstart:colend] <-
                        as.matrix(res$out.list[[m]])
                    nstates.list[[m]][j,i] <- res$nstates.list[[m]]

                }
            else
                stop("hmm fit failed!\n")

        }
        cat("\n")
        
    }
    list(states.hmm = states.list, nstates.hmm = nstates.list)
}

states.hmm.func <-
    function(sample, chrom, dat, datainfo = clones.info, vr = .01,
             maxiter = 100, aic = FALSE, bic = TRUE, delta = 1,
             nlists = 1, eps = .01, print.info = FALSE,
             diag.prob = .99)
{

    obs <- dat[datainfo$Chrom==chrom, sample]
    kb <- datainfo$kb[datainfo$Chrom==chrom]
    ##with current sproc files, data is already ordered by kb's
    obs.ord <- obs[order(kb)]
    kb.ord <- kb[order(kb)]

    ind.nonna <- which(!is.na(obs.ord))

    y <- obs.ord[ind.nonna]
    kb <- kb.ord[ind.nonna]


#####################################

    numobs <- length(y)
    zz <- vector(mode = "list", 5)
    zz[[1]] <-
        list(log.lik =
             sum(dnorm(y, mean = mean(y), sd = sd(y), log = TRUE)))
    for(k in 2:5)
    {
        
        mu <- kmeans(y, k)$centers
        gamma <- matrix((1 - diag.prob) / (k - 1), k, k)
        diag(gamma) <- diag.prob
        zz[[k]] <-
        {

            res <-
                .C("calc_observed_likelihood",
                   as.integer(numobs),
                   as.double(y),
                   as.integer(k),
                   mu = as.double(mu),
                   sigma = as.double(sqrt(vr)),
                   gamma = as.double(gamma),
                   pi = as.double(rep(-log(k), k)),
                   num.iter = as.integer(maxiter),
                   as.double(eps),
                   log.lik = double(1),
                   filtered.cond.probs = double(k * numobs),
                   hidden.states = integer(numobs),
                   as.logical(print.info),
                   PACKAGE = "aCGH")
            res$hidden.states <- res$hidden.states + 1
            res$filtered.cond.probs <-
                matrix(res$filtered.cond.probs, nr = k)
            res$gamma <- matrix(res$gamma, nr = k)
            res
            
        }
        
    }

###############################################3
###############################################3
    ##identify the model with the smallest model selection criteria

    ##now, scroll over all options:

    ##number of states (means) + number of states*(number of states-1) (transitions) #+ 1 (variance)

#    kk <- c(2, 5, 10, 17, 26)
    kk <- (1:5) ^ 2 + 1
    for (nl in 1:nlists)
    {
        if ((aic) && (nl==1))
        {
            ##-loglik+2*k/2
            factor <- 2
        }
        else if (bic)
        {
            ##-loglik+log(n)*k*delta/2
            if (aic)
            {
                factor <- log(numobs)*delta[nl-1]
            }
            else
            {
                factor <- log(numobs)*delta[nl]
            }
        }
        lik <- sapply(zz, function(z) -z$log.lik) + kk * factor / 2
        nstates <- likmin <- which.min(lik)
        z <- zz[[likmin]]

######################################
        ##out rpred and state

        if (nstates > 1) #if non-generic
        {
            ##print(nstates)
            maxstate <- apply(z$filter, 2, which.max)
###            maxstate <- z$hidden.states
            rpred <- as.vector(z$mu %*% z$filter)
            prob <- apply(z$filter, 2, max)
            ##use median for prediction and mad for state dispersions
            maxstate.unique <- unique(maxstate)
            pred <- rep(0, length(y))
            disp <- rep(0, length(y))
            for (m in 1:length(maxstate.unique))
            {
                
                pred[maxstate==maxstate.unique[m]] <-
                    median(y[maxstate==maxstate.unique[m]])
                disp[maxstate==maxstate.unique[m]] <-
                    mad(y[maxstate==maxstate.unique[m]])
                
            }

            ##if (length(z$pshape) == 1)
            ##{
            ##        disp <- rep(z$pshape, length(maxstate))
            ##}
            ##else
            ##{
            ##        disp <- z$pshape[maxstate]
            ##}
            
        }
        else #if generic
        {
            maxstate <- rep(1, length(y))
            ##rpred <- rep(mean(y), length(y))
            rpred <- rep(median(y), length(y))
            prob <- rep(1, length(y))
            ##pred <- rep(mean(y), length(y))
            pred <- rep(median(y), length(y))
            ##disp <- rep(var(y), length(y))
            disp <- rep(mad(y), length(y))
            
        }
        
        out <-
            cbind(matrix(maxstate, ncol=1),
                  matrix(rpred, ncol=1),
                  matrix(prob, ncol=1),
                  matrix(pred, ncol=1),
                  matrix(disp, ncol=1))
        
        out.all <- matrix(NA, nrow=length(kb.ord), ncol=6)
        out.all[ind.nonna,1:5] <- out
        
        out.all[,6] <- obs.ord
        out.all <- as.data.frame(out.all)
        dimnames(out.all)[[2]] <- c("state", "rpred", "prob", "pred", "disp", "obs")
        
        
        if (nl==1)
        {
            out.all.list <- list(out.all)
            nstates.list <- list(nstates)
        }
        else
        {
            out.all.list[[nl]] <- out.all
            nstates.list[[nl]] <- nstates
        }
        
        ##cloneinfo <- as.data.frame(cbind(rep(chrom, length(kb.ord)), kb.ord))
        ##dimnames(cloneinfo)[[2]] <- c("Chrom", "kb")
    }
    list(out.list = out.all.list, nstates.list = nstates.list)
    
}

#####################################################################
#####################################################################

mergeFunc <-
    function(statesres, minDiff = .1)
{
    ##merging states which are too close to each other to avoid showing technical
    ##rather than biological artifacts.
    ##
    ##start with two states closest to each other, if there are close enough , merge them
    ##and continue while treating them as the same state. Stop when no more states
    ##can be merged. write out states.merge file

    ##with merging only states and predicted values need to be changed

    sq.state <- seq(3, ncol(statesres), b = 6)
    sq.pred <- seq(6, ncol(statesres), b = 6)
    sq.obs <- seq(8, ncol(statesres), b = 6)

    chrom <- statesres[,1]

    for (i in 1:length(sq.state)) #over all samples
    {
###        print(i)

        for (j in 1:length(unique(chrom))) ##over all chromosomes
        {
            ##missing values
            statesr <- statesres[ chrom == j, ]
            ind.nonna <- which(!is.na(statesr[ ,sq.obs[i] ]))
            states <- statesr[,sq.state[i] ][ind.nonna] #states
            obs <- statesr[ ,sq.obs[i] ][ind.nonna] #observations
            pred <- statesr[ ,sq.pred[i] ][ind.nonna] #predicted (medians in a state)

            ##if there are K states, there can't be more than (K-1) merges

            num.states <- length(unique(states)) ##number of states:

            ##if more than 1 state
            if (num.states > 1)
            {
                
                for (m in 1:(num.states-1))
                {

#########

                    states.uniq <- unique(states) ##unique list of states
                    pred.states.uniq <- rep(0, length(states.uniq)) ##unique list of predictions for that states
                    for (s in 1:length(states.uniq))
                    {

                        pred[states ==states.uniq[s]] <- median(obs[states ==states.uniq[s]])
                        pred.states.uniq[s] <- (pred[states==states.uniq[s]])[1]
                        
                    }
#########

                    dst <- abs(dist(pred.states.uniq)) ##paiwise distances
                    if (min(dst)  >= minDiff)
                        ##if minimum difference is large enough
                    {
                        ##record the merged version
                        statesres[chrom==j, sq.state[i]][ind.nonna] <-
                            states
                        statesres[chrom==j, sq.pred[i]][ind.nonna] <- pred
                        
                        break
                    }
                    else
                        ##if closest are close enough
                    {
                        pred.dist.matr <- as.matrix(dst)
                        ##find the closest two states, in the case of the tie take the first one
                        for (s1 in 1:(nrow(pred.dist.matr)-1))
                        {
                            for (s2 in (s1+1):ncol(pred.dist.matr))                                                        {
                                if (pred.dist.matr[s1,s2] == min(dst))
                                    ##s1 and s2 are  the first two states that are too close to each other
                                {
                                    states[states==states.uniq[s2]] <- states.uniq[s1] ##update states
                                    pred[states==states.uniq[s1]] <- median(obs[states==states.uniq[s1]])
                                    break
                                }
                            }
                        }
                    }
                }
                
            }
            statesres[chrom==j, sq.state[i]][ind.nonna] <- states
            statesres[chrom==j, sq.pred[i]][ind.nonna] <- pred
        }
    }
    list(states.hmm = statesres)
    
}

#################################################################################
#################################################################################

computeSD.func <-
    function(statesres, maxmadUse = .2, maxmedUse = .2,
             maxState=3, maxStateChange = 100, minClone=20,
             maxChrom=22)
{
    ##1. maxmadUse : use state only if its mad <= maxmadUse (to avoid cases where
    ##points HMM does not seprate reults into several states (e.g. when
    ##individual clones fall out)
    ##maxmedUse: use state only if median is less < maxmedUse (states
    ##with chnages may have higher sd beause not all clones acquare a given
    ##aberration
    ##2. maxState: use only those chromosomes that contain fewer than maxState
    ##states (may modify to state transitions later
    ##maxStateChange show chromosomes with fewer than ceratin number of changes
    ##3. minClone: use only those states that contain greater than minClone
    ##observations
    ##4. maxChrom: use up to maxChrom (e.g. avoid using X and Y)
    chrom <- statesres[,1]
    sq.state <- seq(3, ncol(statesres), b=6)
    sq.obs <- seq(8, ncol(statesres), b=6)
    madChrom <- matrix(NA, nrow=length(unique(chrom)), ncol=length(sq.state))
    madGenome <- rep(NA, length(sq.state))
    for (i in 1:length(sq.obs))
    {

        mad.tmp <- NA
        for (j in 1:maxChrom)
        {
            ind.nonna <- (1:length(statesres[chrom==j, sq.obs[i]]))[!is.na(statesres[chrom==j, sq.obs[i]])]
            mad.tmp1 <- NA
            states.uniq <- unique(statesres[chrom==j, sq.state[i]][ind.nonna])
            states <- statesres[chrom==j, sq.state[i]][ind.nonna]
            obs <- statesres[chrom==j, sq.obs[i]][ind.nonna]
            ##use only chromosomes with <= maxState

            states.change <- diff(states)
            states.change.num <- length(states.change[states.change!=0])
            
            ##if (length(states.uniq) <= maxState)
            if ((length(states.uniq) <= maxState) && (states.change.num <= maxStateChange))
            {
                
                for (k in states.uniq)
                {
                    obs.state <- obs[states==k]
                    md <- mad(obs.state, na.rm = TRUE)
                    med <- median(obs.state, na.rm = TRUE)
                    
                    ##use only states with >=  minClone and mad of state is <= maxmadUse and
                    ##median of state is < maxmedUse
                    if ((length(obs.state) >= minClone) && (md <= maxmadUse) && (abs(med) <= maxmedUse))
                    {
                        
                        mad.tmp1 <- c(mad.tmp1, md)
                        
                    }
                }
            }
            if (length(mad.tmp1[!is.na(mad.tmp1)]) > 0)
            {
                mad.tmp1 <- mad.tmp1[-1]
                mad.tmp <- c(mad.tmp, mad.tmp1)
                madChrom[j,i] <- median(mad.tmp1)
            }
        }
        if (length(mad.tmp[!is.na(mad.tmp)]) > 0)
        {
            mad.tmp <- mad.tmp[-1]
            madGenome[i] <- median(mad.tmp)
        }
        
    }
    if (length(madGenome[is.na(madGenome)]) > 0)
        cat("Warning: MAD could not had ben computed for one of the\
samples\n")
    
    list(madChrom = madChrom, madGenome = madGenome)
    
}
#################################
#################################


findOutliers.func <-
    function(thres, factor=4, statesres)
{
    ##"state", "rpred", "prob", "pred", "disp", "obs"
    
    thres <- thres * factor
    
    chrom <- statesres[,1]
    sq.state <- seq(3, ncol(statesres), b = 6)
    sq.rpred <- seq(4, ncol(statesres), b = 6)
    sq.prob <- seq(5, ncol(statesres), b = 6)
    sq.pred <- seq(6, ncol(statesres), b = 6)
    sq.obs <- seq(8, ncol(statesres), b = 6)

    ##outputs
    ##1 = outlier	
    outlier <- matrix(0, nrow = length(chrom), ncol = length(sq.state))
    states.out <- matrix(0, nrow = length(chrom), ncol = length(sq.state))

    ##predicted values for all (including outliers) with median computed without outliers	
    pred.out <- matrix(0, nrow=length(chrom), ncol=length(sq.state))
    ##predicted values for all but outliers
    pred.obs.out <- matrix(0, nrow=length(chrom), ncol=length(sq.state))
    for (i in 1:length(sq.state))
    {
###        print(i)
        for (j in 1:length(unique(chrom)))
        {
            
            ind.nonna <- which(!is.na(statesres[chrom==j, sq.obs[i]]))
            states <- statesres[chrom == j, sq.state[i]][ind.nonna]
            obs <- statesres[chrom == j, sq.obs[i]][ind.nonna]
            pred <- statesres[chrom == j, sq.pred[i]][ind.nonna]
            pred.obs <- pred
            ##identify outliers
            for (k in 1:length(obs))
            {
                md <-  median(obs[states == states[k]],na.rm = TRUE)
                if ((obs[k] >  md + thres[i]) || (obs[k] < md - thres[i]))
                {
                    outlier[chrom==j, i][ind.nonna][k] <- 1
                    ##assigning observed value for a clone instead of predicted
                    pred.obs[k] <- obs[k]

#####NO longer do this ####################################################
                    ##state is -1 : will reassign state later as well as predicted value
                    ##possibly using pam() clustering and siluette width as a measure

                    ##states[k] <- -1
#####NO longer do this ####################################################
                }
            }
            ##recompute medians (predicted) without outliers -- use medians now	
            ##Note that predicted before indicated means and now indicate median

            ##states other than "-1" -- i.e. statesless
            
            states.uniq <- unique(states)
            for (m in 1:length(states.uniq))
            {
                
                ##predictions for all
                pred[states==states.uniq[m]] <-
                    median(obs[states == states.uniq[m] & outlier[chrom==j, i][ind.nonna] == 0])
                ##predictions for non-outliers only
                pred.obs[states == states.uniq[m] & outlier[chrom==j, i][ind.nonna] == 0] <-
                    median(obs[states==states.uniq[m] & outlier[chrom==j, i][ind.nonna] == 0])
                
            }
            
            pred.obs.out[chrom==j, i][ind.nonna] <- pred.obs
            pred.out[chrom==j, i][ind.nonna] <- pred

###############assign missing:
            outlier[chrom==j, i][-ind.nonna] <- NA
            pred.obs.out[chrom==j, i][-ind.nonna] <- NA
            pred.out[chrom==j, i][-ind.nonna] <- NA
            
        }
    }
    list(outlier=outlier, pred.obs.out=pred.obs.out, pred.out=pred.out)
    
}

#################################
#################################

findAber.func <-
    function(maxClones = 1, maxLen = 1000, statesres)
{
    ##either fewer than maxClones or length <= maxLen. 

    chrom <- statesres[,1]
    kb <- statesres[,2]
    sq.state <- seq(3, ncol(statesres), b=6)
    sq.obs <- seq(8, ncol(statesres), b=6)
    aber <- matrix(0, nrow=length(chrom), ncol=length(sq.state))
    for (i in 1:length(sq.state))
    {
###        print(i)
        for (j in 1:length(unique(chrom)))
        {
            
            st1 <- statesres[chrom==j, sq.obs[i]]
            ind.nonna <- which(!is.na(st1))
            states <- statesres[chrom==j, sq.state[i]][ind.nonna]
            kbnow <- kb[chrom==j][ind.nonna]
            
            abernow <- rep(0, length(kbnow))

            num <- 1
            for (m in 2:length(states))
            {
                if (states[m-1] != states[m])
                {
                    ##first clone is different from 2nd -> it's aberration	
                    if (m == 2)
                    {
                        abernow[1] <- 1
                    }
                    ##2nd clone is dif't from 3rd
                    if (m == 3)
                    {
                        abernow[1:2] <- 1
                    }
                    
                    ##last clone is different from previous -> it's aberration	
                    ##the clones before last may be an aberration as well
                    if (m == length(states))
                    {
                        
                        abernow[length(states)] <- 1
                        
                    }

                    ##clone before last is different from previous -> it's aberration	
                    if (m == (length(states)-1))
                    {
                        
                        abernow[(length(states)-1):(length(states))] <- 1
                        
                    }
                    
                    if (m <= length(states))
                    {	
                        ##take middle distances: if number of clones is small or they are very close together
                        
                        if ((num <= maxClones) || ((kbnow[m-1]-kbnow[m-num]) <= maxLen))
                            
                        {
                            abernow[(m-num):(m-1)] <- 1
                        }
                    }
                    num <- 1
                }
                else
                {
                    num <- num + 1
                }
            }
            aber[chrom==j, i][ind.nonna] <- abernow
            aber[chrom==j, i][-ind.nonna] <- NA
        }
    }
    list(aber=aber)
}

#################################
#################################

findTrans.func <-
    function(outliers, aber, statesres)
{

    ##exclude aberrations but keep outliers in

    chrom <- statesres[,1]
    kb <- statesres[,2]
    sq.state <- seq(3, ncol(statesres), b=6)
    sq.obs <- seq(8, ncol(statesres), b=6)
    ##transition matrix
    trans.matrix <- matrix(0, nrow=length(chrom), ncol=length(sq.state))

    ##lenght of the corresponding stretch matrix: 0 for aberrations and outliers

    translen.matrix <- matrix(0, nrow=length(chrom), ncol=length(sq.state))

    for (i in 1:length(sq.state))
    {
###        print(i)
        for (j in 1:length(unique(chrom)))
        {
            ind.nonna <- (1:length(statesres[chrom==j, sq.obs[i]]))[!is.na(statesres[chrom==j, sq.obs[i]])]
            kbnow <- kb[chrom==j][ind.nonna]
            states <- statesres[chrom==j, sq.state[i]][ind.nonna]
            outliersnow <- outliers[chrom==j,i][ind.nonna]
            abernow <- aber[chrom==j,i][ind.nonna]
            transnow <- rep(0, length(states))
            translennow <- rep(0, length(states))
            states.diff <- diff(states[abernow==0])
            ind <- (1:length(states.diff))[states.diff != 0]
            
            if (length(ind) > 0)
            {
                start <- ind+1
                end <-  ind
                
                transnow[abernow==0][start] <- 1
                transnow[abernow==0][end] <- 2
                
            }
            

#######
            ##compute the length of the corresponding stretches: same number is assigned for all clones
            ##between 1 and 2 including aberrations and outliers. distance to the first end is 
            ##computed from the start and of the last stretch is computed from the last clone to the
            ##end of chromosome.

            st <- c(1,(1:length(transnow))[transnow==1])
            en <- c((1:length(transnow))[transnow==2], length(transnow))
            
            for (m in 1:length(st))
            {
                translennow[st[m]:en[m]] <- kbnow[en[m]]-kbnow[st[m]]
                
            }
            
            translen.matrix[chrom==j,i][ind.nonna] <- translennow
            
            
############

            transnow[abernow==1] <- 3
            trans.matrix[chrom==j,i][ind.nonna] <- transnow

            trans.matrix[chrom==j, i][-ind.nonna] <- NA
            translen.matrix[chrom==j, i][-ind.nonna] <- NA
            
        }
    }
    list(trans.matrix=trans.matrix, translen.matrix=translen.matrix)

}

#################################
#################################

findAmplif.func <-
    	function(absValSingle = 1, absValRegion = 1.5, diffVal1=1,
    			diffVal2 = .5, maxSize =  10000, translen.matr, 
				trans.matr, aber, outliers, pred, pred.obs, statesres)
{
    chrom <- statesres[,1]
    kb <- statesres[,2]
    sq.state <- seq(3, ncol(statesres), b=6)
    sq.obs <- seq(8, ncol(statesres), b=6)

    amplif.matrix <- matrix(0, nrow=length(kb), ncol=length(sq.state))
    
    for (i in 1:length(sq.state))
    {
###        print(i)
        for (j in 1:length(unique(chrom)))
        {
            
            ind.nonna <- (1:length(statesres[chrom==j, sq.obs[i]]))[!is.na(statesres[chrom==j, sq.obs[i]])]
            
            abernow <- aber[chrom==j,i][ind.nonna]
            outliersnow <- outliers[chrom==j,i][ind.nonna]
            ##predicted value for the stretch

            prednow <- pred[chrom==j,i][ind.nonna]
            
            ##predicted value for the stretch, outliers have observed value
            predobsnow <- pred.obs[chrom==j,i][ind.nonna]
            
            obsnow <- statesres[chrom==j,sq.obs[i]][ind.nonna]
            transnow <- trans.matr[chrom==j,i][ind.nonna]
            translennow <- translen.matr[chrom==j,i][ind.nonna]

            amplifnow <- rep(0, length(obsnow))
            
########maybe remove############		
            ##if aberration or outlier and > absValSingle 		
            ##			
            ##			amplifnow[(abernow==1 | outliersnow ==1) & obsnow >= absValSingle] <- 1
##################################
            ##if aberration and greater by diffVal1 than max of the two surrounding  
            ##stretches or
            ##if outlier and greater by diffVal2 than its stretch and > 1
            
            ##outlier is much larger (diffVal) that its stretch
            amplifnow[outliersnow ==1 & ((obsnow - prednow)	>= diffVal1)] <- 1
            ##outlier and > 1 and diffVal2 greater than its stretch

            amplifnow[outliersnow ==1 & ((obsnow - prednow)	>= diffVal2) & obsnow >= absValSingle] <- 1

            ##aberration is much larger than maximum of the two surrounding stretches
            ##need to take spacial care when no stertches to the left or to the right 
            indaber <- (1:length(amplifnow))[abernow==1]
            if (length(indaber) > 0)
            {
                
                indstretch <- (1:length(amplifnow))[abernow==0 & outliersnow==0]
                for (m in 1:length(indaber))
                {
                    stretchleft <-
                        max(0,
                            max(indstretch[indstretch < indaber[m]]),
                            na.rm = TRUE)
                    stretchright <-
                        min((length(amplifnow)+1),
                            min(indstretch[indstretch > indaber[m]]),
                            na.rm = TRUE)
                    ##if no stretches to the left
                    if (stretchleft == 0)
                    {
                        mx <- prednow[stretchright]
                    } #if no stretches to the right:
                    else if (stretchright == (length(amplifnow)+1))
                    {
                        mx <- prednow[stretchleft]
                    }
                    else
                    {
                        mx <- max(prednow[stretchleft], prednow[stretchright])
                    }
                    if (!is.na(mx))
                    {
                        if (((predobsnow[indaber[m]] - mx) >= diffVal1) || ((predobsnow[indaber[m]] - mx) >= diffVal2 && (predobsnow[indaber[m]] >= absValSingle)))
                        {
                            amplifnow[indaber[m]] <- 1
                        }
                    }
                }	
            }			

            ##if part of the stretch and observed value of > absValRegion and 
            ##NOT YET: larger by diffValRegion than max of the surrounding regions regions and 
            ##size of the corresponding stretch <= maxSize
            
            amplifnow[abernow==0 & outliersnow ==0 & obsnow >= absValRegion & translennow <= maxSize] <- 1
            
            amplif.matrix[chrom==j,i][ind.nonna] <- amplifnow
            amplif.matrix[chrom==j, i][-ind.nonna] <- NA
        }
    }
    
    list(amplif = amplif.matrix)
    
}


######################################
######################################

plotChrom.hmm.func <-
    function(sample, chr, statesres, amplif, aber, outliers, trans,
    		pred, yScale = c(-2,2), maxChrom = 23, chrominfo,
    		samplenames, namePSfile = "try.ps", ps = TRUE, plotend = TRUE)
{

    chrom.rat <- chrominfo$length/max(chrominfo$length)
    chrom.start <- rep(0, maxChrom)
    for (i in 2:length(chrom.start))
    {
        chrom.start[i] <- sum(chrominfo$length[1:(i-1)])
    }
    ##
    ##
    ##chrom.mid contains middle positions of the chromosomes relative to
    ##the whole genome (useful for plotting the whole genome)
    ##
    chrom.mid <- rep(0, maxChrom)
    for (i in 1:length(chrom.start))
    {
        chrom.mid[i] <- chrom.start[i]+chrominfo$length[i]/2
    }
###############################################



    chrom <- statesres[,1]

    if (plotend)
    {
        if (ps)
        {
            postscript(namePSfile, paper="letter")
        }
        else
        {
            pdf(namePSfile, width=11, height=8.5)
        }
        par(mfrow=c(2,1))
    }
    par(lab=c(15,6,7), pch=18, cex=1, lwd=1)

    sq.state <- seq(3, ncol(statesres), b=6)
    sq.obs <- seq(8, ncol(statesres), b=6)

    for (j in 1:length(chr))
    {

        ind.nonna <- (1:length(statesres[chrom==chr[j], sq.obs[sample]]))[!is.na(statesres[chrom==chr[j], sq.obs[sample]])]

        kb <- (statesres[chrom==chr[j],2][ind.nonna])/1000
        obs <- statesres[chrom==chr[j], sq.obs[sample]][ind.nonna]
        states <- statesres[chrom==chr[j], sq.state[sample]][ind.nonna]
        nstates <- length(unique(states)) 

        abernow <- aber[chrom==chr[j],sample][ind.nonna]
        outliersnow <- outliers[chrom==chr[j],sample][ind.nonna]
        amplifnow <- amplif[chrom==chr[j],sample][ind.nonna]
        transnow <- trans[chrom==chr[j],sample][ind.nonna]

        ##predicted values when non-aberration of outlier: otherwise observed
        prednow <- obs
        prednow[outliersnow == 0 & abernow==0] <- (pred[chrom==chr[j],sample][ind.nonna])[outliersnow == 0 & abernow==0]

        y.min <- min(yScale[1], min(obs))
        y.max <- max(yScale[2], max(obs))

##################

        ##observed

        plot(kb, obs, xlab="", ylab="", ylim=c(y.min, y.max), type="l", col="blue", xlim=c(0, chrominfo$length[chr[j]]/1000))
        points(kb, obs,col="black")
        title(main=paste("Sample ", sample, " ", samplenames[sample], " - Chr ",chr[j], "Number of states ", nstates), xlab="kb (in 1000's)", ylab="data (observed)")
        abline(h=seq(y.min,y.max, b=.2), lty=3)
        abline(v=chrominfo$centromere[chr[j]]/1000, lty=2, col="red", lwd=3)
        ##start (dotted blue) and end of states (green)


        if (nstates > 1)
        {
            abline(v=kb[transnow==1], col="blue", lwd=2)
            abline(v=kb[transnow==2], col="green", lty=2, lwd=.5)
        }

###########
        ##amplif = red
        ##aber = green
        ##outliers = yellow


        if (length(outliersnow[outliersnow ==1]) > 0)
        {
            points(kb[outliersnow ==1], obs[outliersnow ==1], col="yellow")
        }
        if (length(abernow[abernow ==1]) > 0)
        {
            
            points(kb[abernow ==1], obs[abernow ==1], col="green")
        }
        if (length(amplifnow[amplifnow ==1]) > 0)
        {
            points(kb[amplifnow ==1], obs[amplifnow ==1], col="red")
        }


        
        ##predicted states:
        
        plot(kb, prednow, xlab="", ylab="", ylim=c(y.min, y.max), type="l", col="blue", xlim=c(0, chrominfo$length[chr[j]]/1000))
        points(kb, prednow,col="black")
        title(xlab="kb (in 1000's)", ylab="data (smoothed)")
        abline(h=seq(y.min,y.max, b=.2), lty=3)
        abline(v=chrominfo$centromere[chr[j]]/1000, lty=2, col="red", lwd=3)


        ##start (dotted blue) and end of states (green)
        if (nstates > 1)
        {
            abline(v=kb[transnow==1], col="blue", lwd=2)
            abline(v=kb[transnow==2], col="green", lty=2, lwd=.5)
        }



###########
        ##amplif = red
        ##aber = green
        ##outliers = yellow


        if (length(outliersnow[outliersnow ==1]) > 0)
        {
            points(kb[outliersnow ==1], obs[outliersnow ==1], col="yellow")
        }
        if (length(abernow[abernow ==1]) > 0)
        {
            
            points(kb[abernow ==1], obs[abernow ==1], col="green")
        }
        if (length(amplifnow[amplifnow ==1]) > 0)
        {
            points(kb[amplifnow ==1], obs[amplifnow ==1], col="red")
        }
    } 
    if (plotend)
    {	
	dev.off()
    }
}
##################################
##################################

plotCGH.hmm.func <-
    function (data, datainfo, chrominfo, samplename, sampNm,
    		  yScale = c(-2, 2), namePSfile = "try.ps", ps = TRUE,
    		  statesres, amplif, aber, outliers, trans)
{
################General Comments############################################

#########creating chromFull.info file

    chrom.uniq <- unique(datainfo$Chrom)

    chrominfo <- chrominfo[chrom.uniq,]

    sq.state <- seq(3, ncol(statesres), b=6)
    sq.obs <- seq(8, ncol(statesres), b=6)



#########
    chrom.rat <- chrominfo$length/max(chrominfo$length)  
    ##i.e. for each chromosome it repreesents the fraction of length of the
    ##longest chromosome
    ##
    ##chrom.start contains starting positions of the chromosomes relative to the
    ##whole genome (0 for the first)
    chrom.start <- rep(0, length(chrom.uniq))
    for (i in 2:length(chrom.start))
    {
	chrom.start[i] <- sum(chrominfo$length[1:(i-1)])
    }
    ##
    ##chrom.mid contains middle positions of the chromosomes relative to
    ##the whole genome (useful for plotting the whole genome)
    chrom.mid <- rep(0, length(chrom.uniq))
    for (i in 1:length(chrom.start))
    {
	chrom.mid[i] <- chrom.start[i]+chrominfo$length[i]/2
    }

    chromFull.info <- as.data.frame(cbind(chrominfo, chrom.start, chrom.mid, chrom.rat))
    dimnames(chromFull.info)[[2]] <- c("chr", "length", "centromere", "start", "mid", "rat")

########################################


    ##computing positions in genome for each clone:

    clone.genomepos <- rep(0, length(datainfo$kb))
    for (i in 1:length(chrom.uniq))
    {
	clone.genomepos[datainfo$Chrom==i] <- datainfo$kb[datainfo$Chrom==i]+chromFull.info$start[i]
    }

##########
    ##Now, determine vertical scale for each chromosome:

    y.min <- rep(yScale[1], length(chrom.uniq))
    y.max <- rep(yScale[2], length(chrom.uniq))

##############
    ##figure out the sample
    ##
    smpnames <- sampNm
    if ((samplename >= 1) && (samplename <= ncol(data)))
        ##samplename was the index
    {
	smp <- samplename
	samplename <- smpnames[smp]
        ##so now samplename is a name
    }
    else ##samplename 
    {
	
	smp <- (1:length(smpnames))[smpnames==samplename]
    }
##############
    ##values to plot:

    vals <- data[,smp]

#############
    ##adjust scales of chromosomes that have values outside a fixed scale

    for (i in 1:length(chrom.uniq))
    {
	y.min[i] <- min(c(vals[datainfo$Chrom==i],yScale[1]), na.rm = TRUE)
	y.max[i] <- max(c(vals[datainfo$Chrom==i],yScale[2]), na.rm = TRUE)
    }

    ##set genome scale to the max and min values across chrom's

    ygenome.min <- min(y.min, na.rm = TRUE)
    ygenome.max <- max(y.max, na.rm = TRUE)

#########################
    ##start a postscript file

    postscript(namePSfile, paper="letter", horizontal = FALSE)
    ##just a safety line
    close.screen(all = TRUE)
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
    j.seq <- 1:length(chrom.uniq)
    for (j in j.seq)
    {

        ind.nonna <- (1:length(vals[datainfo$Chrom==j]))[!is.na(vals[datainfo$Chrom==j])]
        screen(scr.seq[j])
        par(cex=.5, pch=20, lab=c(15,4,7), tcl=-.2, las=1, oma=c(0,0,0,0), cex.axis=1.3, cex.main=1.3, mgp=c(0,.15,0), lwd=.5)
        par(pin=c(chromFull.info$rat[j]*fact, .65))
        
        plot((datainfo$kb[datainfo$Chrom==j][ind.nonna])/1000, vals[datainfo$Chrom==j][ind.nonna], ylim=c(y.min[j],y.max[j]), xlab="", ylab="", type="l", col="blue", xlim=c(0, chromFull.info$length[j]/1000))
        points((datainfo$kb[datainfo$Chrom==j][ind.nonna])/1000, vals[datainfo$Chrom==j][ind.nonna], col="black")
        
        ##if (j < 23)
        ##{
        ##title(main=paste("Chr",j), line=.1)
        ##}
        ##else
        ##{
        ##	title(main="Chr. X", line=.1)
        ##}

        title(main=paste("Chr",j), line=.1)

        abline(h=seq(y.min[j],y.max[j], b=.5), lty=3)
        abline(v=0, lty=2)		
        abline(v=chromFull.info$centromere[j]/1000, lty=2, col="red")

####################
        ##plotting transitions and aberrations

        kb <- (datainfo$kb[datainfo$Chrom==j][ind.nonna])/1000
        chrom <- datainfo$Chrom	
        
        obs <- vals[chrom==j][ind.nonna]

        states <- statesres[chrom==j, sq.state[smp]][ind.nonna]
        nstates <- length(unique(states)) 

        abernow <- aber[chrom==j,smp][ind.nonna]
        outliersnow <- outliers[chrom==j,smp][ind.nonna]
        amplifnow <- amplif[chrom==j,smp][ind.nonna]
        transnow <- trans[chrom==j,smp][ind.nonna]


##################

        if (nstates > 1)
        {
            abline(v=kb[transnow==1], col="blue", lwd=1)
            abline(v=kb[transnow==2], col="green", lty=2, lwd=.25)
        }

###########
        ##amplif = red
        ##aber = orange
        ##outliers = yellow


        if (length(outliersnow[outliersnow ==1]) > 0)
        {
            points(kb[outliersnow ==1], obs[outliersnow ==1], col="yellow")
        }
        if (length(abernow[abernow ==1]) > 0)
        {
            points(kb[abernow ==1], obs[abernow ==1], col="orange")
        }
        if (length(amplifnow[amplifnow ==1]) > 0)
        {
            points(kb[amplifnow ==1], obs[amplifnow ==1], col="red")
        }
	



####################




        
    }
    
    ##plot genome:
    screen(9 )

    par(cex=.5, pch=20, lab=c(1,4,7), tcl=-.2, las=1, cex.axis=1.3, mgp=c(0,.15,0), cex.main=1.3, xaxs="i")
    par(pin=c(7.8, .55))
    plot(clone.genomepos/1000, vals, ylim=c(ygenome.min,ygenome.max), xlab="", ylab="", xlim=c(min(clone.genomepos[clone.genomepos>0], na.rm = TRUE)/1000, clone.genomepos[length(clone.genomepos[clone.genomepos>0])]/1000), col="black", type="l", lwd=1)
    title(main="Whole Genome (not to horizontal scale)",line=.1)
    for (i in seq(1,21,b=2))
    {	
        mtext(paste("", i), side = 1, at = (chromFull.info$mid[i]/1000), line=.3, col="red", cex.main=.5)
    }
    mtext("X", side = 1, at = (chromFull.info$mid[nrow(chromFull.info)]/1000), line=.3, col="red",cex.main=.5)
    abline(v=c(chromFull.info$start/1000, (chromFull.info$start[23]+chromFull.info$length[nrow(chromFull.info)])/1000), lty=1)
    abline(h=seq(ygenome.min,ygenome.max, b=.5), lty=3)
    abline(v=(chromFull.info$centromere+chromFull.info$start)/1000, lty=3, col="red")

    mtext(paste("Sample ", samplename, " ", smp, "Log2Ratio of Intensities vs Position in 1000's kb"), outer = TRUE, line=-1.2, cex=.8)
    dev.off()

####################
    

}

plotChrom.samples.func <-
    function(nr, nc, sample, chr, statesres, amplif, aber,
    		 outliers, trans, pred, yScale = c(-2, 2), 
			 maxChrom = 23, chrominfo = human.chrom.info.Jul03,
             samplenames)
{

    par(mfrow=c(nr,nc))

    chrom.rat <- chrominfo$length/max(chrominfo$length)
    chrom.start <- rep(0, maxChrom)
    for (i in 2:length(chrom.start))
    {
        chrom.start[i] <- sum(chrominfo$length[1:(i-1)])
    }
    ##
    ##
    ##chrom.mid contains middle positions of the chromosomes relative to
    ##the whole genome (useful for plotting the whole genome)
    ##
    chrom.mid <- rep(0, maxChrom)
    for (i in 1:length(chrom.start))
    {
        chrom.mid[i] <- chrom.start[i]+chrominfo$length[i]/2
    }
###############################################



    chrom <- statesres[,1]

    par(lab=c(15,6,7), pch=18, cex=1, lwd=1)

    sq.state <- seq(3, ncol(statesres), b=6)
    sq.obs <- seq(8, ncol(statesres), b=6)


    for (j in 1:length(chr))
    {

        ind.nonna <- (1:length(statesres[chrom==chr[j], sq.obs[sample[j]]]))[!is.na(statesres[chrom==chr[j], sq.obs[sample[j]]])]

        kb <- (statesres[chrom==chr[j],2][ind.nonna])/1000
        obs <- statesres[chrom==chr[j], sq.obs[sample[j]]][ind.nonna]
        states <- statesres[chrom==chr[j], sq.state[sample[j]]][ind.nonna]
        nstates <- length(unique(states)) 

        abernow <- aber[chrom==chr[j],sample[j]][ind.nonna]
        outliersnow <- outliers[chrom==chr[j],sample[j]][ind.nonna]
        amplifnow <- amplif[chrom==chr[j],sample[j]][ind.nonna]
        transnow <- trans[chrom==chr[j],sample[j]][ind.nonna]

        ##predicted values when non-aberration of outlier: otherwise observed
        prednow <- obs
        prednow[outliersnow == 0 & abernow==0] <- (pred[chrom==chr[j],sample[j]][ind.nonna])[outliersnow == 0 & abernow==0]


        y.min <- min(yScale[1], min(obs))
        y.max <- max(yScale[2], max(obs))

##################

        ##observed

        plot(kb, obs, xlab="", ylab="", ylim=c(y.min, y.max), type="l", col="blue", xlim=c(0, chrominfo$length[chr[j]]/1000))
        points(kb, obs,col="black")
        title(main=paste(samplenames[sample[j]], " - Chr ",chr[j], "Number of states ", nstates), xlab="kb (in 1000's)", ylab="data (observed)")
        ##abline(h=seq(y.min,y.max, b=.2), lty=3)
        abline(v=chrominfo$centromere[chr[j]]/1000, lty=2, col="red", lwd=3)
        ##start (dotted blue) and end of states (green)


        if (nstates > 1)
        {
            abline(v=kb[transnow==1], col="blue", lwd=2)
            abline(v=kb[transnow==2], col="green", lty=2, lwd=.5)
        }

###########
        ##amplif = red
        ##aber = orange
        ##outliers = yellow


        if (length(outliersnow[outliersnow ==1]) > 0)
        {
            points(kb[outliersnow ==1], obs[outliersnow ==1], col="yellow")
        }
        if (length(abernow[abernow ==1]) > 0)
        {
            points(kb[abernow ==1], obs[abernow ==1], col="orange")
        }
        if (length(amplifnow[amplifnow ==1]) > 0)
        {
            points(kb[amplifnow ==1], obs[amplifnow ==1], col="red")
        }

    } 

}

plotChrom.grey.samples.func <-
    function(nr, nc, sample, chr,  statesres, amplif, aber,
    		outliers, trans, pred, yScale = c(-2, 2), 
			maxChrom = 23, chrominfo = human.chrom.info.Jul03,
            samplenames)
{

    par(mfrow=c(nr,nc))

    chrom.rat <- chrominfo$length/max(chrominfo$length)
    chrom.start <- rep(0, maxChrom)
    for (i in 2:length(chrom.start))
    {
        chrom.start[i] <- sum(chrominfo$length[1:(i-1)])
    }
    ##
    ##
    ##chrom.mid contains middle positions of the chromosomes relative to
    ##the whole genome (useful for plotting the whole genome)
    ##
    chrom.mid <- rep(0, maxChrom)
    for (i in 1:length(chrom.start))
    {
        chrom.mid[i] <- chrom.start[i]+chrominfo$length[i]/2
    }
###############################################



    chrom <- statesres[,1]

    par(lab=c(15,6,7), pch=18, cex=1, lwd=1)

    sq.state <- seq(3, ncol(statesres), b=6)
    sq.obs <- seq(8, ncol(statesres), b=6)


    for (j in 1:length(chr))
    {

        ind.nonna <- (1:length(statesres[chrom==chr[j], sq.obs[sample[j]]]))[!is.na(statesres[chrom==chr[j], sq.obs[sample[j]]])]

        kb <- (statesres[chrom==chr[j],2][ind.nonna])/1000
        obs <- statesres[chrom==chr[j], sq.obs[sample[j]]][ind.nonna]
        states <- statesres[chrom==chr[j], sq.state[sample[j]]][ind.nonna]
        nstates <- length(unique(states)) 

        abernow <- aber[chrom==chr[j],sample[j]][ind.nonna]
        outliersnow <- outliers[chrom==chr[j],sample[j]][ind.nonna]
        amplifnow <- amplif[chrom==chr[j],sample[j]][ind.nonna]
        transnow <- trans[chrom==chr[j],sample[j]][ind.nonna]

        ##predicted values when non-aberration of outlier: otherwise observed
        prednow <- obs
        prednow[outliersnow == 0 & abernow==0] <- (pred[chrom==chr[j],sample[j]][ind.nonna])[outliersnow == 0 & abernow==0]


        y.min <- min(yScale[1], min(obs))
        y.max <- max(yScale[2], max(obs))

##################

        ##observed

        plot(kb, obs, xlab="", ylab="", ylim=c(y.min, y.max), type="l", col="grey50", xlim=c(0, chrominfo$length[chr[j]]/1000))
        points(kb, obs,col="black")
        title(main=paste(samplenames[sample[j]], " - Chr ",chr[j], "Number of states ", nstates), xlab="kb (in 1000's)", ylab="data (observed)")
        ##abline(h=seq(y.min,y.max, b=.2), lty=3)
        abline(v=chrominfo$centromere[chr[j]]/1000, lty=2, col="grey50", lwd=3)
        ##start (dotted blue) and end of states (green)


        if (nstates > 1)
        {
            abline(v=kb[transnow==1], col="black", lwd=2)
            abline(v=kb[transnow==2], col="black", lty=2, lwd=.5)
        }

###########
        ##amplif = red
        ##aber = orange
        ##outliers = yellow


        if (length(outliersnow[outliersnow ==1]) > 0)
        {
            points(kb[outliersnow ==1], obs[outliersnow ==1], col="grey80")
        }
        if (length(abernow[abernow ==1]) > 0)
        {
            points(kb[abernow ==1], obs[abernow ==1], col="grey50")
        }
        if (length(amplifnow[amplifnow ==1]) > 0)
        {
            points(kb[amplifnow ==1], obs[amplifnow ==1], col="grey30")
        }

    } 

}
