mergeLevels <- function(vecObs,vecPred,pv.thres=0.0001,ansari.sign=0.05,thresMin=0.05,thresMax=0.5,verbose=1,scale=TRUE){

# Check if supplied thresholds are valid
if(thresMin>thresMax){cat("Error, thresMax should be equal to or larger than thresMin\n");return()}


# Initializing threshold and threshold vector for keeping track of thresholds
thresAbs=thresMin
sq<-numeric()

#initializing threshold index (threshold count)
j=0

#initializing ansari p-values to keep track of ansari p-values for each threshold in sq
ansari=numeric()

# Initialize levels count
lv=numeric()

# Set backtracking flag. Start with flag=0 indicating significance not yet reached, backtracking not begun
flag=0

# If thresMin=thresMax, fixed threshold is used and we set flag=2, only one run of the algoritm with initial thresMin
if(thresMin==thresMax){flag=2}

# Evaluate optimum steps for algorithm
else {
 l.step <- signif((thresMax-thresMin)/10,1)
 s.step <- signif((thresMax-thresMin)/200,1)
}

while (1){

  # Print current threshold if verbose is 1 or larger
  if(verbose>=1){cat("\nCurrent thresAbs: ",thresAbs,"\n")}

  j=j+1

  # Save current threshold
  sq[j]<-thresAbs

  # temporary predicted values (to be updated)
  vecPredNow=vecPred

  #unmissing unique segment medians
  mnNow=unique(vecPred)
  mnNow=mnNow[!is.na(mnNow)]

  #continuing indicator otherwise get out of the loop
  cont=0

  while(cont==0 & length(mnNow)>1) {

        mnNow=sort(mnNow)  #currennt sorted vector of means
        n <- length(mnNow)  # number of means in mnNow

        # Print current number of levels (n) if verbose is 2 or larger
        if(verbose>=2){ cat("\r",n,":",length(unique(vecPred)),"\t")}

        # Get distances translated to copy number differences
        # Only distances to closest levels
        if(scale){d<-(2*2^mnNow)[-n]-(2*2^mnNow)[-1]}
        else{d<-(mnNow)[-n]-(mnNow)[-1]}

        #order distance between means with the closest on top and corresponding indices
        dst<-cbind(abs(d)[order(abs(d))],(2:n)[order(abs(d))],(1:(n-1))[order(abs(d))])

        #for each pair of means
        for (i in 1:nrow(dst))  {
                #set continuity index to "NOT continue" (=1)
                cont=1
                #test for combining of the two segment means
                out=combine.func(diff=dst[i,1],vecObs, vecPredNow, mnNow, mn1=mnNow[dst[i,2]], mn2=mnNow[dst[i,3]], pv.thres=pv.thres, thresAbs=if(scale){2*2^thresAbs-2}else{thresAbs})
                #if combine?
                if (out$pv > pv.thres) {

                       #set continuity index to "YES" (=0) and break out of the current pairs loop
                       cont=0

                       #update predicted values and segments
                       vecPredNow=out$vecPredNow
                       mnNow=out$mnNow
                       break
                 }                
          }               
 }

### When done merging for a given threshold, test for significance ####
        ansari[j]=ansari.test(sort(vecObs-vecPredNow), sort(vecObs-vecPred))$p.value
  if(is.na(ansari[j])){ansari[j]=0} # If too many numbers for test to be performed, a 0 is returned, resulting in no merging (please use fixed threshold to get any merging)
  lv[j]=length(mnNow) # get number of levels

### If backtracking flag=2, the merging is stopped at this thresMax (or fixed threshold) ###
  if(flag==2){ break }

  # If p.value is less than the significance threshold, set backtracking flag=1 (backtracking on)
  if(ansari[j]<ansari.sign){
                        flag=1
  }

        
### If backtracking is on, a smaller threshold is attempted ####
        if (flag){

        # Stop if backtracking is on and p.value is higher than sign threshold or threshold is less or equal to thresMin
        if (ansari[j]>ansari.sign | thresAbs == thresMin){

#        # Don't merge at all if all tested threshold including thresMin is significant
#                         if (ansari[j] <= ansari.sign) {
#                                 vecPredNow=vecPred
#                                 mnNow=unique(vecPred)
#                                 mnNow=mnNow[!is.na(mnNow)]
#                         }
                                                  
        break
        }

      # Attempt smaller threshold
        else {
        thresAbs=signif(thresAbs-s.step,3)

        # Set threshold to thresMin as a minimum
        if (thresAbs <= thresMin){ thresAbs = thresMin }
      }
        }
        

### Increase threshold if backtracking is not on ###
        else {thresAbs=thresAbs+l.step}

#### Control step so function won't keep running, max threshold = thresMax and if sign not reached, threshold = thresMax ###
          if (thresAbs >= thresMax){
        thresAbs=thresMax
                    flag=2
          }

} # End while


# Return list of results
return(list(vecMerged=vecPredNow,mnNow=mnNow,sq=sq,ansari=ansari))
}

#################################


combine.func <- function(diff,vecObs, vecPredNow, mnNow, mn1, mn2, pv.thres=0.0001, thresAbs=0)
{ 
  #observed values in the first segment
        vec1=vecObs[which(vecPredNow==mn1)]
  #observed values in the second segment
        vec2=vecObs[which(vecPredNow==mn2)]
        
  #if difference between segment medians does not exceed thresAbs, then set pv=1
        if (diff<=thresAbs) {
                pv=1
        }
  #otherwise test for difference in mean based on observed values
        else {
                if((length(vec1) > 10 & length(vec2) > 10) | sum(length(vec1),length(vec2))>100){
                        pv=wilcox.test(vec1,vec2)$p.value
                }
                else{pv=wilcox.test(vec1,vec2,exact=T)$p.value  }       #/10^max(mn1,mn2)
                if(length(vec1) <= 3 | length(vec2) <= 3){pv=0}         
        }
        index.merged<-numeric()
  #if p-value exceeds pv.thres
        if (pv > pv.thres)      {
    #combine observed values
                vec=c(vec1,vec2)
    # Index values to be updated
                index.merged=which((vecPredNow==mn1) | (vecPredNow==mn2))               
    #update predicted values by median of the observed values
                vecPredNow[index.merged]=median(vec, na.rm=TRUE)
    #update segment medians  median of the observed values and remove one of the duplicates
                mnNow[which((mnNow==mn1) | (mnNow==mn2))]=median(vec, na.rm=TRUE)
                mnNow=unique(mnNow)
        }
        list(mnNow=mnNow, vecPredNow=vecPredNow, pv=pv)
}

#########################################