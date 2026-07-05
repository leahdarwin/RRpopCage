#
cat("# Loaded: poolFreqDiffTest.R\n")
# Source the G-test
cat("# Looking in: ",currdir," for G_test.R\n")
source(paste(currdir,"/G_test.R",sep=""))
#
#
# FUNCTION: Woolf-test
# The script comes from the help page for the mantelhaen.test()
# ?mantelhaen.test()
woolf.test <- function(x) {
  x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  woolf <- sum(w * (log(or) - weighted.mean(log(or), w)) ^ 2)
  df <- k-1
  p <- 1 - pchisq(woolf, df)
  dat <- c(woolf,df,p)
  names(dat) <- c("Woolf", "df", "p-value")
  dat
}

# FUNCTION: Get GLM results from array ##
#get GLM data-set from k-way table
get_glm_dat <- function(array,zeroes=1){
  if(zeroes == 1){
    if(any(array == 0)){
      array <- array+1
    }
  }
  A_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  Tot_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  tr_l <- vector(length = dim(array)[3]*dim(array)[1])
  rep <- vector(length = dim(array)[3]*dim(array)[1])
  j <-1
  for(k in seq(1,dim(array)[3],1)){
    for(i in seq(1,dim(array)[1],1)){
      #      print(c(i,j,k))
      A_Cnt[j]<-array[i,1,k]
      Tot_Cnt[j]<-sum(array[i,,k])
      tr_l[j] <- as.character(i)
      rep[j] <- as.character(k)
      j <- j + 1
    }
  }
  d<-data.frame("A_Cnt"=A_Cnt,"Tot_Cnt"=Tot_Cnt,"tr_l"=tr_l,"rep"=rep)
  mod <- anova(glm(
    cbind(d$A_Cnt,d$Tot_Cnt-d$A_Cnt)~d$rep+d$tr_l+d$tr_l:d$rep,
    family = "binomial"),test="LRT")
  return(mod)
}
# FUNCTION: convert array to data.frame
get_dat <- function(array,zeroes=1){
  if(zeroes == 1){
    if(any(array == 0)){
      array<-array+1
    }
  }
  A_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  Tot_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  #tr_l <- vector(length = dim(array)[3]*dim(array)[1])
  #rep <- vector(length = dim(array)[3]*dim(array)[1])
  j <-1
  for(k in seq(1,dim(array)[3],1)){
    for(i in seq(1,dim(array)[1],1)){
      #      print(c(i,j,k))
      A_Cnt[j]<-array[i,1,k]
      Tot_Cnt[j]<-sum(array[i,,k])
      #tr_l[j] <- as.character(i)
      #rep[j] <- as.character(k)
      j <- j + 1
    }
  }


  ##CODE FOR REBECCAS'S RR POP CAGES
  ##-----------------------------------------------

  ##hard code treatment ordering based on /users/drand/data/RR_popcage_poolseq/aligned_reads_6.32/sync_files.txt
  ## {BEI,YAK,ZIM} ~ {1,2,3}
  ## mito = treatment = tr_l
  tr_l = c(rep(1,20),rep(2,20),rep(3,20))
  tr_l = as.character(tr_l)

  ##hard code replicate ordering
  ## {B1,...,B5,Y1,...,Y5,Z1,...Z5} ~ {1,...5,6,...,10,11,...15}
  rep = c(rep(1:15, each = 4))
  rep = as.character(rep)

  ##hard code time point ordering
  ## {F2,F10} ~ {2,10}
  time = c(rep(c(10,15,25,2), times = 15))
  #time = as.character(time)

  ##-----------------------------------------------

  d<-data.frame("A_Cnt"=A_Cnt,"Tot_Cnt"=Tot_Cnt,"tr_l"=tr_l,"rep"=rep, "time"=time)

  return(d)
}
