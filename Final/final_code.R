setwd("~/Courses/Psych 891/Psych891/Final")
rm(list = ls())

#The experiment is a modification of the classic SRT task. Participants learn a continuous
#sequence by moving the cursor on the corners of the screeen. This sequence will have a 
#pattern that the participants won't know when the experiment start.
#Classic SRT results suggests, participants learn the sequence implicitly; while they're
#not aware of the underlying sequence, decreasing RTs suggest they know the sequence. 

#Assumptions and Constraints
#Right now, all code is only for correct predictions for a sequence that has some pattern.
#This can then be compared with RTs when a sequence is truely random.
#Transition probabilities increase logarithmically. This will be eventually easier to fit
#using the RL model.
#Transition probabilities are collected from the data as the distance travelled from the 
#earlier stimulus towards the next correct one. 
#For the purpose of this submission, I've assumed that there is a definitive way to figure
#out the mapping between the distance travelled and boundary separation. That is, the 
#simulated data is my prediction of the boundary parameters

#On code
#Data is generated from a python file. Transposed in excel. Sorry, it was just easier for me
#that way. 

#Each trial consists of a 12 sub-trials as mapped on to 12 movements the participant has
#to make to go through 12 length sequence of stimuli. During a trial, the boundary parameter
#will be constant. Across trials, however, the boundary parameter changes with learning.
#Thus, DDMs are fit to each trial (so just 12 RTs) -- a little unrealisitc, perhaps.
#With a fixed starting parameter, all that needs to be optimized is the drift rate. 



#Data contains, in the second column, how far participants are on the screen
#from the next stimulus after they've started getting predictive. This corresponds to the
#boundary separation parameter. 

data = read.csv(file = "Predicted_start_points.csv")
head(data)

#total_rts = data[1] + data[2]
#plot(seq(1:100), total_rts[,1])

#quantile(total_rts[,1], c(0.1, 0.3, 0.5, 0.7, 0.9))
source("wald.R")
source("TruncNorm.R")

#RTs during each trial. Where a trial is a sequence of length 12. So there'd be 12 RTs. 
#Assume 100 trials
n_trials = 100
v = rtnorm(n_trials, mu = 0.3, sig = 0.1, lo = 0.1, hi = 0.5)

trial_rt = c()

#12 RTs generated from wald distribution for a given starting point and sampled drift rate.
#Quantiles are returned and put in a list. Unpacked in a matrix in the line that follows.
#
for (i in 1:n_trials){
  trial_rt = c(trial_rt, waldSimq(data[,2][n_trials - i + 1]/100, v[i], 12))
}

trial_rt_quantiles = matrix(trial_rt, nrow = n_trials, byrow = TRUE)
#Expected frequencies
f = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)*12
#Fit parameters to save
f_v = c()
f_wald_gsq = c()
cnt = 0
upfreq = 50


#cons=function(cpar){
#  penf=0
#  pmn=.01
#  if(cpar[1]<pmn){
#   penf=penf+(pmn-cpar[1]) 
#    cpar[1]=pmn
#  }
#  pmn=.1
#  if(cpar[2]<pmn){
#    penf=penf+(pmn-cpar[2]) 
#    cpar[2]=pmn
#  }
  #pmn=.01
  #if(cpar[3]<pmn){
  #  penf=penf+(pmn-cpar[3]) 
  #  cpar[3]=pmn
  #}
  #ls=list(cpar,penf)
 # return(ls)
#}

##############################################
#Functions taken from the Midterm code.

waldPred <- function(a, par, q){
  
  #Parameters
  #a = par[1]
  v = par[1]
  #print (a)
  #print (v)
  #Cuts for quantiles
  #print (q)
  cuts = c(0, q, Inf)
  p = c()
  for (i in (1:(length(cuts)-1))){
    p[i] = pwald(cuts[i+1], a, v) - pwald(cuts[i], a, v)
    #print (pwald(cuts[i], a, v))
  }
  return(p)
}

gsqfun_wald <- function(par, a, q, freq){
  
  #ls = cons(par)
  #cpar = ls[[1]]
  #pen = ls[[2]]
  #p = waldPred(a, cpar, q)
  p = waldPred(a, par, q)
  #print (p)
  p[p<.0001]=.0001
  p = p/sum(p)
  #print (p)
  
  pf = p*sum(freq)
  
  
  gslev=1:length(freq)
  gslev=gslev[freq>0]
  
  gsq=2*sum(freq[gslev]*log(freq[gslev]/pf[gslev]))
  #gsq = gsq + pen*gsq
  return (gsq)
}

#Since there's only one parameter to adjust (the drift rate), I'm doing a linear search
#between 0.01 and 1 with a stepsize of 0.01 and storing the value that gives minimum g-square
#error

linear_search = function(boundary, a, q, f){
  par = min(boundary)
  
  min_gsq = gsqfun_wald(par, a, q, f)
  #print ("here")
  for (i in seq(from = min(boundary), to = max(boundary), by = 0.01)){
    gsq = gsqfun_wald(i, a, q, f)
    #print (gsq)
    if (gsq < min_gsq){
      min_gsq = gsq
      par = i
    }
  }
  return(par)
  }

#linear_search(c(0.01, 0.5), a, q_wald, f_wald)
wald_fit = function(a, q, f){
  #par = c(0.11)
  
  #comp = gsqfun_wald(par, a, q, f)
#  while(TRUE){
#    fit = optimize(gsqfun_wald, interval = c(0, 10), a, q, f)
#    new_param = bisect(par, bound)
#    fit=optim(par,gsqfun_wald,a = a, freq=f,q=q)
#    par=fit$par
#    if((comp-fit$value)<=.001) break
#    comp=fit$value
#  }
  par = linear_search(c(0.01, 5), a, q, f)
  wald_par = par
  gsq = gsqfun_wald(a, wald_par, q, f)
  
  ls = list(wald_par, gsq)
  return (ls)
}

si= 1
cnt = 0
upfreq = 20

#Run loops to fit the drift rate for each trial (12 subtrials)
for(si in 1:n_trials){
  
  
  q_wald = trial_rt_quantiles[si,]
  a = data[,2][si]/100
  #runs the sim function above
  #q_wald=waldSimq(a[si], v[si], N)
  f_wald = f
  
  #Creates a vector of observed frequency
  #  counts in each RT bin separated by
  #  quantiles from the exgsim function
  #print (si)
  
  ls=wald_fit(a, q_wald, f_wald)
  #return fit parameters and gsq values. 
  
  fpar=ls[[1]]
  gval=ls[[2]]
  
#  f_a[si]=fpar[1]
  f_v[si]=fpar[1]
  
  f_wald_gsq[si]=gval
  
  
  
  
  #show progress updates
  if(cnt%%upfreq==0) print(cnt)
  cnt=cnt+1
}
f_wald_gsq
f_v
