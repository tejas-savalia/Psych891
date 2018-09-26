kappa = 10
theta = 10
N = 1000
gammaSim <- function(N, kappa, theta){
  x <- rgamma(N, shape = kappa, scale = theta)
  hist(x, freq = FALSE, breaks = 30, col = "grey")
  
  #Plotting lines
  xfit = seq(min(gam), max(gam), length = 2000)
  yfit = dgamma(xfit, shape = kappa, scale = theta)
  
  lines(xfit, yfit)
  
  q=quantile(x,c(.1,.3,.5,.7,.9))
  abline(v=q)
  
  return (q)
}
#par(mfrow=(1, 1))
q = gammaSim(N, kappa, theta)
gammaPred <- function(par, q){

  #Parameters
  kappa = par[1]
  theta = par[2]
  
  #Cuts for quantiles
  cuts = c(0, q, Inf)
  p = c()
  
  for (i in (1:length(cuts)-1)){
    p[i] = pgamma(i+1, kappa, scale = theta) - pgamma(i, kappa, scale = theta)
  }
  return(p)
}
gammaPred(c(kappa, theta), q)

xfit = seq(min(gam), max(gam), length = 2000)
yfit = dgamma(xfit, shape = kappa, scale = theta)

lines(xfit, yfit)
library(e1071)

#Q: 3
mean = c()
median = c()
iqr = c()
range = c()
skew = c()

N = 10000
kappa = 1
for (i in c(1:10)){
  gam <- gammaSim(N, kappa, i)
  summary <- summary(gam)
  mean <- c(mean, summary[4])
  median <- c(median, summary[3])
  iqr <- c(iqr, summary[5] - summary[2])
  range <- c(range, summary[6] - summary[1])
  skew <- c(skew, skewness(gam))
}

mean
median
iqr
range
skew

#The mean and the median of the distribution increases as the scale parameter increases. 
#The spread, denoted by IQR and Range also increases
#The skewness of the distribution remains almost constant.


#Q 4)

mean = c()
median = c()
iqr = c()
range = c()
skew = c()

N = 10000
kappa = 1

for (i in c(1:10)){
  gam <- gammaSim(10000, 1, i)
  summary <- summary(gam)
  mean <- c(mean, summary[4])
  median <- c(median, summary[3])
  iqr <- c(iqr, summary[5] - summary[2])
  range <- c(range, summary[6] - summary[1])
  skew <- c(skew, skewness(gam))
}

mean
median
iqr
range
skew

#The effect on the distribution is similar to the previous one.

#Q 5)

mean = c()
median = c()
iqr = c()
range = c()
skew = c()


N = 10000
theta = 5
kappa_mul = c(1, 2, 3, 4, 5, 10, 20)

for (i in kappa_mul){
  gam <- gammaSim(10000, i, 5)
  summary <- summary(gam)
  mean <- c(mean, summary[4])
  median <- c(median, summary[3])
  iqr <- c(iqr, summary[5] - summary[2])
  range <- c(range, summary[6] - summary[1])
  skew <- c(skew, skewness(gam))
}

mean
median
iqr
range
skew


#The mean and median increase proportional to the increase in the shape parameter. 
#Same goes for IQR and range
#The skew decreases, the distribution shifting leftwards (But not left skewed yet)