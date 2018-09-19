kappa = 10
theta = 10
N = 1000
gammaSim <- function(N, kappa, theta){
  x <- rgamma(N, shape = kappa, scale = theta)
  return (x)
}
N
gam = gammaSim(N, kappa, theta)
hist(gam, freq = FALSE, breaks = 30)
xfit = seq(min(gam), max(gam), length = 2000)
yfit = dgamma(xfit, shape = kappa, scale = theta)

lines(xfit, yfit)
