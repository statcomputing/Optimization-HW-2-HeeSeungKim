
## problem 3

# (a)
## ------------example------------
x <- 1:10
y <- 2*x + 3                            # perfect fit
yeps <- y + rnorm(length(y), sd = 0.01) # added noise
nls(yeps ~ a + b*x, start = list(a = 0.12345, b = 0.54321))
## -------------------------------

graphics.off() # clear all previous plots
rm(list = ls()) # clear workspace
cat("\014") # clear console

x <- data.frame(
t = c(0,8,28,41,63,79,97,117,135,154),
x = c(2,47,192,256,768,896,1120,896,1184,1024))

temp <- rep(NA, nrow(x))


nls(x ~((x[1]*K)/(x[1]+(K-x[1])*exp(-r*t))), start=list(K=200,r=0.2),data=x,trace=TRUE)



