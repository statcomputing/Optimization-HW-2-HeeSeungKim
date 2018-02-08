

graphics.off() # clear all previous plots
rm(list = ls()) # clear workspace
cat("\014") # clear console




## problem 2
data2 <- c(3.91,4.85,2.28,4.06,3.7,4.04,5.46,3.53,2.28,1.96,2.53,3.88,2.22,3.47,4.82,2.46,2.99,2.54,0.52);


## (a) Graph the log-likelihood function for a equally spaced numbers
a.value = seq(-pi, pi, length.out = 200);

m <- length(a.value)
y <- rep(0,m)

for (i in 1:m) {
  y[i] <- -length(data2)*log(2*pi) + sum(log(1-cos(data2-a.value[i])))
}

plot(a.value, y, xlab="parameter", ylab="log-likelihood", type="l");


## (b)
mme <- asin(sum(data2)/length(data2)-pi)


## (c)-(d)
mle.trig <- function(sample, int.value,  lower=-10000, upper=10000) {
  n = length(sample);
  
  neg.log.lik <- function(a) {
    return( n*log(2*pi) - sum(log(1-cos(sample-a))));
  }
  
  minus.score <- function(a) {
    return( sum (sin(sample-a) / (1-cos(sample-a))  ));
  }
  
  hess.f <- function(a) {
    return( matrix(sum (1/ (1-cos(sample-a))  ),nrow=1));
  }
  
  ## MLE  
  a.est = nlminb(start = int.value, neg.log.lik, gradient=minus.score, hessian = hess.f, lower=lower,upper=upper);
  
  
  return(a.est);
}


trig_MME <- mle.trig(data2,mme)
trig_Z1 <- mle.trig(data2,-2.7)
trig_Z2 <- mle.trig(data2,2.7)

trig_par <- c(trig_MME$par,trig_Z1$par,trig_Z2$par)
trig_obj <- -1*c(trig_MME$objective,trig_Z1$objective,trig_Z2$objective)
trig_iter <- c(trig_MME$iteration,trig_Z1$iteration,trig_Z2$iteration)

TRIG <- rbind(trig_par,trig_obj,trig_iter)
TRIG

 
## (e)
mle_trigs <- rep(0,200)
for (i in 1:200) {
  mle_trigs[i] <- mle.trig(data2,a.value[i])$par
}

plot(a.value,mle_trigs,xlab="Initial Value", ylab="MLE")

mle_trigs <- floor(mle_trigs*1000000)/1000000

SetAtt <- cbind(aggregate(a.value, list(mle_trigs), min),aggregate(a.value, list(mle_trigs), max)[,2])
colnames(SetAtt) <- c("Unique Outcome","From","To")
SetAtt
