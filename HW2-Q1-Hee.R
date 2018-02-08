
graphics.off() # clear all previous plots
rm(list = ls()) # clear workspace
cat("\014") # clear console


## problem 1
# (b) Graph the log-likelihood function for a between -15 and 40
data <- c(1.77,-0.23,2.76,3.80,3.47,56.75,-1.34,4.24,-2.44,3.29,3.71,-2.40,4.53,-0.07,-1.05,-13.87,-2.53,-1.75);

a.value = seq(-15, 40, by=0.1);

m <- length(a.value)
y <- c()

for (i in 1:m) {
  y[i] <- -(length(data)*log(pi) + sum(log(1+(a.value[i]-data)^2))) #-log
}

plot(a.value, y, xlab="parameter", ylab="log-likelihood", type="l");


mle.cauchy_NR <- function(sample, int.value, lower=-10000, upper=10000) {
  n = length(sample);
  neg.log.lik <- function(a) {
    return(n*log(pi) + sum(log(1+(a-sample)^2)));
  }
  minus.score <- function(a) {
    return(2*sum((a-sample)/(1+(a-sample)^2)));
  }
  hess.f <- function(a) {
    return( matrix(2*sum( (1-(a-sample)^2) / (1+(a-sample)^2)^2 ),nrow=1));
  }
  ## MLE  
  a.est = nlminb(start = int.value, neg.log.lik, gradient=minus.score, hessian = hess.f, lower=lower,upper=upper);
  return(a.est);
}

St_points <- c(-11,-1,0,1.5,4,4.7,7,8,38)
NR1 <- mle.cauchy_NR(data,St_points[1])
NR2 <- mle.cauchy_NR(data,St_points[2])
NR3 <- mle.cauchy_NR(data,St_points[3])
NR4 <- mle.cauchy_NR(data,St_points[4])
NR5 <- mle.cauchy_NR(data,St_points[5])
NR6 <- mle.cauchy_NR(data,St_points[6])
NR7 <- mle.cauchy_NR(data,St_points[7])
NR8 <- mle.cauchy_NR(data,St_points[8])
NR9 <- mle.cauchy_NR(data,St_points[9])
NR10 <- mle.cauchy_NR(data,sum(data)/length(data))


NR_par <- c(NR1$par,NR2$par,NR3$par,NR4$par,NR5$par,NR6$par,NR7$par,NR8$par,NR9$par,NR10$par)
NR_obj <- -1*c(NR1$objective,NR2$objective,NR3$objective,NR4$objective,NR5$objective,NR6$objective,NR7$objective,NR8$objective,NR9$objective,NR10$objective)
NR_iter <- c(NR1$iterations,NR2$iterations,NR3$iterations,NR4$iterations,NR5$iterations,NR6$iterations,NR7$iterations,NR8$iterations,NR9$iterations,NR10$iterations)


NR <- rbind(NR_par,NR_obj,NR_iter)


# (c) 
St_points <- c(-11,-1,0,1.5,4,4.7,7,8,38,sum(data)/length(data))
data <- c(1.77,-0.23,2.76,3.80,3.47,56.75,-1.34,4.24,-2.44,3.29,3.71,-2.40,4.53,-0.07,-1.05,-13.87,-2.53,-1.75);
alpha <- c(1, 0.64, 0.25)
initial=-1;

Fu <- function(t){
  grad_l <- -2*sum((t-data)/(1+(t-data)^2))
  return(grad_l)
}

nlminb_1c <- function(x0,ap){
  X <- array()
  X[1] <- x0
  i=1
  difference <- 1
  while(abs(difference)>= 0.001 & i<100){
    X[i+1] <- X[i]+Fu(X[i])*ap
    difference <- X[i+1]-X[i]
    i <- i+1
  }
  return(X[i])
}

result <- array(dim=c(length(St_points),length(alpha)))

for (j in 1:length(St_points)){
  for (k in 1:length(alpha)){
  result[j,k] <- nlminb_1c(St_points[j],alpha[k])
  }
}
result


# (d) 
mle.cauchy_FS <- function(sample, int.value, lower=-10000, upper=10000) {
  n = length(sample);
  
  neg.log.lik <- function(a) {
    return(n*log(pi) + sum(log(1+(a-sample)^2)));
  }
  
  minus.score <- function(a) {
    return( 2*sum((a-sample)/(1+(a-sample)^2)));
  }
  
  fisher.inf <- function(a) {
    return(matrix(-n/2,nrow=1));
  }
  
  hess.f <- function(a) {
    return( matrix(2*sum( (1-(a-sample)^2) / (1+(a-sample)^2)^2 ),nrow=1));
  }

  ## MLE  
  a.est = nlminb(start = int.value, neg.log.lik, gradient=minus.score, hessian = fisher.inf, control = list(iter.max=3), lower=lower,upper=upper);
  b.est = nlminb(start = a.est$par, neg.log.lik, gradient=minus.score, hessian = hess.f, lower=lower,upper=upper);
  b.est$iterations <- a.est$iterations + b.est$iterations;
  return(b.est);
}

FS1 <- mle.cauchy_FS(data,St_points[1])
FS2 <- mle.cauchy_FS(data,St_points[2])
FS3 <- mle.cauchy_FS(data,St_points[3])
FS4 <- mle.cauchy_FS(data,St_points[4])
FS5 <- mle.cauchy_FS(data,St_points[5])
FS6 <- mle.cauchy_FS(data,St_points[6])
FS7 <- mle.cauchy_FS(data,St_points[7])
FS8 <- mle.cauchy_FS(data,St_points[8])
FS9 <- mle.cauchy_FS(data,St_points[9])
FS10 <- mle.cauchy_FS(data,sum(data)/length(data))

FS_par <- c(FS1$par,FS2$par,FS3$par,FS4$par,FS5$par,FS6$par,FS7$par,FS8$par,FS9$par,FS10$par)
FS_obj <- -1*c(FS1$objective,FS2$objective,FS3$objective,FS4$objective,FS5$objective,FS6$objective,FS7$objective,FS8$objective,FS9$objective,FS10$objective)
FS_iter <- c(FS1$iterations,FS2$iterations,FS3$iterations,FS4$iterations,FS5$iterations,FS6$iterations,FS7$iterations,FS8$iterations,FS9$iterations,FS10$iterations)

FS <- rbind(FS_par,FS_obj,FS_iter)

