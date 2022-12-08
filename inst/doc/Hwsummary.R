## ----setup11, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
#Import data
data0 <- read.table("E:/Study/Graduate/Statistical_computing/Hw/Hw0/data.csv", sep=",", header=T)
summary(data0)

y <- as.matrix(data0[,1])#consumption
x <- as.matrix(data0[,2])#income

## ----setup12, fig.height=4, fig.width=10, echo=T, eval=T, results='asis'------
#Show the first six sets of data
knitr::kable(head(data0),align='c',caption="A Table of Consumption and Income")

Table1 <- xtable::xtable(head(data0))
print(Table1,type="html")

## ----setup13, fig.height=4, fig.width=10, echo=T, eval=T----------------------
#Draw the scatter plot of y(consumption)
plot(y,ylab="Consumption",main="The scatter plot of Consumption")
#Draw the scatter plot of x(income)
plot(x,ylab="Income",main="The scatter plot of Income")

## ----setup14, fig.height=4, fig.width=10, echo=T, eval=T----------------------
lm1 <- lm(y~x)
summary(lm1)
par(mfrow=c(2,2))
plot(lm1)

# Clean the memory of the variables
rm(list=ls())

## ----setup23, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
n <- 1000
u <- runif(n)#generate random numbers U~U(0,1)
x <- 2/(1-u)^{1/2}
#Graph the histogram of the random sample
hist(x, prob = TRUE, main="Pareto(2,2)", xlim=c(0,50))
box()
#Superimpose the theoretical density
y <- seq(2,100,0.1)
lines(y, 8/y^3, col="blue")

# Clean the memory of the variables
rm(list=ls())

## ----setup27, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
#The function to generate random numbers from Beta(a,b)
generate_beta = function(a,b){
  if(b==1){
    u <- runif(1)
    x <- u^(1/a) #Inverse transform method
  }else if(a==1&b!=1){
    u <- runif(1)
    x <- 1-u^(1/b) #Inverse transform method
  }else if(a>1&b>1){
    repeat{
      u1 <- runif(1)
      u2 <- runif(1)
      x <- u1
      a1 <- a-1;b1 <- b-1;ab <- a+b-2
      y <- u2*(gamma(a+b)*a1^a1*b1^b1)/(gamma(a)*gamma(b)*ab^ab)
      if(y<=dbeta(x,a,b))
        break
    }
  }else if(a<1&b>1){
    repeat{
      u1 <- runif(1)
      y <- runif(1)
      x <- u1^(1/a)
      if(y<=(1-x)^(b-1))
        break
    }
  }else if(a>1&b<1){
    repeat{
      u1 <- runif(1)
      y <- runif(1)
      x <- 1-u1^(1/b)
      if(y<=x^(a-1))
        break
    }
  }else{
    repeat{
      u1 <- runif(1)
      y <- runif(1)
      if(u1<=0.5){
        x <- (u1/2^(a+b-2))^(1/(a+b-1))
        if(y<=(1/x-1)^(b-1))
          break
      }else{
        x <- 1-((1-u1)/2^(a+b-2))^(1/(a+b-1))
        if(y<=(1/x-1)^(1-a))
          break
      }
    }
  }
  return(x)
}

## ----setup271, fig.height=4, fig.width=10, echo=T, eval=T---------------------
n <- 1000
xbeta <- numeric(n)
for(i in 1:n){
  xbeta[i] <- generate_beta(3,2)
}
#Graph the histogram of the random sample
hist(xbeta, prob=TRUE, main="Beta(3,2)")
#Superimpose the theoretical Beta(3,2) density
curve(dbeta(x,3,2),col="blue",add=TRUE)

# Clean the memory of the variables
rm(list=ls())

## ----setup212, fig.height=4, fig.width=10, echo=T, eval=T---------------------
set.seed(22058)
n <- 1000;r <- 4;beta <- 2
lambda <- rgamma(n,r,beta)
x <- rexp(n,lambda)
#Graph the histogram of the random sample
hist(x, prob = TRUE, main="The histogram of the random sample from Exponential-Gamma mixture")
box()

# Clean the memory of the variables
rm(list=ls())

## ----setup213, fig.height=4, fig.width=10, echo=T, eval=T---------------------
set.seed(22058)
n <- 1000;r <- 4;beta <- 2
lambda <- rgamma(n,r,beta)
x <- rexp(n,lambda)
#Graph the histogram of the random sample
hist(x, prob = TRUE, main="The histogram of the random sample with the Pareto density curve superimposed")
y <- seq(0,100,0.1)
#Superimpose the theoretical density
lines(y, 64/(2+y)^5, col="blue")
box()

# Clean the memory of the variables
rm(list=ls())

## ----setup31, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
# The fast sorting function
fast_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(fast_sort(lower),a,fast_sort(upper)))}
}

# Apply the fast sorting function
n1=1e4;n2=2e4;n3=4e4;n4=6e4;n5=8e4
test1<-sample(1:n1)
test2<-sample(1:n2)
test3<-sample(1:n3)
test4<-sample(1:n4)
test5<-sample(1:n5)

# 100 simulations
y1 <- y2 <- y3 <- y4 <- y5 <- numeric(100)
for(i in 1:100){
  y1[i] <- system.time(fast_sort(test1))[1]
  y2[i] <- system.time(fast_sort(test2))[1]
  y3[i] <- system.time(fast_sort(test3))[1]
  y4[i] <- system.time(fast_sort(test4))[1]
  y5[i] <- system.time(fast_sort(test5))[1]
}
# Calculate the average of computation time
a1 <- mean(y1)
a2 <- mean(y2)
a3 <- mean(y3)
a4 <- mean(y4)
a5 <- mean(y5)
a <- c(a1,a2,a3,a4,a5)
# Calculate the theoretical computation time
t1 <- n1*log(n1)
t2 <- n2*log(n2)
t3 <- n3*log(n3)
t4 <- n4*log(n4)
t5 <- n5*log(n5)
t <- c(t1,t2,t3,t4,t5)

# Regression about Simulated computation time and Theoretical computation time
lm1 <- lm(a~t)
summary(lm1)

## ----setup311, fig.height=4, fig.width=10, echo=T, eval=T---------------------
# Draw the scatter plot
plot(t,a,main="Regression about Simulated computation time and Theoretical computation time",xlab="Theoretical computation time t=nlogn",ylab="Simulated computation time")
# Superimpose the regression line
abline(lm1,col="red")

# Clean the memory of the variables
rm(list=ls())

## ----setup361, fig.height=4, fig.width=10, echo=T, eval=T---------------------
set.seed(22058)
# The theoretical value of Cov(e^U,e^(1-U))
cov1 <- -exp(2)+3*exp(1)-1
cat("The theoretical covariance is",cov1,"\n")
# The theoretical value of Var(e^U+e^(1-U))
var1 <- -3*exp(2)+10*exp(1)-5
cat("The theoretical variance is",var1,"\n")

## ----setup362, fig.height=4, fig.width=10, echo=T, eval=T---------------------
N <- 10000
U1 <- runif(N/2)
y11 <- exp(U1)
y12 <- exp(1-U1)
cat("The empirical estimate of covariance is",cov(y11,y12))
y1 <- y11+y12
cat("The empirical estimate of variance",var(y1))

## ----setup363, fig.height=4, fig.width=10, echo=T, eval=T---------------------
# The theoretical percent reduction in variance
reducper <- (-6*exp(1)+2*exp(2)+2)/(-exp(2)+4*exp(1)-3)
cat("The theoretical percent reduction in variance is",reducper,"\n")

# Clean the memory of the variables
rm(list=ls())

## ----setup37, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
# The antithetic variate method
antithetic <- function(N){
  U <- runif(N/2)
  y <- c(exp(U),exp(1-U))
  return(mean(y))
}
# The simple MC method
simpleMC <- function(N){
  U <- runif(N)
  y <- exp(U)
  return(mean(y))
}

N <- 10000
m <- 100
# m simulations
Est1 <- Est2 <- numeric(m)
for(i in 1:m){
  Est1[i] <- antithetic(N)
  Est2[i] <- simpleMC(N)
}

cat(" The antithetic variable estimator is",mean(Est1),"\n","The simple MC estimator is",mean(Est2),"\n")


## ----setup371, fig.height=4, fig.width=10, echo=T, eval=T---------------------
# The variance of the antithetic variate method
var1 <- var(Est1)
# The variance of the simple MC
var2 <- var(Est2)
# The empirical estimate of the percent reduction in variance
cat("The empirical estimate of the percent reduction in variance is",(var2-var1)/var2,"\n")

# Clean the memory of the variables
rm(list=ls())

## ----setup41, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
x <- seq(1,5,0.01)
g <- x^2*exp(-x^2/2)/(2*pi)^0.5
f1 <- (2/pi)^0.5*exp(-(x-1)^2/2)
f2 <- x*exp((1-x^2)/2)
plot(x,g,ylim=c(0,1),lty=1,type="l",ylab="",main="Figure1: The curves of g(x), f1(x) and f2(x)")
lines(x,f1,col="red",lty=1)
lines(x,f2,col="blue",lty=1)
legend("topright",c("g(x)","f1(x)","f2(x)"),lty=c(1,1,1),col=c("black","red","blue"))

## ----setup411, fig.height=4, fig.width=10, echo=T, eval=T---------------------
set.seed(22058)
# The function of g(x)
g <- function(x){
  return(x^2*exp(-x^2/2)/(2*pi)^0.5)
}

# The importance function f1
f1 <- function(x){
  return((2/pi)^0.5*exp(-(x-1)^2/2))
}
# The function to generate random numbers from f1
f1sample <- function(N){
  x <- rnorm(N)
  z <- abs(x)+1
  return(z)
}
# The importance sampling function of f1
importancef1 <- function(N){
  x <- f1sample(N)
  y <- g(x)/f1(x)
  return(mean(y))
}

# The importance function f2
f2 <- function(x){
  return(x*exp((1-x^2)/2))
}
# The function to generate random numbers from f2
f2sample <- function(N){
  U <- runif(N,0,1)
  x <- (1-2*log(U))^0.5
  return(x)
}
# The importance sampling function of f2
importancef2 <- function(N){
  x <- f2sample(N)
  y <- g(x)/f2(x)
  return(mean(y))
}

n <- 100 # The number of simulations
N <- 1000 # The size of random sample
y1 <- y2 <- numeric(n)
for(i in 1:n){
  y1[i] <- importancef1(N)
  y2[i] <- importancef2(N)
}
cat(" The estimator using f1 is",mean(y1),";","\n",
    "The variance is",var(y1),".")
cat(" The estimator using f2 is",mean(y2),";","\n",
    "The variance is",var(y2),".")

# Clean the memory of the variables
rm(list=ls())

## ----setup42, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
M <- 1000 # The number of replicates
k <- 5 # The number of stratum
r <- M/k # Replicates per stratum
n <- 10 # The number of simulation

g <- function(x){
  return(exp(-x)/(1+x^2))
}
# The importance function for each stratum
fj <- function(x,j){
  return(exp(j/5-x)/(1-exp(-1/5)))
}
# The function to generate random numbers from fj
fjsample <- function(N,j){
  U <- runif(N,0,1)
  x <- j/5-log(1-(1-exp(-1/5))*U)
  return(x)
}
# The importance sampling for the jth subinterval
importancefj <- function(N,j){
  x <- fjsample(N,j)
  y <- g(x)/fj(x,j)
  return(mean(y))
}

# The stratified importance sampling
Est <- numeric(n)
Estj <- numeric(k)
for(i in 1:n){
  for (j in 1:k){
    Estj[j] <- importancefj(r,j-1)
  }
  Est[i] <- sum(Estj)
}

cat(" The estimator is",mean(Est),";","\n",
    "The standard deviation is",sd(Est),".")

# Clean the memory of the variables
rm(list=ls())

## ----setup51, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
# Data generation
Generadata <- function(n){
  return(rnorm(n,1,2))
}

# Data analysis
# Calculate the 95% confidence interval
CL <- function(n,alpha){
  y <- Generadata(n)
  LCL <- mean(y)-sd(y)*qt(1-alpha/2,n-1)/(n)^0.5
  UCL <- mean(y)+sd(y)*qt(1-alpha/2,n-1)/(n)^0.5
  return(c(LCL,UCL))
}

# Result reporting
m <- 1000 # The number of simulation
n <- 50 # The size of random sample
alpha <- 0.05
CLvalue <- matrix(0,m,2) # Storage for 95% confidence interval
for(i in 1:m){
  CLvalue[i,] <- CL(n,alpha)
}
# Count the number of intervals that contain μ=1
count=0
for(i in 1:m){
  if(CLvalue[i,1]<1&&CLvalue[i,2]>1)
    count=count+1
}
# Calculate the empirical estimate of the confidence level
cat("The empirical estimate of the confidence level is",count/m)

# Clean the memory of the variables
rm(list=ls())

## ----setup52, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
# Data generation
# Generate samples under H1
sigma1 <- 1
sigma2 <- 1.5
m <- 1000 # The number of simulation
n <- c(20,200,2000) # The size of random sample

# Data analysis
# The function of Count Five test 
count5test <- function(x, y) {
X <- x-mean(x)
Y <- y-mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(as.integer(max(c(outx, outy))>5))
}

# Result reporting
count5power <- Fpower <- numeric(3)
for(i in 1:3){
  count5power[i] <- mean(replicate(m, expr={
    x <- rnorm(n[i], 0, sigma1)
    y <- rnorm(n[i], 0, sigma2)
    count5test(x, y)
    }))
}

for(i in 1:3){
  Fpvalues <- replicate(m, expr={
    x <- rnorm(n[i], 0, sigma1)
    y <- rnorm(n[i], 0, sigma2)
    Ftest <- var.test(x,y,1-0.055)
    Ftest$p.value
    })
  Fpower[i] <- mean(Fpvalues<=0.055)
}

Countdata <- data.frame(n,count5power)
knitr::kable(Countdata,align='c',
             caption="The powers of Count Five test for small, medium, and large sample sizes",
             col.names=c("Sample size","The power of Count Five test"))

Fdata <- data.frame(n,Fpower)
knitr::kable(Fdata,align='c',
             caption="The powers of F test for small, medium, and large sample sizes",
             col.names=c("Sample size","The power of F test"))

# Clean the memory of the variables
rm(list=ls())

## ----setup53,echo=F-----------------------------------------------------------
A <- matrix(c("","0","1","0","a","b","1","c","d"),3,byrow=T)
colnames(A) <- c("","Method1","")
rownames(A) <- c("","Method2","")
knitr::kable(A,align='c')

# Clean the memory of the variables
rm(list=ls())

## ----setup61, fig.height=4, fig.width=10, echo=T, eval=T----------------------
# Data generation
set.seed(22058)
data1 <- c(3,5,7,18,43,85,91,98,100,130,230,487)

# Data analysis
# Bootstrap
library(boot);library(MASS)
boot.lambda <- function(x,i) 1/mean(x[i])
obj <- boot(data=data1,statistic=boot.lambda,R=1e4)

# Result reporting
round(c(original=obj$t0,bias=mean(obj$t)-obj$t0,se=sd(obj$t)),5)

# Clean the memory of the variables
rm(list=ls())

## ----setup62, fig.height=4, fig.width=10, echo=T, eval=T----------------------
# Data generation
set.seed(22058)
data1 <- c(3,5,7,18,43,85,91,98,100,130,230,487)
ci.norm<-ci.basic<-ci.perc<-ci.bca<-numeric(2)

# Data analysis
library(boot)
boot.mean <- function(x,i) mean(x[i])
obj1 <- boot(data=data1,statistic=boot.mean,R=1000)
ci <- boot.ci(obj1,type=c("norm","basic","perc","bca"))
ci.norm<-ci$norm[2:3]
ci.basic<-ci$basic[4:5]
ci.perc<-ci$percent[4:5]
ci.bca<-ci$bca[4:5]

## ----setup621, fig.height=4, fig.width=10, echo=T, eval=T---------------------
# Result reporting
cat(' standard normal:','[',ci.norm,']','\n',
    'basic:','[',ci.basic,']','\n',
    'percentile:','[',ci.perc,']','\n',
    'BCa:','[',ci.bca,']','\n')

# Clean the memory of the variables
rm(list=ls())

## ----setup63, fig.height=4, fig.width=10, echo=T, eval=T----------------------
# Data generation
n <- 20 # The size of random sample
m <- 1000 # The number of simulations

# Data analysis
set.seed(22058)
library(boot)
boot.mean <- function(x,i) mean(x[i])
ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)
for(i in 1:m){
  x <- rnorm(n)
  obj2 <- boot(data=x,statistic=boot.mean,R=1000)
  ci <- boot.ci(obj2,type=c("norm","basic","perc"))
  ci.norm[i,]<-ci$norm[2:3]
  ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]
}

## ----setup631, fig.height=4, fig.width=10, echo=T, eval=T---------------------
# Result reporting
mu <- 0
# The empirical coverage rates
cat('norm =',mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),
    'basic =',mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),
    'perc =',mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu))

## ----setup632, fig.height=4, fig.width=10, echo=T, eval=T---------------------
cat('norm =',mean(ci.norm[,1]>mu),
    'basic =',mean(ci.basic[,1]>mu),
    'perc =',mean(ci.perc[,1]>mu))

## ----setup633, fig.height=4, fig.width=10, echo=T, eval=T---------------------
cat('norm =',mean(ci.norm[,2]<mu),
    'basic =',mean(ci.basic[,2]<mu),
    'perc =',mean(ci.perc[,2]<mu))

# Clean the memory of the variables
rm(list=ls())

## ----setup71, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
library(bootstrap)

# Data analysis
# The function to compute the sample estimate of theta
thetaEst <- function(x){
  lambda.hat <- eigen(cov(x))$values
  return(lambda.hat[1]/sum(lambda.hat))
}
theta.hat <- thetaEst(scor) # The original statistic

# Jackknife
n <- nrow(scor)
theta.jack <- numeric(n)
for(i in 1:n){
  # The Jackknife statistic
  theta.jack[i] <- thetaEst(scor[-i,])
}
# The jackknife estimate of bias
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
# The jackknife estimates of standard error
se.jack <- sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2))

# Result reporting
round(c(original=theta.hat,
        bias.jack=bias.jack,
        se.jack=se.jack),3)

# Clean the memory of the variables
rm(list=ls())

## ----setup72, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
library(DAAG)
attach(ironslag)
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- matrix(NA,n*(n-1)/2,2)

# The function of four models
Jm <- function(J,x1,y1,s){
  if(s==1){# The first model
    yhat <- J1$coef[1] + J1$coef[2]*x1}
  else if(s==2){# The second model
    yhat <- J2$coef[1] + J2$coef[2]*x1 +J2$coef[3]*x1^2}
  else if(s==3){# The third model
    logyhat3 <- J3$coef[1] + J3$coef[2]*x1
    yhat <- exp(logyhat3)}
  else if(s==4){# The fourth model
    logyhat4 <- J4$coef[1] + J4$coef[2]*log(x1)
    yhat <- exp(logyhat4)}
  return(y1-yhat)
}

# Fit models on leave-two-out samples
k <- 1
for(i in 1:(n-1)){
  for(j in (i+1):n){
    y <- magnetic[-c(i,j)];x <- chemical[-c(i,j)]
    
    # The first model
    J1 <- lm(y ~ x)
    e1[k,] <- Jm(J1,chemical[c(i,j)],magnetic[c(i,j)],1)
    # The Second model
    J2 <- lm(y ~ x + I(x^2))
    e2[k,] <- Jm(J2,chemical[c(i,j)],magnetic[c(i,j)],2)
    # The third model
    J3 <- lm(log(y) ~ x)
    e3[k,] <- Jm(J3,chemical[c(i,j)],magnetic[c(i,j)],3)
    # The fourth model
    J4 <- lm(log(y) ~ log(x))
    e4[k,] <- Jm(J4,chemical[c(i,j)],magnetic[c(i,j)],4)
    
    k <- k+1
  }
}

# Result reporting
c(Model1=mean(e1^2), Model2=mean(e2^2), Model3=mean(e3^2), Model4=mean(e4^2))

## ----setup721, fig.height=4, fig.width=10, echo=T, eval=T---------------------
M2 <- lm(magnetic ~ chemical + I(chemical^2))
summary(M2)

## ----setup722, fig.height=4, fig.width=10, echo=T, eval=T---------------------
detach(ironslag)
# Clean the memory of the variables
rm(list=ls())

## ----setup73, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)
# Data generation
attach(chickwts)
x <- as.vector(weight[feed == "sunflower"])
y <- as.vector(weight[feed == "linseed"])
detach(chickwts)

# Data analysis
R <- 999
K1 <- length(x);K2 <- length(y)
cors <- numeric(R)
# The permutation
cor0 <- cor(x,y,method="spearman")# The original statistic
for (i in 1:R) {
  k1 <- sample(1:K1, size = K1, replace = FALSE)
  k2 <- sample(1:K2, size = K2, replace = FALSE)
  x1 <- x[k1];y1 <- y[k2]
  # The permutation statistic
  cors[i] <- cor(x1,y1,method="spearman")
}
# The significance level of permutation
p <- mean(abs(c(cor0,cors))>=abs(cor0))

# Result reporting
c(permutation=p,
  cortest=cor.test(x,y,method="spearman")$p.value)

# Clean the memory of the variables
rm(list=ls())

## ----setup81, fig.height=10, fig.width=10, echo=T, eval=T---------------------
set.seed(22058)
# Target function: standard Laplace distribution
target <- function(x){
  return(0.5*exp(-abs(x)))
}

# Random walk Metropolis sampler
# Using a normal distribution as a proposal distribution
rw.Metropolis <- function(N,sigma,x0) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N){
    y <- rnorm(1, x[i-1], sigma)
    A <- target(y)/target(x[i-1])
    if(u[i] <= A){
      x[i] <- y # Accept the new value
      k <- k+1
    }else{
      x[i] <- x[i-1] # Keep the old value
    }
  }
  return(list(x=x, k=k/N))
}

N <- 10000 # Length of chain
b <- 1000 # The Burn-in length
k <- 4 # The number of chain
sigma <- c(0.05,0.5,2,16)
x0 <- c(5,10,25,35)
rw1 <- rw.Metropolis(N,sigma[1],x0[3])
rw2 <- rw.Metropolis(N,sigma[2],x0[3])
rw3 <- rw.Metropolis(N,sigma[3],x0[3])
rw4 <- rw.Metropolis(N,sigma[4],x0[3])

# Result reporting
# The acceptance rates of each chain
print(c(sigma0.05=rw1$k,sigma0.5=rw2$k,sigma2=rw3$k,sigma16=rw4$k))

## ----setup812, fig.height=10, fig.width=10, echo=T, eval=T--------------------
# The Gelman-Rubin method
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi);k <- nrow(psi)
  psi.means <- rowMeans(psi)
  B <- n*var(psi.means) # between variance est.
  psi.w <- apply(psi,1,"var") # within variances
  W <- mean(psi.w) # within est.
  v.hat <- W*(n-1)/n+(B/n) # upper variance est.
  r.hat <- v.hat/W # G-R statistic
  return(r.hat)
}

# Compute diagnostic statistics
chains1 <- chains2 <- matrix(0,k,N)
for(i in 1:k){
  chains1[i,] <- rw.Metropolis(N,sigma[2],x0[i])$x
  chains2[i,] <- rw.Metropolis(N,sigma[3],x0[i])$x
}
psi1 <- t(apply(chains1, 1, cumsum))
psi2 <- t(apply(chains2, 1, cumsum))
for(i in 1:nrow(psi1)){
  psi1[i,] <- psi1[i,]/(1:ncol(psi1))
}
for(i in 1:nrow(psi2)){
  psi2[i,] <- psi2[i,]/(1:ncol(psi2))
}
print(c(sigma0.5=Gelman.Rubin(psi1),sigma2=Gelman.Rubin(psi2)))

# Plot psi for the four chains
par(mfrow=c(2,2))
# Sigma=0.5
for(i in 1:k){
  plot(psi1[i,(b+1):N],type="l",xlab=i,ylab=bquote(psi))
}
# Sigma=2
for(i in 1:k){
  plot(psi2[i,(b+1):N],type="l",xlab=i,ylab=bquote(psi))
}
par(mfrow=c(1,1)) #restore default

## ----setup813, fig.height=4, fig.width=10, echo=T, eval=T---------------------
#plot the sequence of R-hat statistics
rhat1 <- rhat2 <- rep(0,N)
for (j in (b+1):N){
  rhat1[j] <- Gelman.Rubin(psi1[,1:j])
  rhat2[j] <- Gelman.Rubin(psi2[,1:j])
}
par(mfrow=c(1,2))
# Sigma=0.5
plot(rhat1[(b+1):N], type="l", main="sigma=0.5", ylab="R")
abline(h=1.2, lty=2, col="red")
# Sigma=2
plot(rhat2[(b+1):N], type="l", main="sigma=2", ylab="R")
abline(h=1.2, lty=2, col="red")
par(mfrow=c(1,1))

# Clean the memory of the variables
rm(list=ls())

## ----setup82, fig.height=4, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)

# Function to generate the chain
bivarchain <- function(X,N){
  mu1 <- mu2 <- 0;sigma1 <- sigma2 <- 1;rho <- 0.9
  for (i in 2:N){
    x2 <- X[i-1,2]
    m1 <- mu1+rho*(x2-mu2)*sigma1/sigma2
    s1 <- sqrt(1-rho^2)*sigma1
    X[i,1] <- rnorm(1,m1,s1)
    x1 <- X[i,1]
    m2 <- mu2+rho*(x1-mu1)*sigma2/sigma1
    s2 <- sqrt(1-rho^2)*sigma2
    X[i,2] <- rnorm(1,m2,s2)
  }
  return(X)
}

# Initialize constants and parameters
N <- 5000 # Length of chain
burn <- 1000 # Burn-in length
X <- matrix(0,N,2) # The chain, a bivariate sample

mu1 <- mu2 <- 0 # The initial value
X[1,] <- c(mu1, mu2)
X <- bivarchain(X,N)
b <- burn+1
x <- X[b:N,1];y <- X[b:N,2]
# Draw the scatter plot
plot(x,y,main="The scatter plot",cex=.5,xlab="X", ylab="Y", ylim=range(y))

# Fit a simple linear regression model
fit <- lm(y~x)
abline(fit,col="red")
summary(fit)

## ----setup821, fig.height=5, fig.width=10, echo=T, eval=T---------------------
plot(fit)
# The cor.test between x and the absolute residuals
abse<-abs(fit$residuals)
cor.test(x,abse,alternative="two.sided",method="spearman",conf.level=0.95)

## ----setup822, fig.height=10, fig.width=10, echo=T, eval=T--------------------
# Monitor the convergence of the chains
# The Gelman-Rubin method
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi);k <- nrow(psi)
  psi.means <- rowMeans(psi)
  B <- n*var(psi.means) # between variance est.
  psi.w <- apply(psi,1,"var") # within variances
  W <- mean(psi.w) # within est.
  v.hat <- W*(n-1)/n+(B/n) # upper variance est.
  r.hat <- v.hat/W # G-R statistic
  return(r.hat)
}

# Initialize constants and parameters
mu1 <- mu2 <- c(-5,-1,0,5) # Set the different initial value
x1 <- matrix(0,4,N)
y1 <- matrix(0,4,N)

for(i in 1:4){
  X1 <- matrix(0,N,2)# The chain, a bivariate sample
  X1[1,] <- c(mu1[i], mu2[i])
  X1 <- bivarchain(X1,N)
  x1[i,] <- X1[,1]
  y1[i,] <- X1[,2]
}

# Compute diagnostic statistics
psix <- t(apply(x1, 1, cumsum))
psiy <- t(apply(y1, 1, cumsum))
for(i in 1:nrow(psix)){
  psix[i,] <- psix[i,]/(1:ncol(psix))
}
for(i in 1:nrow(psiy)){
  psiy[i,] <- psiy[i,]/(1:ncol(psiy))
}
print(c(X=Gelman.Rubin(psix),Y=Gelman.Rubin(psiy)))

# Plot psi for the four chains
par(mfrow=c(2,2))
for(i in 1:4){
  plot(psix[i,b:N],type="l",xlab=i,ylab=bquote(psi[x]))
}
for(i in 1:4){
  plot(psiy[i,b:N],type="l",xlab=i,ylab=bquote(psi[y]))
}
par(mfrow=c(1,1)) #restore default

## ----setup823, fig.height=5, fig.width=10, echo=T, eval=T---------------------
# Plot the sequence of R-hat statistics
rhatx <- rhaty <-rep(0,N)
for(j in b:N){
  rhatx[j] <- Gelman.Rubin(psix[,1:j])
  rhaty[j] <- Gelman.Rubin(psiy[,1:j])
}
par(mfrow=c(1,2))
plot(rhatx[b:N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2, col="red")
plot(rhaty[b:N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2, col="red")
par(mfrow=c(1,1))

# Clean the memory of the variables
rm(list=ls())

## ----setup911, fig.height=10, fig.width=10, echo=T, eval=T--------------------
set.seed(123)

# The function to generate the random sample
RSample <- function(n,alpha,beta){
  X <- runif(n,10,20)
  gamma <- 1;aM <- 0.5;aY <- 1
  M <- aM+alpha*X+rnorm(n)
  Y <- aY+beta*M+gamma*X+rnorm(n)
  return(list(X,M,Y))
}

# The function of test statistics computation
Ttest <- function(X,M,Y){
  fit1 <- summary(lm(M~X))
  fit2 <- summary(lm(Y~X+M))
  a <- fit1$coefficients[2,1]
  sea <- fit1$coefficients[2,2]
  b <- fit2$coefficients[3,1]
  seb <- fit2$coefficients[3,2]
  return(a*b/((a*seb)^2+(b*sea)^2)^0.5)
}

# The function to implement the test hypothesis
Imptest <- function(N,n,X,M,Y,T0){
  T1 <- T2 <- T3 <- numeric(N)
  # Condition 1
  for(i in 1:N){
    n1 <- sample(1:n, size=n, replace=FALSE)
    n2 <- sample(1:n, size=n, replace=FALSE)
    X1 <- X[n1];M1 <- M[n2];Y1 <- Y[n2]
    T1[i] <- Ttest(X1,M1,Y1)
  }
  # Condition 2
  for(i in 1:N){
    n1 <- sample(1:n, size = n, replace = FALSE)
    n2 <- sample(1:n, size = n, replace = FALSE)
    X2 <- X[n1];M2 <- M[n1];Y2 <- Y[n2]
    T2[i] <- Ttest(X2,M2,Y2)
  }
  # Condition 3
  for(i in 1:N){
    n1 <- sample(1:n, size = n, replace = FALSE)
    n2 <- sample(1:n, size = n, replace = FALSE)
    M3 <- M[n1];X3 <- X[n2];Y3 <- Y[n2]
    T3[i] <- Ttest(X3,M3,Y3)
  }
  # The p-value of Condition1
  p1 <- mean(abs(c(T0,T1))>abs(T0))
  # The p-value of Condition2
  p2 <- mean(abs(c(T0,T2))>abs(T0))
  # The p-value of Condition3
  p3 <- mean(abs(c(T0,T3))>abs(T0))
  return(c(p1,p2,p3))
}

N <- 1000 # The number of simulation
n <- 100 # The number of random sample
T0 <- numeric(3)
p <- matrix(0,3,3)
# The real values of parameters
alpha <- c(0,0,1);beta <- c(0,1,0)

for(i in 1:3){
  result <- RSample(n,alpha[i],beta[i])
  X <- result[[1]]
  M <- result[[2]]
  Y <- result[[3]]
  # The original value of test statistics
  T0[i] <- Ttest(X,M,Y)
  p[i,] <- Imptest(N,n,X,M,Y,T0[i])
}

## ----setup912, fig.height=10, fig.width=10, echo=T, eval=T--------------------
# Result reporting
colnames(p) <- c("Condition 1","Condition 2","Condition 3")
rownames(p) <- c("alpha=0,beta=0","alpha=0,beta=1","alpha=1,beta=0")
p

# Clean the memory of the variables
rm(list=ls())

## ----setup92, fig.height=5, fig.width=10, echo=T, eval=T----------------------
set.seed(22058)

# The function to solve the alpha
solve <- function(N,b1,b2,b3,f0){
  x1 <- rpois(N,1)
  x2 <- rexp(N,1)
  x3 <- rbinom(N,1,0.5)
  g <- function(alpha){
    tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
  }
  solution <- uniroot(g,c(-50,0))
  return(round(unlist(solution),5)[1])
}

N <- 1e6;b1 <- 0;b2 <- 1;b3 <- -1
f0 <- c(0.1,0.01,0.001,0.0001)
alpha <- numeric(length(f0))
for(i in 1:length(f0)){
  alpha[i] <- solve(N,b1,b2,b3,f0[i])
}
result <- rbind(f0,alpha)
rownames(result) <- c("f0","alpha")
result

par(mfrow=c(1,2))
# Draw the scatter plot of f0 and alpha
plot(f0,alpha,main="The scatter plot of f0 and alpha")
# Draw the scatter plot of log(f0) and alpha
plot(log(f0),alpha,main="The scatter plot of log(f0) and alpha")
par(mfrow=c(1,1))

# Clean the memory of the variables
rm(list=ls())

## ----setup101, fig.height=10, fig.width=10, echo=T, eval=T--------------------
set.seed(22058)
# The data of sample
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- c(12,9,28,14,17,1,24,11,25,3)
n <- length(u)
# MLE
L <- function(lambda){
  tmp <- (v*exp(-lambda*v)-u*exp(-lambda*u))/(exp(-lambda*u)-exp(-lambda*v))
  sum(tmp)
}
lambdaMLE <- uniroot(L,c(0,10))

# EM algorithm
# Set the initial value
lambda0 <- 0
lambdaEM <- 1
while(lambdaEM!=lambda0){
  lambda0 <- lambdaEM
  lambdaEM <- n/(n/lambda0+sum((u*exp(-lambda0*u)-v*exp(-lambda0*v))/(exp(-lambda0*u)-exp(-lambda0*v))))
}

# Result reporting
print(c("MLE"=unlist(lambdaMLE)[1],
        "EM algorithm"=lambdaEM))

# Clear the memory of the variables
rm(list=ls())

## ----setup102, fig.height=10, fig.width=10, echo=T, eval=T--------------------
a <- list(2,2,0,5,8)
is.atomic(unlist(a))
is.atomic(as.vector(a))
is.list(as.vector(a))
# Clear the memory of the variables
rm(list=ls())

## ----setup103, fig.height=10, fig.width=10, echo=T, eval=T--------------------
1 == "1"
-1 < FALSE
"one" < 2

## ----setup104, fig.height=10, fig.width=10, echo=T, eval=T--------------------
a <- c(1,2,3,4,5)
dim(a)
# Clear the memory of the variables
rm(list=ls())

## ----setup105, fig.height=10, fig.width=10, echo=T, eval=T--------------------
a <- matrix(1,5,7)
is.matrix(a)
is.array(a)
class(a)
# Clear the memory of the variables
rm(list=ls())

## ----setup106, fig.height=10, fig.width=10, echo=T, eval=T--------------------
dfm <- data.frame(x=1:3, y=I(matrix(1:9,nrow=3)))
dfm
attributes(dfm)

## ----setup107, fig.height=10, fig.width=10, echo=T, eval=T--------------------
as.matrix(dfm)
# Clear the memory of the variables
rm(list=ls())

## ----setup108, fig.height=10, fig.width=10, echo=T, eval=T--------------------
a <- data.frame()
nrow(a)
ncol(a)
# Clear the memory of the variables
rm(list=ls())

## ----setup201, fig.height=10, fig.width=10, echo=T, eval=T--------------------
set.seed(22058)
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

# Apply to every column of a data frame
# Using the data frame mtcars
attach(mtcars)
result1 <- lapply(mtcars, scale01)
head(as.data.frame(result1))
detach(mtcars)

# Apply to every numeric column in a data frame
# Using the data frame iris
attach(iris)
dfmclass <- sapply(iris,is.numeric)
dfmclass <- as.numeric(which(dfmclass==TRUE))
result2 <- lapply(iris[dfmclass],scale01)
head(as.data.frame(result2))
detach(iris)

# Clear the memory of the variables
rm(list=ls())

## ----setup202, fig.height=10, fig.width=10, echo=T, eval=T--------------------
set.seed(22058)

# a) In a numeric data frame
# Using the data frame mtcars
attach(mtcars)
vapply(mtcars,sd,numeric(1))
detach(mtcars)

# b) In a mixed data frame
# Using the data frame iris
attach(iris)
dfmclass <- vapply(iris,is.numeric,logical(1))
dfmclass <- as.numeric(which(dfmclass==TRUE))
vapply(iris[dfmclass],sd,numeric(1))
detach(iris)

# Clear the memory of the variables
rm(list=ls())

## ----setup203, fig.height=5, fig.width=10, echo=T, eval=T---------------------
# Rcpp function
library(Rcpp)
cppFunction('NumericMatrix Cbivarchain(double a,double b,int N){
  double mu1,mu2,sigma1,sigma2,rho;
  double x1,x2,m1,s1,m2,s2;
  int i;
  NumericMatrix X(N,2);
  mu1 = mu2 = 0;sigma1 = sigma2 = 1;rho = 0.9;
  X(0,0) = a;X(0,1) = b;
  s1 = sqrt(1-pow(rho,2))*sigma1;
  s2 = sqrt(1-pow(rho,2))*sigma2;
  for(i=1;i<N;i++){
    x2 = X(i-1,1);
    m1 = mu1+rho*(x2-mu2)*sigma1/sigma2;
    X(i,0) = R::rnorm(m1,s1);
    x1 = X(i,0);
    m2 = mu2+rho*(x1-mu1)*sigma2/sigma1;
    X(i,1) = R::rnorm(m2,s2);
  }
  return X;
}')

# R function
Rbivarchain <- function(a,b,N){
  mu1 <- mu2 <- 0;sigma1 <- sigma2 <- 1;rho <- 0.9
  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2
  X <- matrix(0,N,2)
  X[1,] <- c(a,b)
  for (i in 2:N){
    x2 <- X[i-1,2]
    m1 <- mu1+rho*(x2-mu2)*sigma1/sigma2
    X[i,1] <- rnorm(1,m1,s1)
    x1 <- X[i,1]
    m2 <- mu2+rho*(x1-mu1)*sigma2/sigma1
    X[i,2] <- rnorm(1,m2,s2)
  }
  return(X)
}

# Generate random numbers
N <- 5000 # Length of chain
burn <- 1000 # Burn-in length
start <- burn+1
a <- b <- 0 # The initial value
XR <- Rbivarchain(a,b,N)
xR <- XR[start:N,1];yR <- XR[start:N,2]
XC <- Cbivarchain(a,b,N)
xC <- XC[start:N,1];yC <- XC[start:N,2]

# Draw Q-Q plots
par(mfrow=c(1,2))
qqplot(xR,xC,main="Q-Q plot of X")
qqline(xR,col="red")
qqplot(yR,yC,main="Q-Q plot of Y")
qqline(yR,col="red")
par(mfrow=c(1,1)) #restore default

## ----setup2031, fig.height=5, fig.width=10, echo=T, eval=T--------------------
# Compare the computation time
library(microbenchmark)
microbenchmark(Rbivarchain=Rbivarchain(a,b,N),
               Cbivarchain=Cbivarchain(a,b,N))
# Clear the memory of the variables
rm(list=ls())

