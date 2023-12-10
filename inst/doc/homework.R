## ----echo=FALSE---------------------------------------------------------------
name=c("a","b","c","d","e","f")
v1=matrix(1,ncol=1,nrow=6)
v2=c(1,3,5,7,9,11)
v3=seq(0.1, 3, 0.5)
v4=1:6
datajf <- data.frame(name, v1, v2, v3,v4)
datajf  #运行R代码后展示的结果
library(knitr)
knitr::kable(datajf,col.names=c('name', 'v1', 'v2', 'v3','v4'),format='latex',align=rep('c', 5)) #利用knitr语句展示的结果, 使用latex格式, align=rep('c', 5)是让表格内容居中

## ----echo=FALSE---------------------------------------------------------------
y=c(4,2.9,5.6,6,3.9,4,5,5,4,5,6)
x=c(0,1,6,7.8,2,3.5,6,7.6,4.5,7,8)
model=lm(y~x)
summary(model)


## ----pressure, echo=FALSE-----------------------------------------------------
plot(y~x)
abline(model,col="red") 

## -----------------------------------------------------------------------------
n <- 1000
u1 <- runif(n,0,0.5)
u2 <- runif(n,0.5,1)
x1=log(2*u1)
x2=-log(2*u2-1)
x=c(x1,x2)
hist(x, prob = TRUE, main = expression("f(x)=1/2e^{-|x|}"),ylim = c(0, 0.55))
y <- seq(-100, 100, 0.01)
lines(y, 0.5*exp(-abs(y)))

## -----------------------------------------------------------------------------
rejacc_beta=function(n,a,b){
j=0; k=0; y=NULL;
while (k < n) {
j=j+1;
u=runif(1) #random variate u from U(0,1)
x=runif(1) #random variate x from g(.)
if (x^2*(1-x)>u) {
k=k+1
y[k]=x}
}
return(y)
}

## ----echo=FALSE---------------------------------------------------------------
n=1000; a=3; b=2
y=rejacc_beta(n,a,b)
hist(y, prob = TRUE, main = expression("Beta distribution"),xlim=c(0,1),ylim = c(0, 2),xlab='x',ylab='Probability density function')
set.seed(123)
x<-seq(0,1,length.out=1000)
lines(x,dbeta(x,3,2),col='red',lwd=2)

## -----------------------------------------------------------------------------
n=100000  # a large simulated random sample
genefe=function(n){
u1=runif(n,-1,1)
u2=runif(n,-1,1)
u3=runif(n,-1,1)
x=NULL;
for (i in 1:n){
if ((abs(u3[i])>=abs(u2[i])) && (abs(u3[i])>=abs(u1[i]))) {
  x=cbind(x,u2[i])
} else {
  x=cbind(x,u3[i])}
}
return(x)
}

## ----echo=FALSE---------------------------------------------------------------
x=genefe(n)
hist(x, prob = TRUE, main = expression("Rescaled Epanechnikov kernel"),xlim=c(-1,1),ylim = c(0, 1),xlab='x',ylab='Probability density function')  ## the histogram density
set.seed(123)
x<-seq(-1,1,length.out=n)
lines(x,3/4*(1-x^2),col='red',lwd=2)

## -----------------------------------------------------------------------------
library(base)
mysample=function(x, size, prob=NULL){
if (is.null(prob)==TRUE) prob=rep(1/length(x),length(x)) 
             ##Set the default to equal probability sampling
cp = cumsum(prob)  
U = runif(size)
r <- x[findInterval(U,cp)+1] # The inverse transform method
return(r)
}

## -----------------------------------------------------------------------------
x=c(1,2,3,4)
size=10
prob=c(0.1,0.2,0.3,0.4)
mysample(x,size, prob)   #Sampling results
table(mysample(x,size, prob)) # collect sampling information

## -----------------------------------------------------------------------------
x=c('a','b','c')
size=40
prob=c(0.8,0.1,0.1)
mysample(x,size,prob)  #Sampling results
table(mysample(x,size,prob)) # collect sampling information

## -----------------------------------------------------------------------------
K=100
pihat1=0; pihat2=0; pihat3=0;
for (i in 1:K){
set.seed(123+K)

l1=1; l2=1.5; l3=2; d=2; 

rho1=l1/d
rho2=l2/d 
rho3=l3/d ## rho_min

m = 1e6
X = runif(m,0,d/2)
Y = runif(m,0,pi/2)
pihat1 = pihat1 + 2*rho1/mean(l1/2*sin(Y)>X)
pihat2 = pihat2 + 2*rho2/mean(l2/2*sin(Y)>X)
pihat3 = pihat3 + 2*rho3/mean(l3/2*sin(Y)>X)
}
pihat1/K; pihat2/K; pihat3/K

## -----------------------------------------------------------------------------

MCjf = function(R, antithetic) 
{
  u = runif(R/2) ##Generate random numbers Xi from the uniform distrubution U(0,1), i = 1,...,m/2
  if (antithetic==TRUE) {v=1-u} ##calucate Yi = 1 - Xi, i = 1,...,m/2
    else {v=runif(R/2)} ##Generate random numbers Xi from the uniform distrubution U(0,1), i = m/2+1,...,m
  u = c(u, v)
  g = exp(u) 
  return (mean(g))
}
m = 10000  ##the number of simulations
MC1=matrix(0,ncol=1,nrow=m)  ## the simple Monte Carlo method
MC2=matrix(0,ncol=1,nrow=m)  ## the antithetic variate approach
for (i in 1:m) {
  set.seed(123+i)
  MC1[i] = MCjf(R = m, antithetic = FALSE)
  MC2[i] = MCjf(R = m, antithetic = TRUE)
}
var(MC1); var(MC2); (var(MC1)-var(MC2))/var(MC1)


## ----fig.align='center'-------------------------------------------------------
##构建函数f1与f2
f1 = function(x) {
  2*sqrt(1/(2*pi))*exp(-(x-1)^2/2)
}
n=2
f2 = function(x) {
  2*factorial(1/2)/(sqrt(pi*2))*(1+(x-1)^2/2)^(-3/2)
##自由度为n时 2*factorial((n+1)/2-1)/(factorial(n/2-1)*sqrt(pi*n))*(1+(x-1)^2/n)^(-(n+1)/2)
}
intf1=integrate(f1, 1, Inf) #积分结果为1，说明可以作为概率密度函数
intf2=integrate(f2, 1, Inf) #积分结果为1，说明可以作为概率密度函数

x = seq(1, 30, 0.01)
y = x^2 * exp(-x^2/2)/sqrt(2 * pi) 
plot(x, y, type = "l", ylim = c(0, 1))
lines(x, 2*dnorm(x,mean = 1), lty = 2, col='blue')
f2=2*factorial(1/2)/(sqrt(pi*2))*(1+(x-1)^2/2)^(-3/2)
lines(x, f2, lty = 3, col='red')
legend("topright", legend = c("g(x)", "f1", "f2"), col=c("black","blue","red"), lty = 1:3)


## ----fig.align='center'-------------------------------------------------------
x = seq(1, 30, 0.01)
plot(x, y/(2*dnorm(x,mean=1)), type = "l", lty = 3, col='blue',ylab='')
f2=2*factorial(1/2)/(sqrt(pi*2))*(1+(x-1)^2/2)^(-3/2)
lines(x, y/f2, lty = 2, col='red')
legend("topright", inset = 0.02, legend = c("g/f1", "g/f2"), col=c("blue","red"), lty = 2:3)

## -----------------------------------------------------------------------------
set.seed(123)
m = 1e5
est = sd = numeric(2)
g = function(x) {
x^2 * exp(-x^2/2)/sqrt(2 * pi)
}

#f1
x = rnorm(m,mean=1) 
result1 = g(x)/(2*dnorm(x,mean=1))
est[1] = mean(result1)
sd[1] = sd(result1)

#f2
x = runif(m,1,100) 
f2=2*factorial(1/2)/(sqrt(pi*2))*(1+(x-1)^2/2)^(-3/2)
result2 = g(x)/f2
est[2] = mean(result2)
sd[2] = sd(result2)
rbind(est,sd)
integrate(function(x) x^2*dnorm(x), 1, Inf)$value ##g(x)的真实值

## -----------------------------------------------------------------------------
set.seed(123)
 n = 10000  
 k = 5  #分成5层
 m = n/k 
theta = numeric(k)
var = numeric(k)
sd = numeric(k) 
 
 g = function(x) exp(-x)/(1 + x^2)*(x>0)*(x<1)
 f = function(x) k*exp(-x)/(1-exp(-1))
 for (j in 1:k) {
 u = runif(m, (j - 1)/k, j/k)
 x = -log(1 - (1 - exp(-1)) * u) # 由逆变换方法得到对应的值
 result = g(x)/f(x)
 theta[j] = mean(result)
  var[j] = var(result)
  sd[j] = sd(result)
 }
sum(theta)
mean(var)

## -----------------------------------------------------------------------------
set.seed(123)
 n = 10000  
 k = 1  #不分层
 m = n/k 
theta1 = numeric(k)
var1 = numeric(k)
sd1 = numeric(k) 
 
 g = function(x) exp(-x)/(1 + x^2)*(x>0)*(x<1)
 f = function(x) k*exp(-x)/(1-exp(-1))
 for (j in 1:k) {
 u = runif(m, (j - 1)/k, j/k)
 x = -log(1 - (1 - exp(-1)) * u) # 由逆变换方法得到对应的值
 result = g(x)/f(x)
 theta1[j] = mean(result)
  var1[j] = var(result)
  sd1[j] = sd(result)
 }
sum(theta1)
mean(var1)

## -----------------------------------------------------------------------------
m=10000
CIl=numeric(m)
CIr=numeric(m)
 n <- 20  ## sample size
 t0 <- qt(0.975, df = n - 1)  #t_{n-1}(\alpha/2), alpha=0.5
 for (i in 1:m){
x <-rchisq(n, df = 2)
CIl[i] <- mean(x)-t0*sd(x)/sqrt(n)  # Lower bound of the interval
CIr[i] <- mean(x)+t0*sd(x)/sqrt(n)  # Upper bound of the interval
}
 sum(CIl < 2 & CIr > 2)  #The mean of the chi-square distribution χ2(2) is 2
 mean(CIl < 2 & CIr > 2)

## -----------------------------------------------------------------------------
m <-1e4
n<-1e3  
alpha <- 0.05 
pval1 <- matrix(0,1,m)
pval2 <- matrix(0,1,m)
pval3 <- matrix(0,1,m)
for(i in 1:m) {
  set.seed(i+123)
  xchi <- rchisq(n,df=1)
  xuni <- runif(n,0,2)
  xexp <- rexp(n,rate=1)
  t01=abs(mean(xchi)-1)/(sd(xchi)/ sqrt(n))
  t02=abs(mean(xuni)-1)/(sd(xuni)/ sqrt(n))
  t03=abs(mean(xexp)-1)/(sd(xexp)/ sqrt(n))
  pval1[i] <- 2*pt(1-t01,n-1)
  pval2[i] <- 2*pt(1-t02,n-1)
  pval3[i] <- 2*pt(1-t03,n-1)
}
c(mean(pval1<=alpha),mean(pval2<=alpha),mean(pval3<=alpha))

## -----------------------------------------------------------------------------
m=1000 #有1000个原假设
m0=950 #前95%个原假设为真
B=1000 #重复次数
alpha=0.1 #显著性水平
V1=0;  V2=0; FDR_1=0; FDR_2=0; TPR_1=0; TPR_2=0

for (k in 1:B){
p0=runif(m0,0,1)
p1=rbeta(m-m0,shape1=0.1,shape2=1)
p=c(p0,p1)
##Bonferroni校正
pb=p.adjust(p, method = "bonferroni")  #所有p值得放在一起校正
pb0=pb[1:m0]
pb1=pb[(m0+1):m]

##B-H校正, 也就是FDR校正
ph=p.adjust(p, method = "BH")  #所有p值得放在一起校正
ph0=ph[1:m0]
ph1=ph[(m0+1):m]


if (sum(pb0<=alpha)>0) {V1=V1+1}
if (sum(ph0<=alpha)>0) {V2=V2+1}
FDR_1=FDR_1+sum(pb0<=alpha)/(sum(pb0<=alpha)+sum(pb1<=alpha))
FDR_2=FDR_2+sum(ph0<=alpha)/(sum(ph0<=alpha)+sum(ph1<=alpha))
TPR_1=TPR_1+sum(pb1<=alpha)/(m-m0)
TPR_2=TPR_2+sum(ph1<=alpha)/(m-m0)

}
FWER_1=V1/B
FWER_2=V2/B
FDR_1=FDR_1/B
FDR_2=FDR_2/B
TPR_1=TPR_1/B
TPR_2=TPR_2/B

ans=matrix(0,ncol=3,nrow=2)
ans[1,1]=FWER_1; ans[2,1]=FWER_2;
ans[1,2]=FDR_1; ans[2,2]=FDR_2;
ans[1,3]=TPR_1; ans[2,3]=TPR_2;
colnames(ans) <- c("FWER","FDR","TPR")
rownames(ans) <- c("Bonf","BH")
ans


## -----------------------------------------------------------------------------
bootstrapjf=function(n)
{
lambda=2
B=1000  # number of bootstrap replicates
m=1000  # repeated times
bias <- se.boot <-numeric(m)

for (i in 1:m){
set.seed(123+i)
x = rexp(n, rate=lambda)
lambdahatstar <- numeric(B)
lambdahat <- 1/mean(x); #要估计的lambda是1/\bar{x}
for(b in 1:B){
xstar <- sample(x,replace=TRUE)
lambdahatstar[b] <- 1/mean(xstar) #要估计的lambda是1/\bar{x}
}
bias[i]=mean(lambdahatstar)-lambdahat
se.boot[i]=sd(lambdahatstar)
}
return(list(bias=mean(bias),se.boot=mean(se.boot),bias_theotical=lambda/(n-1),se.boot_theotical=lambda*n/((n-1)*sqrt(n-2))))
}
ans1=bootstrapjf(n=5)
ans2=bootstrapjf(n=10)
ans3=bootstrapjf(n=20)

ans=matrix(0,ncol=5,nrow=3)
ans[1,1]=5; ans[2,1]=10; ans[3,1]=20;
ans[1,2]=ans1$bias; ans[1,4]=ans1$se.boot; ans[1,3]=ans1$bias_theotical; ans[1,5]=ans1$se.boot_theotical;
ans[2,2]=ans2$bias; ans[2,4]=ans2$se.boot; ans[2,3]=ans2$bias_theotical; ans[2,5]=ans2$se.boot_theotical;
ans[3,2]=ans3$bias; ans[3,4]=ans3$se.boot; ans[3,3]=ans3$bias_theotical; ans[3,5]=ans3$se.boot_theotical;
colnames(ans) <- c("n","Bias_bootstrap","Bias_theoretical","SE_bootstrap","SE_theoretical")
ans

## -----------------------------------------------------------------------------
library(bootstrap) #for the law data
x=law$LSAT
y=law$GPA
thetahat=cor(x,y)
n=NROW(x)

B=1000  # number of bootstrap replicates
m=1000  # number of  second order bootstrap replicates
t <-numeric(B)
thetahatstar <-numeric(B)
for(b in 1:B){
bh <- sample(1:n,replace=TRUE)
xstar <- x[bh]
ystar <- y[bh]
thetahatstar[b] <- cor(xstar,ystar) 
thetahatstarstar <-numeric(m)
for(i in 1:m){
bhstar <- sample(1:n,replace=TRUE)
xstarstar <- xstar[bhstar]
ystarstar <- ystar[bhstar]
thetahatstarstar[i] <- cor(xstarstar,ystarstar) 
}
t[b]=(thetahatstar[b]-thetahat)/sd(thetahatstarstar) #计算t值
}
alpha=0.05
t_l=quantile(t, probs=1-alpha/2)
t_r=quantile(t, probs=alpha/2)
CL=thetahat-t_l*sd(thetahatstar)
CR=thetahat-t_r*sd(thetahatstar)

## -----------------------------------------------------------------------------
library(boot)
set.seed(123)
data = aircondit$hours
lambda.inverse = function(x,i) mean(x[i])
bootobject = boot(data = data,statistic = lambda.inverse,R=10000)
CI = boot.ci(bootobject, type = c("norm", "perc", "basic", "bca"), conf=0.95)
data.frame("Method"=c("Normal", "Basic", "Percentile", "BCa"),"CI_l"=c(CI$norm[2],CI$basic[4],CI$percent[4],CI$bca[4]),"CI_r"=c(CI$norm[3],CI$basic[5],CI$percent[5],CI$bca[5]))

## ----fig.align='center'-------------------------------------------------------
hist(bootobject$t, prob = TRUE, main = "") 
points(bootobject$t0, 0, cex = 2, pch = 5)

## -----------------------------------------------------------------------------
library(bootstrap)
 data=scor
 n=nrow(data)
 theta_jack = numeric(n)
 lambda = eigen(cov(data))$values
 theta_hat = max(lambda)/sum(lambda)  #\hat{\theta}
 for (i in 1:n) {
 Y = data[-i, ]
 lambda = eigen(cov(Y))$values
 theta_jack[i] = max(lambda)/sum(lambda)  #\hat{\theta}_{(i)}
 }
 bias_jack = (n - 1) * (mean(theta_jack) - theta_hat)
 se_jack = sqrt((n-1)*mean((theta_jack-theta_hat)^2))
 c(theta_hat, bias_jack, se_jack)

## -----------------------------------------------------------------------------
 library(DAAG)
 data=ironslag
 n = NROW(data)
 N = choose(n, 2)
 chemical=data[,1]
 magnetic=data[,2]
 e=matrix(0,nrow=1,ncol=4)
 for (i in 1:(n-1))
   for (j in (i+1):n)
  {
  x = chemical[c(-i,-j)]
  y = magnetic[c(-i,-j)]

 ## Linear: Y = β0 + β1X + epsilon

 model1 = lm(y ~ x)
 magnetic_hat1 = model1$coef[1] + model1$coef[2] * chemical[c(i,j)]
 e[1] = e[1]+sum((magnetic[c(i,j)] - magnetic_hat1)^2)

 ## Quadratic: Y = β0 + β1X + β2X^2 + epsilon

 model2 = lm(y ~ x + I(x^2))
 magnetic_hat2 = model2$coef[1] + model2$coef[2] * chemical[c(i,j)] +model2$coef[3] * chemical[c(i,j)]^2
 e[2] = e[2]+sum((magnetic[c(i,j)] - magnetic_hat2)^2)

 ## Exponential: log(Y) = β0 + β1X +  epsilon

 model3 = lm(log(y) ~ x)
 magnetic_hat3 = exp(model3$coef[1] + model3$coef[2] * chemical[c(i,j)])
 e[3] = e[3]+sum((magnetic[c(i,j)] - magnetic_hat3)^2)

 ## Log-log: log(Y) = β0 + β1log(X) + epsilon
 
 model4 = lm(log(y) ~ log(x))
 magnetic_hat4 = exp(model4$coef[1] + model4$coef[2] * log(chemical[c(i,j)]))
 e[4] = e[4]+sum((magnetic[c(i,j)] - magnetic_hat4)^2)

}
e=e/N
e

## -----------------------------------------------------------------------------
x1=chickwts[1:10,1]  #feed=horsebean
x2=chickwts[11:22,1]  #feed=linseed
x3=chickwts[23:36,1]  #feed=soybean

teststat=function(x,y){
n = length(x)
m = length(y)
R = 200
z=c(x,y)
stat0=0
stat =numeric(R)
 test=function(x,y)
{
 s1=0; s2=0;
 for (i in 1:n) {
 s1=s1+(sum(as.integer(x<=x[i]))/n-sum(as.integer(y<=x[i]))/m)^2
 }
 for (j in 1:m) {
 s2=s2+(sum(as.integer(x<=y[j]))/n-sum(as.integer(y<=y[j]))/m)^2
  }
 return ((n*m)/(n+m)^2*(s1+s2))
}
 stat0 =test(x,y) 

 for (k in 1:R){
 kk = sample(1:(n+m))
 zr = z[kk]
 xr = zr[1:n]
 yr = zr[(n+1):(n+m)]
 stat[k] =test(xr,yr)
}
 stat1 = c(stat, stat0)
 statistic = stat0
 p_value = mean(stat1 >= stat0)
 return (list=c(statistic,p_value))
}
ans1=teststat(x1,x2)
ans2=teststat(x1,x3)
ans3=teststat(x2,x3)
data.frame("Objects"=c("horsebean-linseed","horsebean-soybean","linseed-soybean"),"statistic"=c(ans1[1],ans2[1],ans3[1]),"p_value"=c(ans1[2],ans2[2],ans3[2]))

## -----------------------------------------------------------------------------

 outliers = function(x, y) {
 X = x - mean(x)
 Y = y - mean(y)
 outlierx = sum(X > max(Y)) + sum(X < min(Y))
 outliery = sum(Y > max(X)) + sum(Y < min(X))
 return(max(c(outlierx, outliery))) }

 test = function(x, y, R = 300) {
 z = c(x, y)
 nx = length(x)
 ny = length(y)
 stats=numeric(R)
 for (i in 1:R)
 {
 k = sample(1:(nx+ny))
 kx = k[1:nx]
 ky = k[(nx + 1):(nx+ny)]
 stats[i]=outliers(z[kx], z[ky]) 
 }
 stat0 = outliers(x, y)
 statp = c(stats, stat0)
 return(list(estimate = stat0, p = mean(statp >= stat0)))
 }
 
 
 set.seed(123)
 x = rnorm(500, 2, 1)
 y = rnorm(1000, 3, 1)
 test(x, y)

 
 set.seed(123)
 x = rnorm(500, 2, 3)
 y = rnorm(1000, 3, 4)
 test(x, y)

## -----------------------------------------------------------------------------
set.seed(123)
N = 1e6; b1 = 0; b2 = 1; b3=-1; 
x1 = rpois(N, 1 ) 
x2 = rexp(N, 1)
x3=rbinom(N, 1, 0.5)
g1 = function(a){mean(1/(1+exp(a+b1*x1+b2*x2+b3*x3))) - 0.1}
g2 = function(a){mean(1/(1+exp(a+b1*x1+b2*x2+b3*x3))) - 0.01}
g3 = function(a){mean(1/(1+exp(a+b1*x1+b2*x2+b3*x3))) - 0.001}
g4 = function(a){mean(1/(1+exp(a+b1*x1+b2*x2+b3*x3))) - 0.0001}
solution1 = uniroot(g1,c(0,10))
solution2 = uniroot(g2,c(0,10))
solution3 = uniroot(g3,c(0,10))
solution4 = uniroot(g4,c(0,10))
c(solution1$root,solution2$root,solution3$root,solution4$root)

## -----------------------------------------------------------------------------
 rwjf = function(N, x0, sigma) {
 x=numeric(N)
 x[1]=x0
 u=runif(N)
 k = 0
 for (i in 2:N) { 
 y = rnorm(1, x[i-1], sigma)
 if (u[i] <= exp(abs(x[i-1]) - abs(y))) {x[i]=y}
 else { x[i] = x[i-1];  k=k+1; }
 }
 return(list(x=x, k=k))
 }
 N=5000
 rw1=rwjf(N, rnorm(1,0, 0.2), 0.2)
 rw2=rwjf(N, rnorm(1,0, 0.5), 0.5)
 rw3=rwjf(N, rnorm(1,0, 1.7), 1.7)
 rw4=rwjf(N, rnorm(1,0, 3), 3)
 c(rw1$k, rw2$k, rw3$k, rw4$k)
 burn=500
 y1=rw1$x[(burn+1):N]
 y2=rw2$x[(burn+1):N]
 y3=rw3$x[(burn+1):N]
 y4=rw4$x[(burn+1):N]
 plot(y1, type = "l")
 plot(y2, type = "l")
 plot(y3, type = "l")
 plot(y4, type = "l")

## -----------------------------------------------------------------------------
N=5000
burn=1000  ##discard a suitable burn-in sample
data=matrix(0,N,2)
rho=0.9
mux=0
muy=0
sigma_x=1
sigma_y=1
s1=sqrt(1-rho^2)*sigma_x
s2=sqrt(1-rho^2)*sigma_y
data[1,]=c(mux,muy)
for(i in 2:N){
yy=data[i-1,2]
m1=mux+rho*sigma_x/sigma_y*(yy-muy)
data[i,1]=rnorm(1,m1,s1)
xx=data[i,1]
m2=muy+rho*sigma_y/sigma_x*(xx-mux)
data[i,2]=rnorm(1,m2,s2)
}
X=data[(burn+1):N,1]
Y=data[(burn+1):N,2]
model=lm(Y~X)  ##fit a simple linear regression model
model
plot(model)

## -----------------------------------------------------------------------------
 Gelman_Rubin=function(phi) {
 phi=as.matrix(phi)
 n=ncol(phi)
 B=n * var(rowMeans(phi)) #between variance est.
 W=mean(apply(phi, 1, "var"))  #within est.
 v_hat=W * (n - 1)/n + (B/n)
 r_hat=v_hat/W #G-R statistic
 return(r_hat)
 }

 f = function(x, sigma) {
 if (x < 0)  {return(0)}
 return((x/sigma^2) * exp(-x^2/(2 * sigma^2)))
}
 Rayleigh_MH = function(sigma, m, x0) {
 x = numeric(m);  x[1] = x0;  u = runif(m);
 for (i in 2:m) {
 xt = x[i-1]
 y = rchisq(1, df = xt)
 num = f(y, sigma) * dchisq(xt, df = y)
 den = f(xt, sigma) * dchisq(y, df = xt)
 alpha=num/den
 if (u[i] <= alpha)
 x[i] = y
 else x[i] = xt
 }
 return(x)
 }
 sigma=2
 x0 = c(1/sigma^3, 1/sigma, sigma^2, sigma^4)
 k = 4
 m = 2000
 X = matrix(0, nrow = k, ncol = m)
 for (i in 1:k) 
   X[i, ] = Rayleigh_MH(sigma, m, x0[i])
 phi = t(apply(X, 1, cumsum))
 for (i in 1:nrow(phi)) 
   phi[i, ] = phi[i, ]/(1:ncol(phi))
 rhat=Gelman_Rubin(phi)
 rhat

## -----------------------------------------------------------------------------
 library(coda)
 X1=as.mcmc(X[1, ])
 X2=as.mcmc(X[2, ])
 X3=as.mcmc(X[3, ])
 X4=as.mcmc(X[4, ])
 Y=mcmc.list(X1, X2, X3, X4)
 gelman.plot(Y, col = c(1, 1))

## -----------------------------------------------------------------------------
##########   MLE
###Interval data
u=c(11,8,27,13,16,0,23,10,24,2)
v=c(12,9,28,14,17,1,24,11,25,3)

##Initialize
eps=10e-5
MLElambda=0.5
minus=1
while (abs(minus)>eps)
{
  MLElambda1=MLElambda-(sum(v-u*exp(MLElambda)))/sum(-u*exp(MLElambda))
  minus=MLElambda1-MLElambda
  MLElambda=MLElambda1
}

MLElambda 

## -----------------------------------------------------------------------------
########## EM
##Interval data
u=c(11,8,27,13,16,0,23,10,24,2)
v=c(12,9,28,14,17,1,24,11,25,3)
n=10 #sample size

##Initialize
eps=10e-5
EMlambda=0.5
minus=1

while (abs(minus)>eps)
{
  EMlambda1=n*EMlambda/(n-EMlambda*sum((v-u*exp(EMlambda))/(exp(EMlambda)-1)))
  minus=EMlambda1-EMlambda
  EMlambda=EMlambda1
}
EMlambda

## -----------------------------------------------------------------------------
solvejf=function(A){
min.A=min(A);A=A-min.A
max.A=max(A);A=A/max(A)
m=nrow(A);n=ncol(A)
it=n^3
a=c(rep(0,m),1)
A1=-cbind(t(A),rep(-1,n))
b1=rep(0,n)
A3=t(as.matrix(c(rep(1,m),0)))
b3=1
sx=simplex(a=a,A1=A1,b1=b1,A3=A3,b3=b3,maxi=TRUE,n.iter=it)
a=c(rep(0,n),1)
A1=cbind(A,rep(-1,m)); A3=t(as.matrix(c(rep(1,n),0)))
b1=rep(0,m); b3=1
sy=simplex(a=a,A1=A1,b1=b1,A3=A3,b3=b3,maxi=FALSE,n.iter=it)
soln=list(A=A*max.A+min.A,x=sx$soln[1:m],y=sy$soln[1:n],v=sx$soln[m+1]*max.A+min.A)
soln
}
##pay off matrix
A=matrix(c(0,-2,-2,3,0,0,4,0,0,2,0,0,0,
-3,-3,4,0,0,2,0,0,3,0,0,0,-4,-4,-3,
0,-3,0,4,0,0,5,0,0,3,0,-4,0,-4,0,5,
0,0,3,0,0,4,0,-5,0,-5,-4,-4,0,0,0,
5,0,0,6,0,0,4,-5,-5,0,0,0,6,0,0,4,
0,0,5,-6,-6,0),9,9)
library(boot)
B=A+2
s=solvejf(B)
 s$v  # the value of the game
 round(cbind(s$x, s$y), 4) ## optimal strategies for each player

## -----------------------------------------------------------------------------
x = c(1,2)
y = c("a","b")
new = data.frame(col1 = x,col2 = y)
mat = as.matrix(new) ##apply as.matrix()
is.matrix(mat)

## -----------------------------------------------------------------------------
dim(data.frame(x=1)[0,,drop=FALSE])  #0 rows
dim(data.frame(x=1)[,0,drop=FALSE])  #0 columns

## -----------------------------------------------------------------------------
scale01 = function(x) {
rng = range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
data1 = data.frame(matrix(data = 1:9,nrow=3,ncol=3))
data.frame(lapply(data1, function(x) scale01(x)))
data2 = data.frame(a = c(1,2,3),b = c(TRUE,TRUE,FALSE), c = c("a","b","c"))
data.frame(lapply(data2, function(x) if (is.numeric(x)) scale01(x) else x))

## -----------------------------------------------------------------------------
sdjf = function(data){ 
  a = vapply(data,FUN = function(x) sd(x), numeric(1))
  return(data.frame(a))
}
data1 = data.frame(matrix(data = 1:9,nrow=3,ncol=3))
sdjf(data1)

## -----------------------------------------------------------------------------
sdmixjf = function(data){ 
  data = data[vapply(data,FUN = is.numeric,logical(1))]
  a = vapply(data,FUN = function(x) sd(x), numeric(1))
  return(data.frame(a))
}
data2 = data.frame(a = c(1,2,3),b = c(TRUE,TRUE,FALSE), c = c("a","b","c"))
sdmixjf(data2)

## -----------------------------------------------------------------------------
library("Rcpp")
 library(microbenchmark)
set.seed(123)
jfR = function(N,n,a,b){
 x=y=rep(0, N)
 x[1] = rbinom(1, prob = 0.2, size = n)
 y[1] = rbeta(1, x[1] + a, n - x[1] + b)
 for (i in 2:N) {
 x[i] = rbinom(1, prob = y[i-1], size = n)
 y[i] = rbeta(1, x[i]+a, n-x[i]+b)
 }
 return(list(x=x,y=y))
}
sourceCpp(code=' 
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix jfC(int N, double n, double a, double b) {
  NumericMatrix aa(N,2);
  aa(0,0) = R::rbinom(n,0.2);
  aa(0,1) = R::rbeta(aa(0,0)+a, n-aa(0,0)+b);
 for(int i = 1; i < N; i++) {
   aa(i,0) = R::rbinom(n,aa(i-1,1));
   aa(i,1) = R::rbeta(aa(i,0)+a, n-aa(i,0)+b);
 }
 return(aa);
}
')
 N = 1000;  a = 4;  b=3; n = 6
 ## compare time
 t = microbenchmark(ar=jfR(N,n,a,b), aC=jfC(N,n,a,b))
 summary(t)[,c(1,3,5,6)]

