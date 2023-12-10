#' @import stats
#' @import quadprog
#' @importFrom mvtnorm rmvnorm
#' @import glmnet
#' @import microbenchmark
#' @import bootstrap
#' @import boot
#' @import DAAG
#' @import coda 
#' @importFrom Rcpp evalCpp cppFunction
#' @importFrom graphics grid par title
#' @useDynLib SA23204163
NULL

#' @title TLMMA
#' @name TLMMA
#' @description TLMMA is a method to improve the estimation of target model by using multi-source domain samples under the background of transfer learning.
#' @param h Measure the degree of similarity between the source model and the target model parameters in set A
#' @param p Parameter dimension
#' @param s The The number of non-zero elements in the parameter
#' @param n0 The sample size of target model
#' @param nk The sample size of source models
#' @param M The The number of source models
#' @param sigbeta One of the initialization values of the parameter
#' @param sigx The covariance matrix of X
#' @param key 1 or 2 represents different data settings
#' @return MSE pMSE weight betahat
#' @examples
#' \dontrun{
#' p=20;  h=2; s=11;  n0=100;  nk=50;  M=10;  sigbeta=0.4; sigx=diag(1,p); key=1
#' TLMMA(h,p,s,n0,nk,M,sigbeta,sigx,key)
#' }
#' @export
TLMMA <-  function(h,p,s,n0,nk,M,sigbeta,sigx,key)
{
  ## Initialize the data
  cs=100;# 100 independent repetitions
  betahat0=0;
  sizeA0  =  0  # The number of elements of set A
  n_vec = c(n0,rep(nk,M)) # Sample size of target sample and source samples 
  MSE = matrix(0,nrow = M+1,ncol = 1)# MSE
  pMSE = matrix(0,nrow = M+1,ncol = 1)# PMSE
  weight = matrix(0,nrow = M+1,ncol = M+1)# The weight assigned to each candidate model
  
  
  for (sizeA0 in 0:M)  
  {
    A0 = 1: sizeA0
    cppFunction(" NumericMatrix gdata(int h, int p,int s, int M, double sigbeta, int key, int sizeA0, int nk)
 {
   // Initialize the data
   NumericMatrix B(p,M+1);
   double c0=p/5;
   double d0=trunc(M*2/3);
   double e0;
  
   int k,i,j;
   for(i=0; i<p; i++) {
     if (i<s) {B(i,0)=sigbeta;}
     else {B(i,0)=0;}
     for(j=1; j<=M; j++) 
     {B(i,j)=B(i,0);}  
   }
   
   if (key==1) {
     for(i=0; i<trunc(c0); i++) {
       for(k=1; k<=M; k++) 
       { 
       if (k<=sizeA0) {e0=(double)1.0*h/c0; }
          else {e0=(double)2.0*s/c0; }
	    B(i,k)= B(i,k)+R::rnorm(0,e0);
		}
     }
   } 
   
   e0=(double)1.0/nk;
   if (key==2) {
     for(i=0; i<p; i++){
       for(k=1; k<=M; k++) 
       { 
       
	   B(i,k) =B(i,k) + R::rnorm(0,e0);
	   if (k>d0) {B(i,k)=B(i,k)+sigbeta/2;} 
	   }
        // source models that perform particularly poorly
     }
   }
   return(B);
 }")
    
    
    for(t in 1:cs)
    {
    beta=gdata(h,p,s,M,sigbeta,key,sizeA0,nk)
    beta0=beta[,1]
    B=beta[,2:(M+1)]
    
    d = matrix(0,nrow = M+1, ncol = 1)
    D = matrix(0,nrow = n0, ncol = M+1)
    D0 = matrix(0,nrow = M+1, ncol = M+1)
    b =  NULL
    x0 = rmvnorm(n0, rep(0, p), sigx)   # x0 is the target sample
    y0 = x0 %*% beta0 + rnorm (n0, 0, 1) # y0 is the target sample
    xc = rmvnorm(n0, rep(0, p), sigx)   # xc is the test sample 
    
    ##Step1 least squares estimation
    fit = lm(y0~x0+0)         
    bb = coefficients(fit)      
    bb[is.na(bb)] = 0 
    b = cbind(b,bb)
    
    A = matrix(rep(y0,M+1),ncol = M+1,nrow = n0, byrow = F)
    for (k in 1:M) {
      xk = rmvnorm(n_vec[k+1], rep(0, p), sigx)  #x1 to xM are auxiliary samples
      yk = xk %*% B[,k] + rnorm (n_vec[k+1], 0, 1) 
      
      fit = lm(yk~xk+0)        
      bb = coefficients(fit)       
      bb[is.na(bb)] = 0 
      b = cbind(b,bb)
    }
    
    ##Step2 Calculate D,D0,d,a3,a4 to prepare for quadratic programming
    
    ehat = y0 - x0 %*% b[,1]          
    sighat = (t(ehat) %*% ehat)/(n0-p)       ## Estimates of the variance matrix
    
    D = A-x0%*%b
    D0 = 2*t(D)%*%D
    d = c(-2*sighat*p,rep(0,M))
    a3 = t(rbind(matrix(1,nrow = 1,ncol = M+1),diag(M+1),-diag(M+1)))       
    a4 = rbind(1,matrix(0,nrow = M+1,ncol = 1),matrix(-1,nrow = M+1,ncol = 1)) 
    
    
    ##Step3 Quadratic programming solution
    QP = solve.QP(Dmat = D0,dvec = d,Amat = a3, bvec = a4,meq = 1,factorized = FALSE)  
    w = QP$solution                        
    w = as.matrix(w)                       
    w = w*(w>0)                           
    w = w/sum(w)                          
    
    betahat = b %*% w                 ## Estimate the weighted parameters
    
    ## Step4 Calculate MSE and PMSE
    
    MSE[sizeA0+1] = MSE[sizeA0+1]+sum((betahat-beta0)^2) 
    pMSE[sizeA0+1] = pMSE[sizeA0+1]+mean((xc%*%(beta0-betahat))^2)  
    weight[,(sizeA0+1)] = weight[,(sizeA0+1)]+w
    betahat0=betahat0+betahat
    }
    
    MSE[sizeA0+1]=MSE[sizeA0+1]/cs
    pMSE[sizeA0+1]=pMSE[sizeA0+1]/cs
    weight[,(sizeA0+1)]=weight[,(sizeA0+1)]/cs
    betahat0=betahat0/cs
  }
  return(list(MSE = MSE,pMSE = pMSE,weight = weight, betahat = betahat0))
}




#' @title plotwto1
#' @name plotwto1
#' @description A Graph which used to verify the theory of the method TLMMA.
#' @param h Measure the degree of similarity between the source model and the target model parameters in set A
#' @param p Parameter dimension
#' @param s The The number of non-zero elements in the parameter
#' @param M The The number of source models
#' @param sigbeta One of the initialization values of the parameter
#' @param sigx The covariance matrix of X
#' @return A graph about the trend of the sum of weights of potential auxiliary models as the sample size increases
#' @examples
#' \dontrun{
#' h=2;p=20;s=14;M=10;sigbeta=0.4;sigx=diag(1,p);
#' plotwto1(h,p,s,M,sigbeta,sigx)
#' }
#' @export
plotwto1 <- function(h,p,s,M,sigbeta,sigx){
  n0=c(100,500,1000,5000,10000,20000,50000) 
  nk=n0
  wright=matrix(0,ncol=1,nrow=NROW(n0))
  for (i in 1:NROW(n0)){
    result1=TLMMA(h,p,s,n0[i],nk[i],M,sigbeta,sigx,2)
    ww=result1$weight[,(M+1)] # when sizeA0=M
    wright[i]=sum(ww[1:(floor(M*2/3)+1)]) # The sum of weights of the potential auxiliary models
  }
  plot(n0,wright,"b",lty=3,lwd=2,pch=1,cex=1.5,main="The sum of weights of the potential auxiliary models varies with sample size",xlab="Sample size",ylab="The sum of weights of the potential auxiliary models")
  grid(col='orange',lty=1)
}


#' @title plotwto2
#' @name plotwto2
#' @description Graphs are used to verify the performance of the method TLMMA.
#' @param h Measure the degree of similarity between the source model and the target model parameters in set A
#' @param p Parameter dimension
#' @param s The The number of non-zero elements in the parameter
#' @param n0 The sample size of target model
#' @param nk The sample size of source models
#' @param M The The number of source models
#' @param sigbeta One of the initialization values of the parameter
#' @param sigx The covariance matrix of X
#' @param key 1=Calculation method result 2=Verify relevant theories
#' @return A graph about MSE and PMSE of the TLMMA method
#' @examples
#' \dontrun{
#' h=6;p=20;s=11;n0=100;nk=100;M=20;sigbeta=0.3;sigx=diag(1,p);key=1;
#' plotwto2(h,p,s,n0,nk,M,sigbeta,sigx,key)
#' }
#' @export
plotwto2 <- function(h,p,s,n0,nk,M,sigbeta,sigx,key){
  result2=TLMMA(h,p,s,n0,nk,M,sigbeta,sigx,key)
  MSE=result2$MSE
  pMSE=result2$pMSE
  par(mfrow=c(1,2))
  plot(1:(M+1),MSE,"b",col='red',lty=3,lwd=2,pch=1,cex=1.5,main ='MSE varies with the The number of source models in set A',xlab="The number of source models in set A",ylab="MSE")
  grid(col='orange',lty=1)
  plot(1:(M+1),pMSE,"b",col='blue',lty=3,lwd=2,pch=1,cex=1.5,main ='pMSE varies with the The number of source models in set A',xlab="The number of source models in set A",ylab="pMSE")
  grid(col='orange',lty=1)
}


