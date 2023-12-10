#include <Rcpp.h>
using namespace Rcpp;
//' @title gdata
//' @name gdata
//' @description A function that generates parameters of source models and target model 
//' @param h Measure the degree of similarity between the source model and the target model parameters in set A
//' @param p Parameter dimension
//' @param s The number of non-zero elements in the parameter
//' @param M The number of source models
//' @param sigbeta One of the initialization values of the parameter
//' @param key 1 or 2 represents different data settings
//' @param sizeA0 Number of source models in set A
//' @param nk The sample size of source models
//' @return B Generated parameters of target model and source models 
//' @examples
//' \dontrun{
//' h=2; p=20; s=11; M=10; sigbeta=0.4; key=1;  sizeA0=4; nk=50
//' B <- gdata(h,p,s,M,sigbeta,key,sizeA0,nk)
//' }
//' @export
// [[Rcpp::export]]
 NumericMatrix gdata(int h, int p,int s, int M, double sigbeta, int key, int sizeA0, int nk)
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
 }
 
