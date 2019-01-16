#' Generate random draws using MLHS algorithm
#'
#' Generate random draws using the Modified Latin Hypercube Sampling algorithm.
#'
#' Internal use only.
#' Algorithm described in
#' Hess, S., Train, K., and Polak, J. (2006) Transportation Research 40B, 147 - 163.
#' @param  N The number of draws to generate in each dimension
#' @param  d The number of dimensions to generate draws in
#' @param  i The number of individuals to generate draws for
#' @return A (N*i) x d matrix with random draws
apollo_mlhs=function(N,d,i){
  shuffle=function(inv){
    out=inv[rank(stats::runif(length(inv)))];
    out}

  temp=seq(0,N-1)/N;
  out=matrix(0,N*i,d);
  j=1;
  k=1;
  while(j<i+1){
    k=1;
    while(k<d+1){
      out[(1+N*(j-1)):(N*j),k]=shuffle(temp+stats::runif(1)/N);
      k=k+1
    }
    j=j+1
  }
  out
}
