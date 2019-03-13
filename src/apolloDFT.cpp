// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]


#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace stats;
#define Rcpp__stats__random_rnorm_h

// [[Rcpp::export]]
arma::mat RCPPinv(arma::mat M) {
  return  arma::inv(M);
}

// [[Rcpp::export]]
arma::mat RCPPpower(arma::mat S, const double& Sn) {
  arma::vec  Sva;
  arma::mat  Sve;
  arma::eig_sym(Sva, Sve, S);
  arma::mat SB = diagmat(pow(Sva,Sn));
  return Sve*SB*Sve.t();
}

// [[Rcpp::export]]
arma::mat RCPPeta(arma::mat S, const double& Sn, const double& o,arma::vec mu,arma::vec P0) {
  arma::mat St = RCPPpower(S,Sn);
  arma::mat OS = eye<mat>(o,o);
  arma::mat S1 = RCPPinv(OS-S)*((OS-St)*mu);
  S1 = S1+ St*P0;
  return S1;
}

// [[Rcpp::export]]
arma::mat RCPPZ(arma::mat S) {
  arma::mat Z = kron(S,S);
  return Z;
}

// [[Rcpp::export]]
arma::mat RCPPomega(arma::mat Z, const double& Zn, const double& o,const double& o2,arma::vec phiB) {
  arma::mat Zt = RCPPpower(Z,Zn);
  arma::mat OS = eye<mat>(o2,o2);
  arma::mat Z1 = RCPPinv(OS-Z)*((OS-Zt)*phiB);
  Z1.reshape(o,o);
  return Z1;
}

// [[Rcpp::export]]
arma::mat RCPPphiB(arma::mat phi, const double& o) {
  phi.reshape(o,1);
  return phi;
}

// [[Rcpp::export]]
arma::mat RCPPphi(arma::mat C,arma::mat M, arma::mat psi, const double& s) {
  arma::mat phi = C*M*psi*M.t()*C;
  phi.diag() += s;
  return phi;
}

// [[Rcpp::export]]
List DFTprob(arma::vec attr, arma::vec P0, arma::vec wts, double choice, double erv, double ts, double phi1, double phi2){
  double Alts = P0.n_elem;
  double AltsW = attr.n_elem;
  double W = wts.n_elem;
  double check1 = W*Alts;
  // throw error if matrix cannot be built from attributes
  if (check1!=AltsW) {
    throw std::invalid_argument ("The number of alternative attribute values passed in does not equal the number of alternatives multiplied by the number of attributes");
  }
  // rescale weights
  double wtsum = sum(wts);
  wts = wts/wtsum;
  arma::mat M = attr;
  M.reshape(Alts,W);
  //get distances:
  arma::mat N = eye<mat>(Alts,Alts);
  for (int i = 1; i <= Alts; ++i) {
    for (int j = 1; j <= Alts; ++j) {
    N(i-1,j-1) = sum(square(M.row(i-1) - M.row(j-1)));
    }
  }
  arma::mat w1 = wts;
  w1.reshape(W,1);
  arma::mat C = eye<mat>(Alts,Alts);
  C = C*(Alts)/(Alts-1);
  C = C- 1/(Alts-1);
  arma::mat mu = C*(M*w1);
  arma::mat psi = eye<mat>(W,W);
  for (int i = 1; i <= W; ++i) {
    psi(i-1,i-1) = wts(i-1);
  }
  w1.reshape(1,W);
  psi = psi-wts*w1;
  arma::mat phi = RCPPphi(C,M,psi,erv);
  arma::mat eta(Alts,1);
  arma::mat omega(Alts,Alts);
  arma::mat S(Alts,Alts);
  if(phi2==0){
    eta = ts*mu + P0;
    omega = ts*phi;
  } else {
    arma::mat OS = eye<mat>(Alts,Alts);

    S = OS - phi2* exp(-phi1*N);
    // check for size for off-diagonals.
    // too small results in invertible matrix
    double Scheck = 0;
    for (int i = 1; i <= Alts; ++i) {
      Scheck = Scheck + sum(abs(S.row(i-1)));
    }
    double Scheck2 = (1-phi2) * Alts;
    // check required amount to add,
    // and whether it is mathematically correct.
    Scheck2 = Scheck2 + 0.0000000001;
    double Scheck3 = 0;
    for (int i = 1; i <= Alts; ++i) {
      Scheck3 = Scheck3 + sum(abs(N.row(i-1)));
    }
    if(Scheck<Scheck2){
      eta = ts*mu + P0;
      omega = ts*phi;
    } else {

      if(Scheck3<0.000001){
      // basically have M=0 -> diff calculation?
        eta = ts*mu + P0;
        omega = ts*phi;
      } else {

        arma::mat phiB = RCPPphiB(phi,Alts*Alts);
        arma::mat Z = RCPPZ(S);
        eta = RCPPeta(S,ts,Alts,mu,P0);
        omega = RCPPomega(Z,ts,Alts,Alts*Alts,phiB);
      }
    }


  }
  arma::mat L = ones<mat>(Alts-1,Alts);
  arma::mat L1 = eye<mat>(Alts-1,Alts-1);
  L1 = -L1;
  int k = 1;
  for (int i = 1; i <= Alts; ++i) {
    if (choice!=i) {
      L.col(i-1)=L1.col(k-1);
      k=k+1;
    }
  }


  arma::mat gam = L*eta;
  arma::mat lambda = L*omega*L.t();
  //arma::vec R = {gam(0),gam(1),lambda(0,0),lambda(0,1),lambda(1,1)};

  double eta_check = sum(abs(eta.row(0)));
  if (eta_check<0.0000000001) {
    lambda.zeros();
  // R then know to assign probability equal to 1/Alt
  }

  return Rcpp::List::create(gam, lambda);
}


// [[Rcpp::export]]
NumericVector cpp_pmnorm(NumericVector x, NumericVector y, Function f) {
  NumericVector res = f(x, y);
  return res;
}

