// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RCPPinv
arma::mat RCPPinv(arma::mat M);
RcppExport SEXP _apollo_RCPPinv(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPinv(M));
    return rcpp_result_gen;
END_RCPP
}
// RCPPpower
arma::mat RCPPpower(arma::mat S, const double& Sn);
RcppExport SEXP _apollo_RCPPpower(SEXP SSEXP, SEXP SnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double& >::type Sn(SnSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPpower(S, Sn));
    return rcpp_result_gen;
END_RCPP
}
// RCPPeta
arma::mat RCPPeta(arma::mat S, const double& Sn, const double& o, arma::vec mu, arma::vec P0);
RcppExport SEXP _apollo_RCPPeta(SEXP SSEXP, SEXP SnSEXP, SEXP oSEXP, SEXP muSEXP, SEXP P0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double& >::type Sn(SnSEXP);
    Rcpp::traits::input_parameter< const double& >::type o(oSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type P0(P0SEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPeta(S, Sn, o, mu, P0));
    return rcpp_result_gen;
END_RCPP
}
// RCPPZ
arma::mat RCPPZ(arma::mat S);
RcppExport SEXP _apollo_RCPPZ(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPZ(S));
    return rcpp_result_gen;
END_RCPP
}
// RCPPomega
arma::mat RCPPomega(arma::mat Z, const double& Zn, const double& o, const double& o2, arma::vec phiB);
RcppExport SEXP _apollo_RCPPomega(SEXP ZSEXP, SEXP ZnSEXP, SEXP oSEXP, SEXP o2SEXP, SEXP phiBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const double& >::type Zn(ZnSEXP);
    Rcpp::traits::input_parameter< const double& >::type o(oSEXP);
    Rcpp::traits::input_parameter< const double& >::type o2(o2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phiB(phiBSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPomega(Z, Zn, o, o2, phiB));
    return rcpp_result_gen;
END_RCPP
}
// RCPPphiB
arma::mat RCPPphiB(arma::mat phi, const double& o);
RcppExport SEXP _apollo_RCPPphiB(SEXP phiSEXP, SEXP oSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const double& >::type o(oSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPphiB(phi, o));
    return rcpp_result_gen;
END_RCPP
}
// RCPPphi
arma::mat RCPPphi(arma::mat C, arma::mat M, arma::mat psi, const double& s);
RcppExport SEXP _apollo_RCPPphi(SEXP CSEXP, SEXP MSEXP, SEXP psiSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPphi(C, M, psi, s));
    return rcpp_result_gen;
END_RCPP
}
// DFTprob
List DFTprob(arma::vec attr, arma::vec P0, arma::vec wts, double choice, double erv, double ts, double phi1, double phi2);
RcppExport SEXP _apollo_DFTprob(SEXP attrSEXP, SEXP P0SEXP, SEXP wtsSEXP, SEXP choiceSEXP, SEXP ervSEXP, SEXP tsSEXP, SEXP phi1SEXP, SEXP phi2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type attr(attrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type P0(P0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wts(wtsSEXP);
    Rcpp::traits::input_parameter< double >::type choice(choiceSEXP);
    Rcpp::traits::input_parameter< double >::type erv(ervSEXP);
    Rcpp::traits::input_parameter< double >::type ts(tsSEXP);
    Rcpp::traits::input_parameter< double >::type phi1(phi1SEXP);
    Rcpp::traits::input_parameter< double >::type phi2(phi2SEXP);
    rcpp_result_gen = Rcpp::wrap(DFTprob(attr, P0, wts, choice, erv, ts, phi1, phi2));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pmnorm
NumericVector cpp_pmnorm(NumericVector x, NumericVector y, Function f);
RcppExport SEXP _apollo_cpp_pmnorm(SEXP xSEXP, SEXP ySEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pmnorm(x, y, f));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_apollo_RCPPinv", (DL_FUNC) &_apollo_RCPPinv, 1},
    {"_apollo_RCPPpower", (DL_FUNC) &_apollo_RCPPpower, 2},
    {"_apollo_RCPPeta", (DL_FUNC) &_apollo_RCPPeta, 5},
    {"_apollo_RCPPZ", (DL_FUNC) &_apollo_RCPPZ, 1},
    {"_apollo_RCPPomega", (DL_FUNC) &_apollo_RCPPomega, 5},
    {"_apollo_RCPPphiB", (DL_FUNC) &_apollo_RCPPphiB, 2},
    {"_apollo_RCPPphi", (DL_FUNC) &_apollo_RCPPphi, 4},
    {"_apollo_DFTprob", (DL_FUNC) &_apollo_DFTprob, 8},
    {"_apollo_cpp_pmnorm", (DL_FUNC) &_apollo_cpp_pmnorm, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_apollo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
