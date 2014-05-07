// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;
// using namespace arma;
// using namespace Eigen;

using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Upper;
using Eigen::VectorXd;
typedef Map<MatrixXd> MapMatd;
typedef Map<MatrixXi> MapMati;
typedef Map<VectorXd> MapVecd;


// myFastLm
List myFastLm(MapMatd X, MapMatd y);
RcppExport SEXP simmen_myFastLm(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< MapMatd >::type X(XSEXP );
        Rcpp::traits::input_parameter< MapMatd >::type y(ySEXP );
        List __result = myFastLm(X, y);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// ddirichlet
NumericVector ddirichlet(NumericVector xx, NumericVector alphaa);
RcppExport SEXP simmen_ddirichlet(SEXP xxSEXP, SEXP alphaaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alphaa(alphaaSEXP );
        NumericVector __result = ddirichlet(xx, alphaa);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rdirichlet
NumericVector rdirichlet(NumericVector alpha);
RcppExport SEXP simmen_rdirichlet(SEXP alphaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        NumericVector __result = rdirichlet(alpha);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// lik_random_weight
NumericVector lik_random_weight(List W_list, List z_list, NumericVector alpha);
RcppExport SEXP simmen_lik_random_weight(SEXP W_listSEXP, SEXP z_listSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type W_list(W_listSEXP );
        Rcpp::traits::input_parameter< List >::type z_list(z_listSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        NumericVector __result = lik_random_weight(W_list, z_list, alpha);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// updateWeight
List updateWeight(List W_list, List z_list, List x_list, NumericVector alpha, NumericVector beta, NumericVector sigma, NumericVector y, NumericVector tau, NumericVector uniform_rv);
RcppExport SEXP simmen_updateWeight(SEXP W_listSEXP, SEXP z_listSEXP, SEXP x_listSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP sigmaSEXP, SEXP ySEXP, SEXP tauSEXP, SEXP uniform_rvSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type W_list(W_listSEXP );
        Rcpp::traits::input_parameter< List >::type z_list(z_listSEXP );
        Rcpp::traits::input_parameter< List >::type x_list(x_listSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type uniform_rv(uniform_rvSEXP );
        List __result = updateWeight(W_list, z_list, x_list, alpha, beta, sigma, y, tau, uniform_rv);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// test
List test(arma::vec x);
RcppExport SEXP simmen_test(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP );
        List __result = test(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// update_e_internal_single
List update_e_internal_single(List location_index_all, List location_index1, List location_index2, arma::vec e, arma::vec new_e, arma::vec xb1, arma::vec xb2, arma::vec delta, arma::vec prob_old1, arma::vec prob_old2, arma::vec ystar1, arma::vec ystar2, arma::vec e_i, arma::vec e_j, arma::vec uniform_rv, arma::vec prob_diff, double sigma2);
RcppExport SEXP simmen_update_e_internal_single(SEXP location_index_allSEXP, SEXP location_index1SEXP, SEXP location_index2SEXP, SEXP eSEXP, SEXP new_eSEXP, SEXP xb1SEXP, SEXP xb2SEXP, SEXP deltaSEXP, SEXP prob_old1SEXP, SEXP prob_old2SEXP, SEXP ystar1SEXP, SEXP ystar2SEXP, SEXP e_iSEXP, SEXP e_jSEXP, SEXP uniform_rvSEXP, SEXP prob_diffSEXP, SEXP sigma2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type location_index_all(location_index_allSEXP );
        Rcpp::traits::input_parameter< List >::type location_index1(location_index1SEXP );
        Rcpp::traits::input_parameter< List >::type location_index2(location_index2SEXP );
        Rcpp::traits::input_parameter< arma::vec >::type e(eSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type new_e(new_eSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type xb1(xb1SEXP );
        Rcpp::traits::input_parameter< arma::vec >::type xb2(xb2SEXP );
        Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type prob_old1(prob_old1SEXP );
        Rcpp::traits::input_parameter< arma::vec >::type prob_old2(prob_old2SEXP );
        Rcpp::traits::input_parameter< arma::vec >::type ystar1(ystar1SEXP );
        Rcpp::traits::input_parameter< arma::vec >::type ystar2(ystar2SEXP );
        Rcpp::traits::input_parameter< arma::vec >::type e_i(e_iSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type e_j(e_jSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type uniform_rv(uniform_rvSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type prob_diff(prob_diffSEXP );
        Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP );
        List __result = update_e_internal_single(location_index_all, location_index1, location_index2, e, new_e, xb1, xb2, delta, prob_old1, prob_old2, ystar1, ystar2, e_i, e_j, uniform_rv, prob_diff, sigma2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// computeNetworkSummary_cxx
List computeNetworkSummary_cxx(arma::mat seq_m, arma::mat D);
RcppExport SEXP simmen_computeNetworkSummary_cxx(SEXP seq_mSEXP, SEXP DSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type seq_m(seq_mSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP );
        List __result = computeNetworkSummary_cxx(seq_m, D);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// genLogDetGamma
NumericVector genLogDetGamma(List W_list, NumericVector lambda);
RcppExport SEXP simmen_genLogDetGamma(SEXP W_listSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type W_list(W_listSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP );
        NumericVector __result = genLogDetGamma(W_list, lambda);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Mahalanobis
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov);
RcppExport SEXP simmen_Mahalanobis(SEXP xSEXP, SEXP centerSEXP, SEXP covSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP );
        Rcpp::traits::input_parameter< arma::rowvec >::type center(centerSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type cov(covSEXP );
        arma::vec __result = Mahalanobis(x, center, cov);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// dmvnorm_arma
arma::vec dmvnorm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool log = false);
RcppExport SEXP simmen_dmvnorm_arma(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP logSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP );
        Rcpp::traits::input_parameter< arma::rowvec >::type mean(meanSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< bool >::type log(logSEXP );
        arma::vec __result = dmvnorm_arma(x, mean, sigma, log);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
