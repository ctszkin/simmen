// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>

using namespace Rcpp;

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

MatrixXd AtA(const MapMatd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
  .rankUpdate(A.adjoint());
}

// [[Rcpp::export]]
List myFastLm(MapMatd X, MapMatd y){
  const int n(X.rows()), p(X.cols());
  const LLT<MatrixXd> llt(AtA(X));
  const VectorXd betahat(llt.solve(X.adjoint() * y));
  const VectorXd fitted(X * betahat);
  const VectorXd resid(y - fitted);
  const int df(n - p);
  const double s(resid.norm() / std::sqrt(double(df)));
  const MatrixXd cov(llt.matrixL().solve(MatrixXd::Identity(p, p)));
  const MatrixXd cov2(MatrixXd(p,p).setZero().selfadjointView<Lower>().rankUpdate(cov.adjoint()) );

  return List::create(
    Named("coefficients") = betahat,
    Named("fitted.values") = fitted,
    Named("residuals") = resid,
    Named("s") = s,
    Named("df.residual") = df,
    Named("rank") = p,
    Named("cov") = cov2*s*s); 
}


// [[Rcpp::export]]
NumericVector ddirichlet(NumericVector xx , NumericVector alphaa){
  
  if (is_true(any(alphaa<=0))){
    return wrap(-1e20);
  }

  if (is_true(any( (xx <=0) ))){
    return wrap(-1e20);
  }

  double alpha_sum = std::accumulate(alphaa.begin(), alphaa.end(), 0.0); 
  
  // this will return alpha_sum = 0 0 0
  //return (alpha_sum);


  NumericVector log_gamma_alpha = lgamma(alphaa);
  NumericVector log_alpha_sum = lgamma(wrap(alpha_sum));


  double sum_log_gamma_alpha = std::accumulate(log_gamma_alpha.begin(), log_gamma_alpha.end(), 0.0); 
  
  double logD = sum_log_gamma_alpha - as<double>(log_alpha_sum);


  NumericVector a = ( alphaa - 1.0 ) * log(xx);
  return wrap( std::accumulate(a.begin(), a.end(), 0.0) - logD ); 
}

// [[Rcpp::export]]
NumericVector rdirichlet(NumericVector alpha){
  int n = alpha.size();
  NumericVector out(n);
  for (int i=0 ; i <n; i++){
    out[i] = as<double>(rgamma(1, alpha(i) ));
  }
  out = out / std::accumulate(out.begin(), out.end(), 0.0);
  return out;
}


// [[Rcpp::export]]
NumericVector lik_random_weight(List W_list, List z_list, NumericVector alpha){
  int n = W_list.size();
  NumericVector lik(n);

  arma::colvec alpha_arma = as<arma::colvec>(alpha);

  for (int i=0 ; i < n ; i++){
    SEXP temp_w = W_list[i];
    SEXP temp_z = z_list[i];
    NumericVector w(temp_w);
    arma::mat z = Rcpp::as<arma::mat>(temp_z);
    arma::colvec za = exp(  z * alpha_arma) ;

    lik(i) = ddirichlet( w, wrap(za) )[0] ;
  }
  if (all(is_finite(lik))){
    return wrap(std::accumulate(lik.begin(), lik.end(), 0.0));
  } else{
    return wrap(-1e20);
  }

}


// [[Rcpp::export]]
List updateWeight(List W_list, List z_list, List x_list, NumericVector alpha, NumericVector beta, NumericVector sigma, NumericVector y, NumericVector tau, NumericVector uniform_rv){

  int n = W_list.size();
  int update=0;

  arma::colvec alpha_arma = as<arma::colvec>(alpha);
  arma::colvec beta_arma = as<arma::colvec>(beta);
  Rcpp::List out(n);

  for (int i=0 ; i < n ; i++){
    SEXP temp_w = W_list[i];
    SEXP temp_z = z_list[i];
    SEXP temp_x = x_list[i];

    arma::mat z = Rcpp::as<arma::mat>(temp_z);
    arma::colvec za = exp( z * alpha_arma) ;
    arma::rowvec w_old = as<arma::rowvec>(temp_w);
    arma::mat x = Rcpp::as<arma::mat>(temp_x);
    arma::colvec xb = w_old * x * beta_arma ;

    NumericVector res(1);
    res(0) = y(i) - xb(0);

    NumericVector lik_old = ddirichlet(wrap(w_old) , wrap(za) ) + dnorm(res , 0, sigma) ;

    NumericVector temp_alpha = wrap(w_old);

    temp_alpha = temp_alpha* tau(0);
    arma::rowvec w_new = as<arma::rowvec>( rdirichlet( temp_alpha ) );

    arma::colvec xb_new = w_new * x * beta_arma ;

    NumericVector res_new(1);
    res_new(0) = y(i) - xb_new(0);

    NumericVector lik_new = ddirichlet(wrap(w_new) , wrap(za) ) +  dnorm(res_new, 0, sigma)  ;

    NumericVector w_new_tau = wrap(w_new);
    NumericVector w_old_tau = wrap(w_old);
    w_new_tau = w_new_tau*tau(0);
    w_old_tau = w_old_tau*tau(0);
   
    NumericVector Q = ddirichlet( wrap(w_old) , w_new_tau) - ddirichlet( wrap(w_new) , w_old_tau);

    NumericVector a = exp(lik_new - lik_old + Q);

   //  return wrap( all(!is_na( a ) ) );

    if ( ::R_finite( as<double>(a) ) ){
      if ( a(0) > uniform_rv(i) ){
          update++;
          out[i] = wrap(w_new );
      } else {
        out[i] = wrap(w_old);
      }      
    } else {
      out[i] = wrap(w_old);
    }
  
  }
 return Rcpp::List::create(Rcpp::Named("W")=out, Rcpp::Named("update")=wrap(update) );
}

// [[Rcpp::export]]
List test(arma::vec x){
  return Rcpp::List::create(Rcpp::Named("e")=x-mean(x));
}

// [[Rcpp::export]]
List update_e_internal_single(
  List location_index_all,
  List location_index1,
  List location_index2,
  arma::vec e,
  arma::vec new_e,
  arma::vec xb1,
  arma::vec xb2,
  arma::vec delta,
  arma::vec prob_old1,
  arma::vec prob_old2,
  arma::vec ystar1,
  arma::vec ystar2,
  arma::vec e_i,
  arma::vec e_j,
  arma::vec uniform_rv,
  arma::vec prob_diff,
  double sigma2
  )
{
  Rcpp::List location_index_all_list(location_index_all);
  int n = location_index_all_list.size();

  //int update_count=0;
  NumericVector update_count(n);
  for (int i=0; i <n ; i++){
    SEXP l1 = location_index_all[i]; 
    SEXP l2 = location_index1[i]; 
    SEXP l3 = location_index2[i]; 

    arma::uvec location_index_all_vec = Rcpp::as<arma::uvec>(l1)-1;
    arma::uvec location_index1_vec = Rcpp::as<arma::uvec>(l2)-1;
    arma::uvec location_index2_vec= Rcpp::as<arma::uvec>(l3)-1;

    arma::vec e_i_new = e_i;
    arma::vec e_j_new = e_j;

    e_i_new(location_index1_vec).fill(new_e(i)) ;
    e_j_new(location_index2_vec).fill(new_e(i)) ; 

    // arma::vec diff_e_new = square((e_i_new(location_index_all_vec)- e_j_new(location_index_all_vec)));
    arma::vec diff_e_new = abs((e_i_new(location_index_all_vec)- e_j_new(location_index_all_vec)));

    //diff_e_new = diff_e_new - mean(diff_e_new);

    arma::vec m1 = ystar1(location_index_all_vec) - xb1(location_index_all_vec) - diff_e_new * delta;
    arma::vec m2 = ystar2(location_index_all_vec) - xb2(location_index_all_vec) - diff_e_new * delta;

    // compute normal density 
    arma::vec prob_new1 = log( 1/sqrt(2*M_PI*sigma2) * exp(-0.5* square(m1)/sigma2 ) ) ; 
    arma::vec prob_new2 = log( 1/sqrt(2*M_PI*sigma2) * exp(-0.5* square(m2)/sigma2 ) ) ; 


    arma::vec prob_old11_subset = prob_old1(location_index_all_vec);
    arma::vec prob_old22_subset = prob_old2(location_index_all_vec);

    double a1 = std::accumulate(prob_new1.begin(),prob_new1.end(), 0.0);
    double a2 = std::accumulate(prob_new2.begin(),prob_new2.end(), 0.0);
    double a3 = std::accumulate(prob_old11_subset.begin(),prob_old11_subset.end(), 0.0);
    double a4 = std::accumulate(prob_old22_subset.begin(),prob_old22_subset.end(), 0.0);


    double alpha = exp( a1+a2-a3-a4 + prob_diff(i)); 


    if (alpha > arma::as_scalar(uniform_rv(i)) ){
      e_i = e_i_new;
      e_j = e_j_new;
      prob_old1(location_index_all_vec) = prob_new1 ;
      prob_old2(location_index_all_vec) = prob_new2;
      e(i) = new_e(i);
      update_count(i) = 1;
    }
  }

  return Rcpp::List::create(Rcpp::Named("e")=e, Rcpp::Named("update_rate")=update_count );
}


// [[Rcpp::export]]
List computeNetworkSummary_cxx(arma::mat seq_m, arma::mat D){
    int nn = D.n_rows;
    arma::mat D00 = arma::zeros(nn,nn);
    arma::mat degreee = arma::zeros(nn,nn);
    arma::mat common_frds_11 = arma::zeros(nn,nn);

    for ( int i=0 ; i<nn*(nn-1)/2; i++ ){
      int index1 = seq_m(i,0) -1 ;
      int index2 = seq_m(i,1) -1;
      if (D(index1,index2)==1) {
        D00(index1,index2) = 1;
        D00(index2,index1) = 1;
      }

      degreee(index1,index2) = sum(D00.col(index1)) ; 
      degreee(index2,index1) = sum(D00.col(index2))  ;
      common_frds_11(index1,index2) = arma::as_scalar(D00.col(index1).t() * D00.col(index2)) ; 
    }
    arma::mat out1 = arma::zeros(nn*(nn-1)/2, 3);
    arma::mat out2 = arma::zeros(nn*(nn-1)/2, 3);


    int k = 0;
    for ( int j=0 ; j < nn ; j++){
      for ( int i=j+1  ; i< nn ; i++){
        out1(k,0) = arma::as_scalar( degreee(i,j) );
        out1(k,1) = arma::as_scalar( degreee(i,j)*degreee(i,j) );
        out1(k,2) = arma::as_scalar( common_frds_11(i,j) );

        out2(k,0) = arma::as_scalar( degreee(j,i) );
        out2(k,1) = arma::as_scalar( degreee(j,i)*degreee(j,i) );
        out2(k,2) = arma::as_scalar( common_frds_11(i,j) );

        k++;
      }
    }

    return Rcpp::List::create(Rcpp::Named("self")=out1, Rcpp::Named("friends")=out2);
}

// [[Rcpp::export]]
NumericVector genLogDetGamma(List W_list, NumericVector lambda){
  int m = W_list.size();

  SEXP ll = W_list[0];
  NumericMatrix temp(ll);

  int n = temp.ncol();
  arma::mat out;
  out.zeros(n,n);

  for (int i=0; i < m ; i++ ){
    SEXP ll = W_list[i];
    arma::mat temp = as<arma::mat>(ll) ;

    out -=  temp * lambda[i];
  }

  arma::colvec v ;
  v.ones(n,1);

  out.diag() = v;
  double value,sign;

  arma::log_det(value,sign,out);

  return wrap( value );
}

// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov){
    int n = x.n_rows;
    arma::mat x_cen;
    x_cen.copy_size(x);
    for (int i=0; i < n; i++) {
        x_cen.row(i) = x.row(i) - center;
    }
    return sum((x_cen * cov.i()) % x_cen, 1);    
}

// [[Rcpp::export]]
arma::vec dmvnorm_arma(arma::mat x,  arma::rowvec mean,  arma::mat sigma, bool log = false) { 
    arma::vec distval = Mahalanobis(x,  mean, sigma);
    double logdet = sum(arma::log(arma::eig_sym(sigma)));
    double log2pi = std::log(2.0 * M_PI);
    arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
    
    if (log){ 
        return(logretval);
    }else { 
        return(exp(logretval));
    }
}
