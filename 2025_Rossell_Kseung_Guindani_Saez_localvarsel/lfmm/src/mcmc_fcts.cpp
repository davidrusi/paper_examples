#include <RcppArmadillo.h>
#include <progress.hpp>
#include <R_ext/Utils.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include <rgen.h> 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(rgen)]]


//' Dirichlet Distribution
//'
//' Dirichlet distribution random number generator.
//'
//' @param deltas Vector of scale parameters.
//' @return One sample from the Dirichlet distribution.
//' @export
// [[Rcpp::export]]
arma::vec rdirichlet(arma::vec deltas){
  arma::vec out = rgen::rdirichlet(deltas);
  
  return (out);
}


//' Count Uniques
//' 
//' Analogous of the table function in R: computes the counts of the unique 
//' elements from 1 to K.
//'
//' @param x Input vector.
//' @param K Number of levels of the vector x.
//' @return Vector of counts.
//' @export
// [[Rcpp::export]]
arma::uvec table_int(arma::uvec x, unsigned int K){
  
  unsigned int n = x.n_elem;
  arma::uvec count(K, arma::fill::zeros); // here we consider that x has values in {1, ..., max(x)}
  
  for (unsigned int i = 0; i < n; i++){
    count[x[i] - 1] += 1;
  }
  return (count);
}


//' Stable Normalization
//' 
//' Normalize a vector on the log scale with the log-sum-exp trick.
//' 
//' @param x Vector (on the log-scale) to be normalized.
//' @return Normalized vector on the exponentiated scale.
//' @export
// [[Rcpp::export]]
arma::vec normalize_log(arma::vec x){
  unsigned int n = x.n_elem;
  
  double c = x.max();
  x -= c;
  
  arma::vec out(n);
  out = arma::exp(x)/arma::accu(arma::exp(x));
  
  return (out);
}


//' Stable Sum
//' 
//' Sum the rows of a matrix on the log scale.
//' 
//' @param Q Matrix (on the log-scale) whose rows have to be summed over.
//' @return Vector containing the row sums.
//' @export
// [[Rcpp::export]]
arma::vec sum_rows_log(arma::mat Q){
  
  unsigned int S = Q.n_cols;
  arma::vec out(S, arma::fill::zeros);
  double m_star; 
  
  for (unsigned hh = 0; hh < S; hh++){
    m_star = max(Q.row(hh));
    if (m_star != - arma::datum::inf){
      out(hh) = m_star + log(arma::accu(exp(Q.row(hh) - m_star)));
    }
    else{
      out(hh) = - arma::datum::inf;
    }
  }
  
  return (out);
}


//' Count transitions
//' 
//' Count the number of transitions from a matrix of labels evolving over time.
//' 
//' @param z Matrix whose rows represent different covariate levels and whose 
//'          columns are the times.
//' @return Matrix counting the number of transitions.
// [[Rcpp::export]]
arma::umat count_assign(arma::umat z){
  
  unsigned int d_1 = z.n_rows;
  unsigned int K = z.n_cols;
  
  arma::umat count(d_1, d_1, arma::fill::zeros); 
  for (unsigned int i = 0; i < d_1; i++){
    for (unsigned int j = 1; j < K; j++){
      count(z(i,j-1) - 1,z(i,j) - 1) += 1;
    }
  }
  return (count);
}


//' Backpropagation
//' 
//' Backpropagate the messages of the HMM.
//' 
//' @param y Vector of observations.
//' @param beta_star Matrix of possible spline coefficients.
//' @param sigma2_y Residual variance.
//' @param mess_old Backward messages for the HMM.
//' @param Q_old Transition probability matrix.
//' @param sigma2_1_beta Smoothness inducing parameter.
//' @return Sampled values for the latent variables.
// [[Rcpp::export]]
arma::mat backprop_mess_lin(arma::vec y, arma::mat beta_star, double sigma2_y, 
                            arma::vec mess_old, arma::mat Q_old, 
                            double sigma2_1_beta){
  
  unsigned int S = beta_star.n_rows;
  arma::mat out(S, S, arma::fill::zeros);
  
  for (unsigned int i = 0; i < S; i++){
    for (unsigned int j = 0; j < S; j++){
      out(i,j) = log(mess_old(j)) -
        0.5/sigma2_y * arma::accu(arma::square(y - beta_star(j,1))) -
        0.5 / sigma2_1_beta * std::pow((beta_star(i,0) - beta_star(j,1)), 2.0) +
        log(Q_old(i,j));
    }
  }
  
  return (out);
}


//' Random Effects Update
//' 
//' Full conditional distribution to update the random effects spline coefficients.
//' 
//' @param res vector of residuals (observations minus main effects)
//' @param B B-spline basis matrix corresponding to the vector of residuals.
//' @param P Covariance of the smoothness inducing prior (penalizing first 
//'          differences in beta)
//' @param ind Vector of participants indicators corresponding to the observations.
//' @param sigma2_y Residual variance.
//' @param sigma2_us Smoothness inducing parameter for the random effects.
//' @param sigma2_ua Random effects variance. 
//' @return 
// [[Rcpp::export]]
Rcpp::List sample_reff(arma::vec res, arma::mat B, arma::mat P, arma::uvec ind, 
                       double sigma2_y, double sigma2_us, double sigma2_ua){
  
  arma::uvec un_ind = arma::unique(ind);
  unsigned int n_ind = un_ind.n_elem;
  unsigned int K = B.n_cols;
  unsigned int n = res.n_elem;
  
  arma::vec y_temp;
  arma::mat beta_ui(n_ind, K, arma::fill::zeros);
  arma::vec B_beta_u(n, arma::fill::zeros);
  arma::uvec idx_i; 
  
  arma::mat I_K(K, K, arma::fill::eye);
  arma::mat cov_post;
  arma::vec m_post;
  
  for (unsigned int i = 1; i <= n_ind; i++){
    idx_i = arma::find(ind == i);
    y_temp = res.elem(idx_i);
    
    cov_post = arma::inv_sympd(B.rows(idx_i).t() * B.rows(idx_i)/sigma2_y +
      P/sigma2_us + I_K/sigma2_ua);
    m_post = cov_post * B.rows(idx_i).t() * y_temp/sigma2_y;
    
    beta_ui.row(i - 1) = arma::mvnrnd(m_post, cov_post).t();
    B_beta_u.elem(idx_i) = B.rows(idx_i) * beta_ui.row(i - 1).t();
  }
  
  return Rcpp::List::create(
    Rcpp::Named("B_beta_u")=B_beta_u,
    Rcpp::Named("beta_ui")=beta_ui
  );
  
}


// [[Rcpp::export]]
arma::uvec match(arma::uvec large_vec, arma::uvec small_vec){
  
  unsigned int N = large_vec.n_elem;
  unsigned int n = small_vec.n_elem;
  arma::uvec out(N, arma::fill::zeros);
  
  for (unsigned int i = 0; i < n; i++){
    out = (out) || (large_vec == small_vec(i));
  }
  return (out);
}


// [[Rcpp::export]]
arma::uvec find_cov(arma::umat X, arma::field<arma::umat> Z_old,
                    arma::uvec cl_cov, unsigned int k){
  
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  arma::uvec idx_cov(n, arma::fill::ones);
  
  for (unsigned int j = 0; j < p; j++){
    idx_cov = idx_cov && match(X.col(j), arma::find(Z_old(j).col(k - 1) == cl_cov[j]) + 1);
  }
  
  return(idx_cov);
}


// // [[Rcpp::export]]
// arma::uvec find_cov_2(arma::field<arma::umat> Z_old, arma::uvec cl_cov, 
//                       unsigned int k){
//   
//   unsigned int p = Z_old.n_elem;
//   arma::uvec out(p, arma::fill::ones);
//   
//   for (unsigned int j = 0; j < p; j++){
//     out(j) = arma::find(Z_old(j).col(k - 1) == cl_cov[j]) + 1;
//   }
//   
//   return(out);
// }

// // [[Rcpp::export]]
// void prova_list(arma::field<arma::mat> Z){
//   
//   unsigned int n = Z.n_elem;
//   arma::mat Z_j;
//   for (unsigned int j = 0; j < n; j++){
//     Z_j = Z(j);
//     Rcpp::Rcout << "prova: " << Z_j << std::endl;
//   }
// }



// // [[Rcpp::export]]
// double log_lik_z(arma::vec y, arma::umat X, arma::field<arma::umat> Z,
//                  double sigma_2_y, double mu_0, double sigma_2_mu,
//                  arma::umat C_k, unsigned int k){
//   
//   unsigned int n = y.n_elem;
//   double R_i, R_2_i, n_h;
//   double out = 0.0;
//   arma::uvec idx_cov(n, arma::fill::ones);
//   unsigned int n_un = C_k.n_rows;
//   
//   for (unsigned int h = 0; h < n_un; h++){
//     idx_cov = find_cov(X, Z, C_k.row(h).t(), k);
//     n_h = arma::accu(idx_cov);
//     R_i = arma::accu(y.elem(arma::find(idx_cov == 1)));
//     R_2_i = arma::accu(arma::pow(y.elem(arma::find(idx_cov == 1)), 2.0));
//     
//     out += -0.5 * n_h * log(2.0 * arma::datum::pi * sigma_2_y) +
//       0.5 * (log(sigma_2_y) - log(sigma_2_mu * n_h + sigma_2_y)) -
//       0.5 / (sigma_2_y * sigma_2_mu) * (sigma_2_mu * R_2_i + sigma_2_y * std::pow(mu_0, 2.0) - std::pow(sigma_2_mu * R_i + sigma_2_y * mu_0, 2.0) / (sigma_2_mu * n_h + sigma_2_y));
//     
//   }
//   return (out);
// }


// [[Rcpp::export]]
double log_lik_z(arma::vec y, arma::umat X, arma::field<arma::umat> Z,
                 double sigma_2_y, double mu_0, double sigma_2_mu,
                 arma::umat C_k, arma::uvec z_sl, unsigned int k){
  
  unsigned int n = y.n_elem;
  double R_i, R_2_i, n_h;
  double out = 0.0;
  arma::uvec idx_cov(n, arma::fill::ones);
  arma::umat part_el;
  // unsigned int n_un = C_k.n_rows;
  arma::uvec z_sl_un = arma::unique(z_sl);
  unsigned int n_un = z_sl_un.n_elem;
  
  for (unsigned int h = 0; h < n_un; h++){
    part_el = C_k.rows(arma::find(z_sl == z_sl_un[h]));
    idx_cov = find_cov(X, Z, part_el.row(0).t(), k);
    
    for (unsigned int hh = 1; hh < part_el.n_rows; hh++){
      idx_cov += find_cov(X, Z, part_el.row(hh).t(), k);
    }
    
    n_h = arma::accu(idx_cov);
    R_i = arma::accu(y.elem(arma::find(idx_cov == 1)));
    R_2_i = arma::accu(arma::pow(y.elem(arma::find(idx_cov == 1)), 2.0));
    
    out += -0.5 * n_h * log(2.0 * arma::datum::pi * sigma_2_y) +
      0.5 * (log(sigma_2_y) - log(sigma_2_mu * n_h + sigma_2_y)) -
      0.5 / (sigma_2_y * sigma_2_mu) * (sigma_2_mu * R_2_i + sigma_2_y * std::pow(mu_0, 2.0) - std::pow(sigma_2_mu * R_i + sigma_2_y * mu_0, 2.0) / (sigma_2_mu * n_h + sigma_2_y));
    
  }
  
  return (out);
}
