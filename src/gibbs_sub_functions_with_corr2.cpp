// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <unistd.h>

#define pi 3.1415926

using namespace Rcpp;
using namespace arma;



//' @keywords internal
//' 
mat row_wise_upper_tri_to_full(const vec& Beta_upper_tri, const int& K) {
  mat B = zeros(K,K);
  int n = 0;
  for (int i = 0; i<K; i++) {
    for (int j = i+1; j < K; j++) {
      B(i,j) = Beta_upper_tri(n);
      n++;
    }
  }
  return B;
}


//' @keywords internal
//' 
mat LTstack(const vec& v) {
  const int K = v.n_elem;
  //I don't know if the cast to int is needed, note that K(K-1) is always even
  mat M = zeros( (int)(K*(K-1)/2) , K);
  int rowstart = 0;
  for (int j = 0; j < K-1; j++) {
    int len = (K - j) - 1;
    vec subvec = v.tail(len);
    M( rowstart, j, size(len, 1) ) = subvec;
    rowstart +=  len;
  }
  return M;
}


// //dont: [[Rcpp::export]]
mat LTstackMatrix(const mat& A){
  const int K = A.n_rows;
  const int P = A.n_cols;
  mat A_tilde = mat( K*(K-1)/2, K*P);
  for (int j=0; j < P; j++) {
    A_tilde.cols(j*K, (j+1)*K - 1) = LTstack(A.col(j));
  }
  return A_tilde;
}


//' @keywords internal
//' 
double dmvn1(const vec& x, const vec& mean, const mat& sigma) {
  auto p = x.n_elem;
  //cout << "p " << endl << p << endl;
  
  auto R = chol(sigma);
  //cout << "R " << endl << R << endl;
  
  //cout << "x-m " << endl << x - mean << endl;
  //cout << trimatl(R.t()) << endl;
  vec tmp = solve(trimatl(R.t()), x - mean);
  //cout << "tmp " << endl << tmp << endl;
  
  auto rss = sum(square(tmp));
  //cout << "rss " << rss << endl;
  
  auto logretval = -sum(log(mat(R).diag())) - 0.5 * p * log(2 * pi) - 0.5 * rss;
  return exp(logretval);
}

//' @keywords internal
//' 
vec rmvn(const vec& M, const mat& C) {
  return mvnrnd(M, C);
}


//' @keywords internal
//' 
mat sigma_z(const vec& sigma1, const vec& sigma0, const int& z, const int& K){
  vec out = vec(K);
  for(int i=0; i < K; i++){
    if( z & (1<<i) ){
      out(i) = sigma1(i);
    }
    else{
      out(i) = sigma0(i);
    }
  }
  return(square(diagmat(out)));
}


//' @keywords internal
//' 
void sample_Z_and_A_CPP(const mat& Y, const mat& SD, const mat& B,
                        const vec& sigma0, const vec& sigma1,
                        const vec& priors,
                        mat& A_out,
                        umat& Z_out){
  
  int K = Y.n_rows;
  int P = Y.n_cols;
  umat I;
  I.eye(size(B));
  
  std::vector<mat> Sigma_zs;
  std::vector<mat> Var_first_terms;
  
  for (int z=0; z < (1<<K); z++){
    Sigma_zs.push_back(sigma_z(sigma1, sigma0, z, K));
    Var_first_terms.push_back((I + B) *  Sigma_zs.back() * (I + B).t());
  }
  
  vec likelihood = vec(1<<K);
  
  for (int j=0; j<P; j++){
    mat Delta = diagmat(square(SD.col(j)));
    for (int z=0; z < (1<<K); z++){
      likelihood(z) = dmvn1(Y.col(j), zeros(K), Var_first_terms[z] + Delta);
    }
    //% is element-wise mult
    //multiply by prior (in place) to get posterior
    likelihood %= priors;
    //here we sample from 1...K weighted by likelihood:
    double selector = randu() * sum(likelihood);
    double cum_prob = 0;
    int z=0;
    for (; z < (1<<K); z++){
      cum_prob += likelihood(z);
      if (selector <= cum_prob){
        //we have selected this choice of z, assign it to Z out
        for(int i=0; i < K; i++){
          Z_out(i, j) = (z>>i) & 1;
        }
        break;
      }
    }
    //Compute variance
    mat C =  inv_sympd((I + B).t() *  inv(diagmat(Delta)) * (I + B) + inv(diagmat(Sigma_zs[z])));
    //Compute mean
    vec m = C *  (I + B).t() * inv(diagmat(Delta)) * Y.col(j);
    A_out.col(j)  = mvnrnd(m, C);
  }
}

//' @keywords internal
//' 
void sample_Z_and_A_CPP_with_corr(const mat& Y,
                                  mat& Lambda,
                                  mat& Lambda_inv,
                                  const mat& B,
                                  const vec& sigma0, const vec& sigma1,
                                  const vec& priors,
                                  mat& A_out,
                                  umat& Z_out){
  
  int K = Y.n_rows;
  int P = Y.n_cols;
  umat I;
  I.eye(size(B));
  
  std::vector<mat> Sigma_zs;
  std::vector<mat> Var_first_terms;
  
  for (int z=0; z < (1<<K); z++){
    Sigma_zs.push_back(sigma_z(sigma1, sigma0, z, K));
    Var_first_terms.push_back((I + B) *  Sigma_zs.back() * (I + B).t());
  }
  
  vec likelihood = vec(1<<K);
  
  for (int j=0; j<P; j++){
    for (int z=0; z < (1<<K); z++){
      likelihood(z) = dmvn1(Y.col(j), zeros(K), Var_first_terms[z] +
        Lambda.submat((K * j ), (K*j), (K*j +K-1 ), (K*j +K-1)));
    }
    //% is element-wise mult
    //multiply by prior (in place) to get posterior
    likelihood %= priors;
    //here we sample from 1...K weighted by likelihood:
    double selector = randu() * sum(likelihood);
    double cum_prob = 0;
    int z=0;
    for (; z < (1<<K); z++){
      cum_prob += likelihood(z);
      if (selector <= cum_prob){
        //we have selected this choice of z, assign it to Z out
        for(int i=0; i < K; i++){
          Z_out(i, j) = (z>>i) & 1;
        }
        break;
      }
    }
    //Compute variance
    //cout << j;
    
    
    ////  //////  //////  This is the line where it gets the error ///// ///// //// ////
    ////  //////  //////  This is the line where it gets the error ///// ///// //// ////
    ////  //////  //////  This is the line where it gets the error ///// ///// //// ////
    ////  //////  //////  This is the line where it gets the error ///// ///// //// ////
    ////  //////  //////  This is the line where it gets the error ///// ///// //// ////
    
    
    mat C =  inv_sympd((I + B).t() *   Lambda_inv.submat((K * j ), (K*j ), (K*j +K-1 ), (K*j + K-1))
                         * (I + B) + inv(diagmat(Sigma_zs[z])));//this one
    //Compute mean
    vec m = C *  (I + B).t() *  Lambda_inv.submat((K * j ), (K*j), (K*j +K-1), (K*j +K - 1)) * Y.col(j);
    A_out.col(j)  = mvnrnd(m, C);
  }
}




//' @keywords internal
//' 
inline double square(double x){return x*x;}



//' @keywords internal
//' 
mat sample_B(const mat& Y, const mat& SD, mat& A, double& sigma){
  //this function is actually slower than R but may be useful if we implement the whole sampling loop in C++
  int K = Y.n_rows;
  mat B;
  B.zeros(K, K);
  for(int i=0; i < (K-1); i++){
    mat A_tilde = A.rows(i+1,K-1);
    mat Delta_inv = inv(diagmat(square(SD.row(i))));
    mat C = inv_sympd(A_tilde * diagmat(Delta_inv) * A_tilde.t() + (1/(sigma*sigma)) * eye(K-i-1, K-i-1));
    //mat first = C * A_tilde * diagmat(Delta_inv);
    //vec second = trans(Y.row(i) - A.row(i));
    vec m = C * A_tilde * diagmat(Delta_inv) * trans(Y.row(i) - A.row(i));
    //std::cout << "m " << std::endl << m << std::endl;
    //fill up B transposed for efficiency and traspose return value
    //(probably negligible gain but takes more time to check)
    B(span(i+1, K-1), i) = mvnrnd(m, C);
  }
  return B.t();
}

//' @keywords internal
//' 
mat sample_B_with_corr(const mat& Y, const sp_mat& Lambda_inv, mat& A, double& sigma,
                       int& n){
  //this function is actually slower than R but may be useful if we implement the whole sampling loop in C++
  int K = Y.n_rows;
  int KC2 = K * (K - 1) / 2;
  mat A_tilde = LTstackMatrix(A);
  mat C = inv_sympd(A_tilde * Lambda_inv * A_tilde.t() 
                      + (1/square(sigma)) * eye(KC2, KC2));
  vec m = C * A_tilde * diagmat(Lambda_inv) * vectorise(Y - A);
  vec B_ut = mvnrnd(m, C) ;
  mat B = row_wise_upper_tri_to_full(B_ut, K);
  
  

  
  return B;
}




//' @keywords internal
//' 
mat sample_B_no_corr(const mat& Y, const mat& SD, mat& A, double& sigma){
  //this function is actually slower than R but may be useful if we implement the whole sampling loop in C++
  int K = Y.n_rows;
  int P = Y.n_cols;
  int KC2 = K * (K - 1) / 2;
  mat A_tilde = LTstackMatrix(A);
  
  sp_mat Lambda_inv = sp_mat(K*P, K*P);
  
  for (int j = 0; j < P; j++){
    Lambda_inv(K*j, K*j, size(K, K)) = square(inv(diagmat((SD.col(j)))));
  }
  
  mat C = inv_sympd(A_tilde * Lambda_inv * A_tilde.t() + (1/square(sigma)) * eye(KC2, KC2));
  vec m = C * A_tilde * diagmat(Lambda_inv) * vectorise(Y - A);
  vec B_ut = mvnrnd(m, C);
  mat B = row_wise_upper_tri_to_full(B_ut, K);
  
  return B;
}


//' @keywords internal
//' 
mat sample_B_with_corr_2(const mat& Y, const mat& SD, mat& A, double& sigma,
                         const mat& corr){
  //this function is actually slower than R but may be useful if we implement the whole sampling loop in C++
  int K = Y.n_rows;
  int P = Y.n_cols;
  int KC2 = K * (K - 1) / 2;
  mat A_tilde = LTstackMatrix(A);
  
  sp_mat Lambda_inv = sp_mat(K*P, K*P);
  
  for (int j = 0; j < P; j++){
    Lambda_inv(K*j, K*j, size(K, K)) =
      inv_sympd(diagmat(SD.col(j)) * corr * diagmat(SD.col(j)));
  }
  
  mat C = inv_sympd(A_tilde * Lambda_inv * A_tilde.t() + (1/square(sigma)) * eye(KC2, KC2));
  vec m = C * A_tilde * diagmat(Lambda_inv) * vectorise(Y - A);
  vec B_ut = mvnrnd(m, C);
  mat B = row_wise_upper_tri_to_full(B_ut, K);
  
  return B;
}

//' @keywords internal
//' 
mat sample_B_with_corr_3(const mat& Y, const mat& SD, mat& A, double& sigma,
                         int& n,
                         const mat& Lambda_inv){
  //this function is actually slower than R but may be useful if we implement the whole sampling loop in C++
  int K = Y.n_rows;
  int KC2 = K * (K - 1) / 2;
  mat A_tilde = LTstackMatrix(A);
  
  mat C = inv_sympd(A_tilde * Lambda_inv * A_tilde.t() + 
    (1/square(sigma)) * eye(KC2, KC2));
  vec m = C * A_tilde * Lambda_inv* vectorise(Y - A);
  vec B_ut = mvnrnd(m, C);
  mat B = row_wise_upper_tri_to_full(B_ut, K);
  
  
  return B;
}



//' Perform Gibbs sampling for mediation MR with correlation
//' 
//' Model is that Y = (I+B)A + noise
//' @param Y KxP arma::matrix of GWAS marginal associations
//' @param Sd_hat KxP arma::matrix of GWAS std errors
//' @param trait_corr KxK arma::matrix of correlation of noises.
//' @param N number of samples to obtain
//' @param B initial value of B arma::matrix
//' @param sigma Initial value of prior std deviation for values of B
//' @param sigma1 K-arma::vector of Initial values of prior slab SDs for each trait's pleiotropy
//' @param sigma0 K-arma::vector of Initial values of prior splike SDs for each trait's pleiotropy
//' @param p K-arma::vector of initial values for p, which is slab probability for mixture
//' @param A KxP arma::matrix of initial values for A (ignored)
//' @param Z KxP arma::matrix of initial values for Z (ignored), which is latent varaible indicating if A[i,j] comes from slab
//' @param alpha_B Shape parameter of InvGamma prior for sigma (for B)
//' @param beta_B Scale parameter of InvGamma prior for sigma (for B)
//' @param alpha_0 K-arma::vector Shape parameter of InvGamma prior for sigma0 (Spike)
//' @param alpha_1 K-arma::vector Shape parameter of InvGamma prior for sigma0 (Slab)
//' @param beta_0 K-arma::vector Scale parameter of InvGamma prior for sigma0 (Spike)
//' @param beta_1 K-arma::vector Scale parameter of InvGamma prior for sigma0 (Slab)
//' @param a K-arma::vector First parameter of beta prior for p
//' @param b K-arma::vector Second parameter of beta prior for p
//' @param Lambda3 Precomputed covariance matrix
//' @param Lambda_inv3 Precomputed inverse covariance matrix
//' @export
// [[Rcpp::export]]
List gibbs_sampler_with_corr(arma::mat& Y,
                             arma::mat& Sd_hat,
                             arma::mat& trait_corr,
                             int& N,
                             arma::mat& B,
                             double& sigma,
                             arma::vec& sigma1,
                             arma::vec& sigma0,
                             arma::vec& p,
                             arma::mat& A,
                             arma::umat& Z,
                             double& alpha_B,
                             double& beta_B,
                             arma::vec& alpha_0,
                             arma::vec& alpha_1,
                             arma::vec& beta_0,
                             arma::vec& beta_1,
                             arma::vec& a,
                             arma::vec& b,
                             arma::mat& Lambda3,
                             arma::mat& Lambda_inv3){
  
  
  
  
  const auto P = Y.n_cols;
  const auto K = Y.n_rows;
  
  
  cube B_samples = cube(K, K, N);
  vec sigma_samples = vec(N);
  mat sigma1_samples = mat(K, N);
  mat sigma0_samples = mat(K, N);
  mat p_samples = mat(K, N);
  ucube Z_samples;
  Z_samples.set_size(K, P, N);
  cube A_samples = cube(K, P, N);
  
  //store initialization:
  int n=0;;
  B_samples.slice(n) = B;
  sigma_samples(n) = sigma;
  sigma1_samples.col(n) = sigma1;
  sigma0_samples.col(n) = sigma0;
  p_samples.col(n) = p;
  Z_samples.slice(n) = Z;
  A_samples.slice(n) = A;
  
  
  //this is just for debugging:
  Environment env = Environment::global_env();
  
  
  //generate correlated noise covariance matrices
  std::vector<mat> Lambdas;
  std::vector<mat> Lambda_invs;
  sp_mat Lambda_inv = sp_mat(K*P, K*P);
  
  
  //iterate for each sample:
  for(n=1; n<N; n++){
    
    
    //First update z
    
    //Compute priors
    vec priors =  ones(1<<K);
    //vector of priors for each choice of z
    for (int z=0; z < (1<<K); z++){
      //prior is equal to product of Bernoulli probabilities with
      //parameters p_i
      for (int i=0; i<K; i++){
        priors(z) *= (z & (1<<i)) ? p(i) : (1 - p(i));
      }
    }
    
    
    sample_Z_and_A_CPP_with_corr(Y, Lambda3, Lambda_inv3, B, sigma0, sigma1, priors, A, Z);
    
    
    //Then update p
    for(int i=0; i<K; i++){
      //Sample two Gamma rvs to get Beta rv: B = G_1 / (G_1 + G_2)
      //See en.wikipedia.org/wiki/Beta_distribution#Computational_methods (accessed 3/31/22)
      double X = randg(distr_param(a(i) + sum(Z.row(i)), 1.0));
      double Y = randg(distr_param(b(i) + P - sum(Z.row(i)), 1.0));
      p(i) =  X / (X+Y);
    }
    
    
    //And update sigma1 and sigma0
    for(int i=0; i<K; i++){
      //Compute the parameters for sigma1 (% is element wise mult)
      double alpha1_par = alpha_1(i) + sum(Z.row(i))/2;
      double beta1_par = beta_1(i) +  sum( square(A.row(i) % Z.row(i) ))/2;
      
      //Compute the parameters for sigma0
      double alpha0_par =  alpha_0(i) + (P - sum(Z.row(i)) )/2;
      double beta0_par= beta_0(i) +  sum ( square(A.row(i) % ( 1 - Z.row(i))) )/2;
      
      //Update sigma1 and sigma0
      sigma1(i) = sqrt(1/randg(distr_param(alpha1_par, 1/beta1_par)));
      sigma0(i) = sqrt(1/randg(distr_param(alpha0_par, 1/beta0_par)));
    }
    
    
    ////Then update B
    B = sample_B_with_corr_3(Y, Sd_hat, A, sigma, n, Lambda_inv3);
    
    //Finally update sigma
    sigma = sqrt(1/randg(distr_param(alpha_B + K * (K - 1)/2, 1/(beta_B + square(norm(B, "fro"))/2))));
    
    //Now store samples:
    B_samples.slice(n) = B;
    sigma_samples(n) = sigma;
    sigma1_samples.col(n) = sigma1;
    sigma0_samples.col(n) = sigma0;
    p_samples.col(n) = p;
    Z_samples.slice(n) = Z;
    A_samples.slice(n) = A;
    
    
  }
  
  return Rcpp::List::create(Rcpp::Named("A") = A_samples,
                            Rcpp::Named("Z") = Z_samples,
                            Rcpp::Named("sigma1") = sigma1_samples,
                            Rcpp::Named("sigma0") = sigma0_samples,
                            Rcpp::Named("sigma") = sigma_samples,
                            Rcpp::Named("B") = B_samples,
                            Rcpp::Named("p") = p_samples);
  
}



//' Perform Gibbs sampling for mediation MR without correlation
//' 
//' Model is that Y = (I+B)A + noise
//' @param Y KxP arma::matrix of GWAS marginal associations
//' @param Sd_hat KxP arma::matrix of GWAS std errors
//' @param N number of samples to obtain
//' @param B initial value of B arma::matrix
//' @param sigma Initial value of prior std deviation for values of B
//' @param sigma1 K-arma::vector of Initial values of prior slab SDs for each trait's pleiotropy
//' @param sigma0 K-arma::vector of Initial values of prior splike SDs for each trait's pleiotropy
//' @param p K-arma::vector of initial values for p, which is slab probability for mixture
//' @param A KxP arma::matrix of initial values for A (ignored)
//' @param Z KxP arma::matrix of initial values for Z (ignored), which is latent varaible indicating if A[i,j] comes from slab
//' @param alpha_B Shape parameter of InvGamma prior for sigma (for B)
//' @param beta_B Scale parameter of InvGamma prior for sigma (for B)
//' @param alpha_0 K-arma::vector Shape parameter of InvGamma prior for sigma0 (Spike)
//' @param alpha_1 K-arma::vector Shape parameter of InvGamma prior for sigma0 (Slab)
//' @param beta_0 K-arma::vector Scale parameter of InvGamma prior for sigma0 (Spike)
//' @param beta_1 K-arma::vector Scale parameter of InvGamma prior for sigma0 (Slab)
//' @param a K-arma::vector First parameter of beta prior for p
//' @param b K-arma::vector Second parameter of beta prior for p
//' 
//' @export
// [[Rcpp::export]]
List gibbs_sampler(arma::mat& Y,
                   arma::mat& Sd_hat,
                   int& N,
                   arma::mat& B,
                   double& sigma,
                   arma::vec& sigma1,
                   arma::vec& sigma0,
                   arma::vec& p,
                   arma::mat& A,
                   arma::umat& Z,
                   double& alpha_B,
                   double& beta_B,
                   arma::vec& alpha_0,
                   arma::vec& alpha_1,
                   arma::vec& beta_0,
                   arma::vec& beta_1,
                   arma::vec& a,
                   arma::vec& b){
  
  
  
  
  const auto P = Y.n_cols;
  const auto K = Y.n_rows;
  
  
  cube B_samples = cube(K, K, N);
  vec sigma_samples = vec(N);
  mat sigma1_samples = mat(K, N);
  mat sigma0_samples = mat(K, N);
  mat p_samples = mat(K, N);
  ucube Z_samples;
  Z_samples.set_size(K, P, N);
  cube A_samples = cube(K, P, N);
  
  //store initialization:
  int n=0;;
  B_samples.slice(n) = B;
  sigma_samples(n) = sigma;
  sigma1_samples.col(n) = sigma1;
  sigma0_samples.col(n) = sigma0;
  p_samples.col(n) = p;
  Z_samples.slice(n) = Z;
  A_samples.slice(n) = A;
  
  
  //Get all possible Z vectors in order
  //candidates = permutations(2, K, 0:1, repeats.allowed = TRUE)
  
  //iterate for each sample:
  for(n=1; n<N; n++){
    
    //First update z
    
    //Compute priors
    vec priors =  ones(1<<K);
    //vector of priors for each choice of z
    for (int z=0; z < (1<<K); z++){
      //prior is equal to product of Bernoulli probabilities with
      //parameters p_i
      for (int i=0; i<K; i++){
        priors(z) *= (z & (1<<i)) ? p(i) : (1 - p(i));
      }
    }
    
    
    sample_Z_and_A_CPP(Y, Sd_hat, B, sigma0, sigma1, priors, A, Z);
    
    
    //Then update p
    for(int i=0; i<K; i++){
      //Sample two Gamma rvs to get Beta rv: B = G_1 / (G_1 + G_2)
      //See en.wikipedia.org/wiki/Beta_distribution#Computational_methods
      double X = randg(distr_param(a(i) + sum(Z.row(i)), 1.0));
      double Y = randg(distr_param(b(i) + P - sum(Z.row(i)), 1.0));
      p(i) =  X / (X+Y);
    }
    
    
    //And update sigma1 and sigma0
    for(int i=0; i<K; i++){
      //Compute the parameters for sigma1 (% is element wise mult)
      double alpha1_par = alpha_1(i) + sum(Z.row(i))/2;
      double beta1_par = beta_1(i) +  sum( square(A.row(i) % Z.row(i) ))/2;
      
      //Compute the parameters for sigma0
      double alpha0_par =  alpha_0(i) + (P - sum(Z.row(i)) )/2;
      double beta0_par= beta_0(i) +  sum ( square(A.row(i) % ( 1 - Z.row(i))) )/2;
      
      //Update sigma1 and sigma0
      sigma1(i) = sqrt(1/randg(distr_param(alpha1_par, 1/beta1_par)));
      sigma0(i) = sqrt(1/randg(distr_param(alpha0_par, 1/beta0_par)));
    }
    
    
    ////Then update B
    B = sample_B(Y, Sd_hat, A, sigma);
    
    //Finally update sigma
    sigma = sqrt(1/randg(distr_param(alpha_B + K * (K - 1)/2, 1/(beta_B + square(norm(B, "fro"))/2))));
    
    //Now store samples:
    B_samples.slice(n) = B;
    sigma_samples(n) = sigma;
    sigma1_samples.col(n) = sigma1;
    sigma0_samples.col(n) = sigma0;
    p_samples.col(n) = p;
    Z_samples.slice(n) = Z;
    A_samples.slice(n) = A;
    
    
  }
  
  return Rcpp::List::create(Rcpp::Named("A") = A_samples,
                            Rcpp::Named("Z") = Z_samples,
                            Rcpp::Named("sigma1") = sigma1_samples,
                            Rcpp::Named("sigma0") = sigma0_samples,
                            Rcpp::Named("sigma") = sigma_samples,
                            Rcpp::Named("B") = B_samples,
                            Rcpp::Named("p") = p_samples);
  
}



