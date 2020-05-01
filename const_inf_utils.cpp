// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <queue> 
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector makebridges(NumericVector tau, int m, int p) {
  int nt = tau.size();
  NumericVector M(nt), S(nt);
  NumericVector B(nt * p * m);
  NumericVector Z = rnorm(nt * p * m);
  IntegerVector dims = IntegerVector::create(nt, p, m);
  /* M and S are conditional mean and standard deviation vectors */
  if (tau[0] != 0.0) {
    S[0] = sqrt(tau[0] * (1 - tau[0]));
  }
  for(int i = 1; i < nt; i++) {
    M[i] = (1 - tau[i]) / (1 - tau[i - 1]);
    S[i] = sqrt((1 - tau[i]) * (tau[i] - tau[i - 1]) / (1 - tau[i - 1]));
  }
  for (int k = 0; k < m; k++) {
    for(int j = 0; j < p; j++) {
      B[(nt * p) * k + nt * j] = S[0] * Z[(nt * p) * k + nt * j];
      for (int i = 1; i < nt; i++) {
        B[(nt * p) * k + nt * j + i] = M[i] * B[(nt * p) * k + nt * j + i - 1] + 
          S[i] * Z[(nt * p) * k + nt * j + i];
      }
    }
  }
  B.attr("dim") = dims;
  return B;
}

// Array multiplication thing that I couldn't figure out how to vectorize in R.
// Made for scaling a vector of Brownian bridges by covariance matrices to get 
// processes that have the same distribution as inference processes.
// [[Rcpp::export]]
arma::cube scalebridge(arma::cube bridge, arma::cube scl) {
  arma::cube ref;
  ref.copy_size(bridge);
  for (arma::uword i = 0; i < bridge.n_slices; ++i) {
    for (arma::uword j = 0; j < bridge.n_rows; ++j) {
      ref.slice(i).row(j) = bridge.slice(i).row(j) * scl.slice(j);
    }
  }
  return(ref);
}

// This is for calculating the mixture of tied down squared Bessel processes.
// It uses the matrix flavor of the identity a^2 - b^2 = (a+b)(a-b).
// proc1 is supposed to be the gaussian process corresponding to the alternative
// hypothesis, proc0 is the one for the null.  sig is an array of covariance
// matrices with each slice corresponding to one quantile level.

// [[Rcpp::export]]
arma::mat mixproc(arma::cube proc1, arma::cube proc0, arma::cube sig) {
  int rep = proc1.n_slices;
  int nt = sig.n_slices;
  int p = proc1.n_cols;
  arma::mat sinv(p, p);
  arma::rowvec q1(p);
  arma::rowvec q2(p);
  arma::mat cbp(nt, rep);
  for (arma::uword i = 0; i < rep; ++i) {
    for (arma::uword j = 0; j < nt; ++j) {
      sinv = inv_sympd(sig.slice(j));
      q1 = proc1.slice(i).row(j) + proc0.slice(i).row(j);
      q2 = proc1.slice(i).row(j) - proc0.slice(i).row(j);
      cbp.row(j).col(i) = q1 * sinv * q2.t(); /* since qs are rowvecs */
    }
  }
  return(cbp);
}

// Calculate a square root matrix and the square root of the inverse of 
// the matrix too.
// Assumed that the matrices are all symmetric and positive definite
// [[Rcpp::export]]
List matsqrts(arma::mat A) {
  arma::mat A12 = sqrtmat_sympd(A); 
  arma::mat Ainv = inv_sympd(A);
  arma::mat Am12 = sqrtmat_sympd(Ainv);
  return(List::create(_["A12"] = wrap(A12), _["Am12"] = wrap(Am12)));
}

// code stolen from a stackexchange thread, "indices of the k largest 
// elements in an unsorted length n array" and modified a bit
// https://stackoverflow.com/questions/14902876/indices-of-the-k-largest-elements-in-an-unsorted-length-n-array

// [[Rcpp::export]]
IntegerVector smallind(NumericVector test, int num) {
  std::priority_queue< std::pair<double, int>, 
                        std::vector< std::pair<double, int> >, 
                        std::less <std::pair<double, int> > > q;
  for (int i = 0; i < test.size(); ++i) {
    if(q.size() < num)
      q.push(std::pair<double, int>(test[i], i));
    else if(q.top().first > test[i]){
      q.pop();
      q.push(std::pair<double, int>(test[i], i));
    }
  }
  int k = q.size();
  std::vector<int> res(k);
  for (int i = 0; i < k; ++i) {
    res[k - i - 1] = q.top().second + 1;
    q.pop();
  }
  return(wrap(res));
}

