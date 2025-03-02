// This script implement sketch kernel PCA with feature selection
// Date: 2024-12-31
// Revised log:


#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
//#include<boost/math/tools/minima.hpp>
#include <cstdlib>
#include <vector>
#include <algorithm>

#define INT_MIN1 (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;
// using boost::math::tools::brent_find_minima;


// Define global variables
//double X_count_ij, a_i, invLambda_j,Mu_x_ij;
/*
 * Auxiliary
 */
//' @keywords internal
 //' @noRd
 //'
// [[Rcpp::export]]
arma::sp_mat sketchfun_cpp( const int& m, const int& n,const std::string& type = "subGaussian"){
   mat S;
   if(type == "subGaussian"){
     S = randn(m, n);
   }else if(type == "ROS"){
     arma::mat E = arma::eye<arma::mat>(n, n);  // Identity matrix of size n x n
     uvec ind = arma::conv_to<arma::uvec>::from(randperm(n));
     uvec ind_vec = ind.subvec(0, m-1);
     // Generate a random boolean vector
     arma::vec tmp = arma::randu<arma::vec>(n);
     // Create a vector of -0.5 and set the elements corresponding to true in bool_vec to 0.5
     uvec idx1 = find(tmp>=0.5);
     tmp.fill(-1);
     tmp.elem(idx1).ones();
     S = sqrt(n/m)*E.rows(ind_vec) * (1.0/sqrt(n)*randn(n, n) % repmat(tmp/2.0,1,n)); // R;
   }else if(type == "subSampling"){
     arma::mat E = arma::eye<arma::mat>(n, n);  // Identity matrix of size n x n
     uvec ind = arma::conv_to<arma::uvec>::from(randperm(n));
     uvec ind_vec = ind.subvec(0, m-1);
     arma::mat Praw_ind = E.rows(ind_vec);
     S = sqrt(static_cast<double>(n) / m) * Praw_ind;
   }

   return arma::conv_to<arma::sp_mat>::from(S);
}


List irlbaCpp(const mat& X, const int& q){
  Rcpp::Environment irlba("package:irlba");
  Rcpp::Function f = irlba["irlba"];

  return f(Named("A") = X, Named("nv") = q);

}

// arma::mat eigenR(const mat& X, const int& d, const bool& sym=true){
//
//   Rcpp::Environment base("package:base");
//   Rcpp::Function f = base["eigen"];
//   Rcpp::List res =  f(Named("x") = X, Named("symmetric") = sym);
//   mat A = res["vectors"];
//   return A.cols(0, d-1); // type transfer
// }
// void sorteigen(const vec& eigvals, mat& eigvecs){
//   // sort the eigenvalues in decreasing order
//   // Method 1: faster version
//   uvec index = sort_index(eigvals, "descend");
//   // sorted_eigvals = eigvals(index);
//   eigvecs = eigvecs.cols(index);
// }
// [[Rcpp::export]]
arma::mat kernelMatrix_cpp(const Rcpp::Function& kernel, const arma::mat& x){
  Rcpp::Environment kernlab("package:kernlab");
  Rcpp::Function f = kernlab["kernelMatrix"];
  SEXP K = f(Named("kernel") = kernel, Named("x") = x);
  // mat K = f(Named("kernel") = kernel, Named("x") = x);
  return as<mat>(K); // type transfer
}

// [[Rcpp::export]]
Rcpp::List SKPCA_cpp(const arma::mat& X, const Rcpp::Function& kernel, const int& d,
                   const int& m, const std::string& Stype){

  int n = X.n_rows;
  sp_mat S = sketchfun_cpp(m, n, Stype);
  mat K = kernelMatrix_cpp(kernel, X);
  // K.row(0).print();
  mat SKS = S*K*S.t();
  // Rprintf("good1, \n");
  vec eigval;
  mat eigvec, eigvec2;
  svd(eigvec, eigval, eigvec2, SKS);
  // eig_sym( eigval, eigvec, SKS);
  // sorteigen(eigval, eigvec);
  // Rcpp::List eigX = eigenR(SKS);
  // mat eigvec = eigX["vectors"];
  // vec eigval= eigX["values"];
  // Rprintf("good2, \n");

  mat V = K * S.t()* eigvec * diagmat(1/sqrt(eigval));
  Rcpp::List svdV = irlbaCpp(V, d);
  // mat A1 = svdV["u"];
  // A = A1; A1.reset();
  // K = V * V.t();
  // Rprintf("good3, \n");
  List res = List::create(
    Rcpp::Named("A") = svdV["u"],
    Rcpp::Named("svals") = svdV["d"],
    Rcpp::Named("K") = V * V.t());
  return res;
}

// [[Rcpp::export]]
arma::mat update_w(const arma::mat& X, const arma::mat& A, const arma::mat& K,
                   const int& num_fea,const std::string& Ktype, const double& Kg_h){
  // update w using sorting method:
  int l, n = X.n_rows, p = X.n_cols;
  //vec  w(p, fill::zeros), y(p, fill::zeros);
  vec  w(p, fill::zeros);
  mat KA, wy(p, 2);
  if(Ktype == "linear"){
    mat XtA = X.t()*A;
    wy.col(1) = sum(XtA % XtA, 1); // p-vector
    // wy.t().print();
  }else if(Ktype == "gaussian"){
    for(l=0; l<p; ++l){
      KA = -2/Kg_h* (K % ((repmat(X.col(l), 1, n)-repmat(X.col(l).t(), n, 1)) % (repmat(X.col(l), 1, n)-repmat(X.col(l).t(), n, 1)))) * A;
      wy(l,1) = accu(KA % A);
      // wy.t().print();
    }
  }
  uvec index = sort_index(wy.col(1), "descend");
  w(index.subvec(0,num_fea-1)).fill(1.0);
  wy.col(0) = w;
  return wy;
}
// [[Rcpp::export]]
double objfun_Cpp(const arma::mat& X, const arma::mat& A){
  int n = X.n_rows;
  mat XA = X.t()*A;
  return accu(XA % XA);
}
