// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatcrossprod(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A.adjoint() * B;
  
  return Rcpp::wrap(C);
}
// [[Rcpp::export]]
SEXP eigenMapMattcrossprod(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B.adjoint();
  
  return Rcpp::wrap(C);
}



