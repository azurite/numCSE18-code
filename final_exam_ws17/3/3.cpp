#include <iostream>
#include <Eigen/Dense>
#include "gaussQuad.hpp"


using namespace Eigen;


// routine: evalF
// " Computes F(\alpha) "
// (in)  qr: quad nodes and weights
// (in)  \alpha: model parameters
// (in)  u: vector of moments
// (out) F(\alpha)
VectorXd evalF(const QuadRule& qr, const VectorXd& alpha, const VectorXd& u) {

  int m = alpha.size();
  VectorXd F = VectorXd::Random(m);

  //TODO: Implement here

  return F;
}


// routine: evalJ
// " Computes the Jacobian of $F(\alpha)$ "
// (in)  qr: quad nodes and weights
// (in)  alpha: model parameters
// (out) jacobian: Jacobian
MatrixXd evalJ(const QuadRule& qr, const VectorXd& alpha, const VectorXd& u) {

  int m = alpha.size();
  MatrixXd jacobian(m,m);

  //TODO: Implement here

  return jacobian;
}


// routine: newtonMethod
// " Solves the non-linear system using Newton method "
// (in)  qr: quad nodes and weights
// (in)  u: vector of moments
// (in)  atol: absolute tolerance
// (in)  rtol: relative tolerance
// (in)  maxItr: maximum iterations
// (in/out) alpha\_: model parameters
void newtonMethod(const QuadRule& qr, const VectorXd& u, const double atol, const double rtol, const int maxItr, VectorXd& alpha_) {

  //TODO: Implement here
	
}


// -- DO NOT CHANGE THIS ROUTINE --
// routine: testEvalF
// " test evalF computation "
// (in) alpha: model parameters
// (in) u: vector of moments
// (in) nQuad: number of quad nodes and weights
void testEvalF(const VectorXd& alpha, const VectorXd& u, const int nQuad) {

  QuadRule qr;
  gaussQuad(nQuad, qr);

  VectorXd F = evalF(qr, alpha, u);
  std::cout << "Norm of F: " << F.norm() << ", should be small" << std::endl;
  
}


// routine: main
int main() {

  int nQ = 40; // number of Gauss points

  // -- DO NOT CHANGE THIS TEST SNIPPET --
  // test evalF
  // optimal model parameters corresponding to a moment vector are given
  VectorXd u_test(3), alpha_test(3);
  u_test << 1, 0, 1;
  alpha_test << -0.9189385332, 0, -0.5;
  testEvalF (alpha_test, u_test, nQ);

  // subproblem (h)
  VectorXd u(3);
  u << 1, -1, 2;
  
  //TODO: run your Newton method implementation
  
}


// END OF FILE
