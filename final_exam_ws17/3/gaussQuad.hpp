# include <Eigen/Dense>
# include <Eigen/Eigenvalues>


using namespace Eigen;


struct QuadRule {
  VectorXd nodes, weights;
};


// routine: gaussQuad
// " Computes the Gauss nodes and weights in the domain [-,1,1] "
// (in)  n: number of nodes 
// (out) qr: stores the computed nodes and weights 
void gaussQuad(const unsigned n, QuadRule& qr) {
  
  MatrixXd M = MatrixXd::Zero(n, n);
  for (unsigned i = 1; i < n; ++i) {
    const double b = i/std::sqrt(4.*i*i - 1.);
    M(i,i-1) = b; M(i-1,i) = b;
  }
  SelfAdjointEigenSolver<MatrixXd> eig(M);
  
  qr.nodes = eig.eigenvalues();
  qr.weights = 2*eig.eigenvectors().topRows<1>().array().pow(2);

}


// END OF FILE
