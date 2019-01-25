#include <iostream>

#include <Eigen/Dense>

Eigen::MatrixXd Jacobian(const Eigen::MatrixXd &A, Eigen::VectorXd &z) {
	int n = z.size() - 1;	
	int lambda = z(n);

	Eigen::MatrixXd J(n + 1, n + 1);
	J.topLeftCorner(n , n) = (A - z(n) * Eigen::MatrixXd::Identity(n, n));
	J.row(n) = -z.transpose();
	J.col(n) = -z;

	// overwrite value back to 0
	J(n, n) = 0;

	return J;
}

void eigNewton(const Eigen::MatrixXd &A, double tol, int maxItr, Eigen::VectorXd &z) {

	int n = z.size() - 1;
	Eigen::VectorXd F(n + 1);

	for(int i = 0; i < maxItr; i++) {
		Eigen::VectorXd x = z.head(n);
		F.head(n) = A * x - z(n) * x;
		F(n) = 1 - 0.5 * x.squaredNorm();
		
		if(F.squaredNorm() < tol) {
			return;		
		}

		z += Jacobian(A, z).fullPivLu().solve(-F);
	}
}

int main() {
	int n = 2;
	int m = n + 1;

	Eigen::VectorXd x = Eigen::VectorXd::Ones(n);
	x << 1.0, 0.0;
	double lambda = 0.0;

	Eigen::VectorXd z(m);
	z.head(n) = x;
	z(n) = lambda;

	Eigen::MatrixXd A = Eigen::MatrixXd::Ones(n, n);
	A << 2, 1, 1, 2;

	eigNewton(A, .0001, 10, z);

	std::cout << "eigenvector:" << std::endl;
	std::cout << z.head(n) << std::endl;
	std::cout << "eigenvalue: " << z(n) << std::endl;

	return 0;
}
