#include <iostream>

#include <Eigen/Dense>
#include <Eigen/QR>

/* @brief QR decomposition from Cholesky decomposition
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
void CholeskyQR(const Eigen::MatrixXd &A, Eigen::MatrixXd &R, Eigen::MatrixXd &Q) {
	Eigen::MatrixXd L = (A.transpose() * A).llt().matrixL();
	Q = L.triangularView<Eigen::Lower>().solve(A.transpose()).transpose();
	R = L.transpose();
}

/* @brief Direct QR decomposition
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
void DirectQR(const Eigen::MatrixXd &A, Eigen::MatrixXd &R, Eigen::MatrixXd &Q) {

	int m = A.rows();
	int n = A.cols();

	Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
	Q = qr.householderQ() * Eigen::MatrixXd::Identity(m, n);
	R = (Q.transpose() * A);
}

int main() {
	int m = 3;
	int n = 2;

	Eigen::MatrixXd A(m,n);
	double epsilon = 1.e-8;
	//A << 3, 5, 1, 9, 7, 1;
	A << 1, 1, epsilon, 0, 0, epsilon;
	std::cout << "A =" << std::endl << A << std::endl;

	Eigen::MatrixXd R, Q;

	CholeskyQR(A, R, Q);
	std::cout << "CholeskyQR: ===========" << std::endl;
	std::cout << "R =" << std::endl << R << std::endl;
	std::cout << "Q =" << std::endl << Q << std::endl;

	DirectQR(A, R, Q);
	std::cout << "DirectQR: =============" << std::endl;
	std::cout << "R =" << std::endl << R << std::endl;
	std::cout << "Q =" << std::endl << Q << std::endl;

	return 0;
}
