#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

VectorXd lowBidiagSolve(const MatrixXd &A, const VectorXd &y) {
	int n = y.size();
	VectorXd x(n);

	x(0) = y(0) / A(0, 0);

	for(int i = 1; i < n; i++) {
		x(i) = (y(i) - x(i - 1) * A(i, i - 1)) / A(i, i);
	}

	return x;
}


VectorXd upBidiagSolve(const MatrixXd &A, const VectorXd &y) {
	int n = y.size();
	VectorXd x(n);

	x(n - 1) = y(n - 1) / A(n - 1, n - 1);

	for(int i = n - 2; i >= 0; i--) {
		x(i) = (y(i) - x(i + 1) * A(i, i + 1)) / A(i, i);
	}

	return x;
}


VectorXd tridiagSolve(const MatrixXd &A, const VectorXd &y) {
	int n = y.size();
	VectorXd x(n);

	VectorXd b(n);

	// b_0 = d_0
	b(0) = A(0, 0);

	// b_i = d_i - (l_{i-1} * u_{i-1} / b_{i-1})
	// i = {1, ..., n - 1}
	for(int i = 1; i < n; i++) {
		b(i) = A(i, i) - A(i, i - 1) * A(i - 1, i) / b(i - 1);
	}

	// a_i = l_i / b_i 
	// i = {0, .., n - 2}
	VectorXd a = A.diagonal(-1).cwiseQuotient(b.head(n - 1));

	// c_i = u_i
	// i = {0, ..., n - 2}
	VectorXd c = A.diagonal(1);

	// Ax = y <==> LUx = y
	// Le = y (1)
	// Ux = e (2)


	// (1) the a's are the lower diagonal of L
	VectorXd e(n);
	e(0) = y(0);

	for(int i = 1; i < n; i++) {
		e(i) = (y(i) - e(i - 1) * a(i - 1));
	}

	// (2) the c's are the upper diagonal of U and the b's are the diagonal
	x(n - 1) = e(n - 1) / b(n - 1);

	for(int i = n - 2; i >= 0; i--) {
		x(i) = (e(i) - x(i + 1) * c(i)) / b(i);
	}

	return x;
}


/***** TESTING ******/

int main() {
	const double TOL = 10e-10;

	{ // lowBidiagSolve test1
	MatrixXd A(1,1);
	A << 4;
	VectorXd Y(1);
	Y << 7;

	if ((lowBidiagSolve(A,Y) - A.fullPivLu().solve(Y)).norm() < TOL) {
		std::cout << "lowBidiagSolve test1 successfull. ";
	} else {
		std::cout << "lowBidiagSolve test1 NOT successfull. ";
	}
	}

	{ // lowBidiagSolve test2
	MatrixXd A(4,4);
	A << 4, 0, 0, 0,
		 2, 9.4, 0, 0,
		 0, 7.2, 8, 0,
		 0, 0, -1.1,-4;
	VectorXd Y(4);
	Y << 1, 3, -2, 2;

	if ((lowBidiagSolve(A,Y) - A.fullPivLu().solve(Y)).norm() < TOL) {
		std::cout << "lowBidiagSolve test2 successfull. ";
	} else {
		std::cout << "lowBidiagSolve test2 NOT successfull. ";
	}
	}

	{ // upBidiagSolve test1
	MatrixXd A(1,1);
	A << 4;
	VectorXd Y(1);
	Y << 7;

	if ((upBidiagSolve(A,Y) - A.fullPivLu().solve(Y)).norm() < TOL) {
		std::cout << "upBidiagSolve test1 successfull. ";
	} else {
		std::cout << "upBidiagSolve test1 NOT successfull. ";
	}
	}

	{ // upBidiagSolve test2
	MatrixXd A(4,4);
	A << 4, 2.7, 0, 0,
		 0, 6.4, 3, 0,
		 0, 0, 8, 6,
		 0, 0, 0,-7;
	VectorXd Y(4);
	Y << 1, 3, -2, 2;

	if ((upBidiagSolve(A,Y) - A.fullPivLu().solve(Y)).norm() < TOL) {
		std::cout << "upBidiagSolve test2 successfull. ";
	} else {
		std::cout << "upBidiagSolve test2 NOT successfull. ";
	}
	}

	{ // tridiagSolve test1
	MatrixXd A(1,1);
	A << 4;
	VectorXd Y(1);
	Y << 7;

	if ((tridiagSolve(A,Y) - A.fullPivLu().solve(Y)).norm() < TOL) {
		std::cout << "tridiagSolve test1 successfull. ";
	} else {
		std::cout << "tridiagSolve test1 NOT successfull. ";
	}
	}

	{ // tridiagSolve test2
	MatrixXd A(4,4);
	A << 31.1, 2.7, 0, 0,
		 2, 22.4, 3, 0,
		 0, 11, 15, 6,
		 0, 0, -7.1,-21;
	VectorXd Y(4);
	Y << 1, 3, -2, 2;

	if ((tridiagSolve(A,Y) - A.fullPivLu().solve(Y)).norm() < TOL) {
		std::cout << "tridiagSolve test2 successfull. ";
	} else {
		std::cout << "tridiagSolve test2 NOT successfull. ";
	}
	}
}
