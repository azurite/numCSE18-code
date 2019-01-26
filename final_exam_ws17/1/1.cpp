#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

VectorXd lowBidiagSolve(const MatrixXd &A, const VectorXd &y) {
	int n = y.size();
	VectorXd x(n);

	// TODO: implement lower bidiagonal matrix solver

	return x;
}


VectorXd upBidiagSolve(const MatrixXd &A, const VectorXd &y) {
	int n = y.size();
	VectorXd x(n);

	// TODO: implement upper bidiagonal matrix solver

	return x;
}


VectorXd tridiagSolve(const MatrixXd &A, const VectorXd &y) {
	int n = y.size();
	VectorXd x(n);

	// TODO: implement tridiagonal matrix solver

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
