////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

/* \brief Performs Gram-Schidt orthonormalization
 * Given a matrix $\mathbf{A}$ of linearly independent columns,
 * returns the result of a Gram-Schmidt orthonormalization.
 * Unstable GS algorithm: output is prone to cancellation issues.
 * \param[in] $\mathbf{A}$ Matrix of linearly independent columns
 * \return Matrix with ONB of $span(a_1, \cdots, a_n)$ as columns
 */

MatrixXd gram_schmidt(const MatrixXd &A) {
	// We create a matrix Q with the same size and data of A
	MatrixXd Q(A);

	for(int k = 0; k < A.cols(); k++) {
		VectorXd sum = VectorXd::Zero(A.cols());

		for(int j = 0; j < k; j++) {
			sum += (Q.col(j).dot(A.col(k))) * Q.col(j);
		}

		Q.col(k) = A.col(k) - sum;
		Q.col(k).normalize();
	}

	return Q;
}

int main(void) {

	// Orthonormality test
	unsigned int n = 4;
	MatrixXd A = MatrixXd::Random(n, n);
	MatrixXd Q = gram_schmidt(A);

	// "How far is Q from being orthonormal?"
	MatrixXd I =  Q.transpose() * Q ; 	// for orthogonal matrices this should be the identity
	std::cout << I << std::endl;

	// Error has to be small, but not zero (why?)
	double err = (I - MatrixXd::Identity(n, n)).norm();
	std::cout << "Error is: " << err << std::endl;

	return 0;
}
