#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

//!
//! \brief buildDistanceLSQMatrix Efficiently build the system matrix
//! mapping positions to distances.
//! \param n Number of points (including $x_1$).
//! \return The system matrix $A$.
//!

SparseMatrix<double> buildDistanceLSQMatrix(int n) {
	SparseMatrix<double> A(n*(n-1)/2, n-1);

	std::vector<Triplet<double>> triplets;

	int row = 0;
	for(int i = 0; i < n - 1; i++) {
		for(int j = i; j < n - 1; j++) {
			triplets.push_back(Triplet<double>(row, j, 1));
			if(i > 0) {
				triplets.push_back(Triplet<double>(row, i - 1, -1));			
			}
			row++;
		}
	}
	
	A.setFromTriplets(triplets.begin(), triplets.end());
	A.makeCompressed();
	return A;
}


//!
//! \brief estimatePointsPositions Return positions (without $x_1$).
//! The coordinate $x_1$ is assumed $x_1 = 0$.
//! \param D An $n \times n$ anti-symmetric matrix of distances.
//! \return Vector of positions $x_2, \dots, x_n$.
//!

VectorXd estimatePointsPositions(const MatrixXd& D) {

	VectorXd x;
	int n = D.rows();
	VectorXd b(n * (n - 1) / 2);

	// get matrix A
	SparseMatrix<double> A = buildDistanceLSQMatrix(n);

	// transform matrix D into vector with distance entries
	int idx = 0;	
	for(int i = 0; i < n; i++) {
		for(int j = i + 1; j < n; j++) {
			b(idx++) = -D(i, j);		
		}
	}
	
	// Solve sparse LSQ with normal equations
	SparseLU<SparseMatrix<double>> lu(A.transpose() * A);
	x = lu.solve(A.transpose() * b);	

	return x;
}


int main() {
    
    int n = 5;
    
    // PART 1: build and print system matrix A
    std::cout << "**** PART 1 ****" << std::endl;
    std::cout << "The Matrix A is:"
              << std::endl
              << buildDistanceLSQMatrix(n)
              << std::endl;

    // PART 2: solve the LSQ system and find positions
    std::cout << "**** PART 2 ****" << std::endl;
    // Vector of positions
    n = 5;

    // Build D
    MatrixXd D(n,n);
    D << 0,   -2.1, -3,   -4.2, -5,
         2.1,  0,   -0.9, -2.2, -3.3,
         3,    0.9,  0,   -1.3, -1.1,
         4.2,  2.2,  1.3,  0,   -1.1,
         5,    3.3,  1.1,  1.1,  0;
    std::cout << "The matrix D is:"
              << std::endl
              << D
              << std::endl;

    // Find out positions
    VectorXd x_recovered = estimatePointsPositions(D);
    std::cout << "The positions [x_2, ..., x_n] obtained from the LSQ system are:"
              << std::endl
              << x_recovered
              << std::endl;
}
