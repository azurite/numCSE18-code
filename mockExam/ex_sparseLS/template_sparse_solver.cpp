#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using Triplet = Eigen::Triplet<double>;
using Triplets = std::vector<Triplet>;

using Vector = Eigen::VectorXd;
using Matrix = Eigen::SparseMatrix<double>;

//! \brief Efficiently construct the sparse matrix A given c, i_0 and j_0
//! \param[in] c contains entries c_i for matrix A
//! \param[in] i0 row index i_0
//! \param[in] j0 column index j_0
//! \return Sparse matrix A
Matrix buildA(const Vector & c, unsigned int i0, unsigned int j0) {
    assert(i0 > j0);
    
    unsigned int n = c.size() + 1;
    Matrix A(n,n);
    Triplets triplets;
    
    // TODO: problem 2a, construct and return the matrix A given c, i0 and j0
    
    return A;
}

//! \brief Solve the system Ax = b with optimal complexity O(n)
//! \param[in] c contains entries c_i for matrix A
//! \param[in] b r.h.s. vector
//! \param[in] i0 row index
//! \param[in] j0 column index
//! \return Solution x, s.t. Ax = b
Vector solveLSE(const Vector & c, const Vector & b, unsigned int i0, unsigned int j0) {
    assert(c.size() == b.size()-1 && "Size mismatch!");
    assert(i0 > j0);
    
    // Allocate solution vector
    Vector ret(b.size());
    
    // TODO: problem 2b, solve system Ax = b in O(n)
    
    return ret;
}

int main(int, char**) {
    // Setup data for problem
    unsigned int n = 15; // A is n x n matrix, b has length x
    
    unsigned int i0 = 6, j0 = 4;
    
    Vector b = Vector::Random(n); // Random vector for b
    Vector c = Vector::Random(n-1); // Random vector for c
    
    //// PROBLEM 2a
    std::cout << "*** PROBLEM 2a:" << std::endl;
    
    // Solve sparse system using sparse LU and our own routine
    Matrix A = buildA(c, i0, j0);
    A.makeCompressed();
    Eigen::SparseLU<Matrix> splu;
    splu.analyzePattern(A); 
    splu.factorize(A);
    
    std::cout << "Error: " << std::endl << ( splu.solve(b) - solveLSE(c, b, i0, j0) ).norm() << std::endl;
}
