#include <iostream>
#include <limits>
#include <eigen3/Eigen/Dense>
#include <Eigen/SVD>

using namespace Eigen;


// routine: llsq_gsol
// " Computes the Generalized linear least squares solution "
// (in)  A: model matrix
// (in)  b: rhs vector
// (out)  x: solution vector
double llsq_gsol(const MatrixXd& A, const VectorXd& b, VectorXd& x) {
	
	double llsq_err = 0;
	
	JacobiSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV);

	const int m = A.rows();
	const int n = A.cols();
	int r = 0;

	VectorXd sv = svd.singularValues();

	// calculate rank of A
	for(int i = 0; i < sv.size(); i++) {
		if(std::abs(sv(i)) > 1e-12) {
			r++;
		}
	}

	MatrixXd U = svd.matrixU();
	MatrixXd V = svd.matrixV();

	MatrixXd U1 = U.block(0, 0, m, r);
	MatrixXd U2 = U.block(0, r, m, m - r);
	MatrixXd V1 = V.block(0, 0, n, r);

	// singular values are sorted in decreasing order
	MatrixXd sigma_inv = sv.head(r).cwiseInverse().asDiagonal();
	
	x = V1 * sigma_inv * U1.transpose() * b;
	llsq_err = (U2.transpose() * b).norm();

	return llsq_err;
}


// routine: run_llsq_gsol
// " Executes the Generalized LLSQ solver "
void run_llsq_gsol() {

	int m = 5;
	int n = 4;

	VectorXd k(m);
	VectorXd h_k(m);
	MatrixXd A(m, n);

	k << 11, 32, 77, 121, 152;
	h_k << 9, 10, 12, 14, 15;

	A.col(0) = VectorXd::Ones(m);
	
	for(int i = 0; i < m; i++) {
		A(i, 1) = std::sin(2 * M_PI * k(i) / 366);
		A(i, 2) = std::cos(2 * M_PI * k(i) / 366);
		A(i, 3) = std::sin(M_PI * k(i) / 366) * std::cos(M_PI * k(i) / 366);
	}

	VectorXd x;
	double err = llsq_gsol(A, h_k, x);

	std::cout << "solution: " << x.transpose() << std::endl;
	std::cout << "error: " << err << std::endl;
}


/***** TESTING ******/
// "" Do NOT CHANGE the routines below ""

void runTests() {
    
    bool success = true;
    double TOL = 1e-5;
    
    size_t m = 5;
	size_t n = 4;
    
    VectorXd t_test(m), b_test(m);
	t_test << 11, 32, 77, 121, 152;
	b_test << 10, 11, 12, 15, 16;
    
    MatrixXd A_test(m,n);
    A_test << 1, 0.187718565437199, 0.982222856682841, 0.0938592827185993,\
              1, 0.522132581076825, 0.85286433140216,  0.261066290538412,\
              1, 0.969178065439254, 0.246361274293315, 0.484589032719627,\
              1, 0.87448095783104, -0.485059846195196, 0.43724047891552,\
              1, 0.507415093293846, -0.861701759948068, 0.253707546646923;
    
    std::cout << "A_test =" << std::endl << A_test << "\n" << std::endl;
    
    VectorXd x_test(n), x_ans(n);
    double llsq_err_test, llsq_err_ans=0.67829;
    x_ans << 13.455, -0.238115, -3.21762, -0.119057;
    llsq_err_test = llsq_gsol(A_test, b_test, x_test);
    if((x_ans-x_test).norm() > TOL * x_ans.norm()){
		std::cout << "\nTest llsq_gsol FAILED.\n";
		std::cerr << "Your answer:\n" << x_test << "\n\n" << "Correct answer:\n" << x_ans << "\n\n";
		success = false;
	}
	
	if (success) {
	    std::cout << "\n All tests PASSED.\n" << std::endl;
	}

}


int main(int argc, char** argv) {
    
    int runTests_flag=0;
    
    if (argc==2) {
        runTests_flag = std::stoi(argv[1]);
    }
    
    if (runTests_flag) {
        std::cout << "\nRun tests ...." << std::endl;
        runTests();
    }
    
    std::cout << "\nRun generalized linear least squares solver ...." << std::endl;
    run_llsq_gsol();
    
}
