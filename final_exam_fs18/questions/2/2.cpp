#include <iostream>
#include <limits>
#include <eigen3/Eigen/Dense>

using namespace Eigen;


// routine: llsq_gsol
// " Computes the Generalized linear least squares solution "
// (in)  A: model matrix
// (in)  b: rhs vector
// (out)  x: solution vector
double llsq_gsol(const MatrixXd& A, const VectorXd& b, VectorXd& x) {
	
	double llsq_err = 0;
	
	// TODO: Solve exercise 2.c
	
	return llsq_err;
}


// routine: run_llsq_gsol
// " Executes the Generalized LLSQ solver "
void run_llsq_gsol() {

    // TODO: Solve exercise 2.d
    
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
