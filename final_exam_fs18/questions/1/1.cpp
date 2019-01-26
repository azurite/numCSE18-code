#include <iostream>
#include <vector>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using namespace Eigen;

MatrixXd quadraticSpline(const VectorXd &x, const VectorXd &y) {
	int N = x.size()-1;
	MatrixXd out(N, 3);

	VectorXd delta_x = x.tail(N) - x.head(N);
	VectorXd delta_x2 = (delta_x.array() * delta_x.array()).matrix();

	MatrixXd A = MatrixXd::Zero(3*N, 3*N);

	A.block(0, 2*N, N, N) = MatrixXd::Identity(N, N);
	A.block(N, 2*N, N, N) = MatrixXd::Identity(N, N);
	A.block(N, 0, N, N) = delta_x2.asDiagonal();
	A.block(N, N, N, N) = delta_x.asDiagonal();
	A.block(2*N, 0, N, N) = (2 * delta_x).asDiagonal();

	MatrixXd M = MatrixXd::Identity(N, N);
	
	for(int i = 0; i < N; i++) {
		M(i, (i + 1) % N) = -1;
	}
	
	A.block(2*N, N, N, N) = M;

	VectorXd rhs = VectorXd::Zero(3*N);
	rhs.segment(0, N) = y.head(N);
	rhs.segment(N, N) = y.tail(N);

	VectorXd sol = A.fullPivLu().solve(rhs);

	out.col(0) = sol.segment(0, N);
	out.col(1) = sol.segment(N, N);
	out.col(2) = sol.segment(2*N, N);

	return out;
}

MatrixXd quadraticSplineFast(const VectorXd &x, const VectorXd &y) {
	int N = x.size()-1;
	MatrixXd out(N, 3);

	VectorXd delta_x = x.tail(N) - x.head(N);
	VectorXd delta_x2 = (delta_x.array() * delta_x.array()).matrix();

	std::vector<Triplet<double>> triplets(7*N);

	for(int i = 0; i < N; i++) {
		triplets.push_back(Triplet<double>(i, 2*N + i, 1));
		triplets.push_back(Triplet<double>(N + i, 2*N + i, 1));
		triplets.push_back(Triplet<double>(2*N + i, N + i, 1));
		triplets.push_back(Triplet<double>(2*N + i, N + ((i + 1) % N), -1));
		triplets.push_back(Triplet<double>(N + i, i, delta_x2(i)));
		triplets.push_back(Triplet<double>(N + i, N + i, delta_x(i)));
		triplets.push_back(Triplet<double>(2*N + i, i, 2 * delta_x(i)));
	}

	SparseMatrix<double> A(3*N, 3*N);
	A.setFromTriplets(triplets.begin(), triplets.end());
	A.makeCompressed();

	SparseLU<SparseMatrix<double>> splu;
	splu.compute(A);

	VectorXd rhs = VectorXd::Zero(3*N);
	rhs.segment(0, N) = y.head(N);
	rhs.segment(N, N) = y.tail(N);

	VectorXd sol = splu.solve(rhs);
	out.col(0) = sol.segment(0, N);
	out.col(1) = sol.segment(N, N);
	out.col(2) = sol.segment(2*N, N);
	
	return out;
}




/***** TESTING ******/


int main() {
	const double TOL = 10e-5;
	int testnum;
	std::cin >> testnum;
	for(int i=1; i<testnum; i++){
		int N;
		std::cin >> N;
		VectorXd x(N+1), y(N+1);
		MatrixXd ans(N,3);

		for(int j=0; j<N+1; j++)
			std::cin >> x(j);
		for(int j=0; j<N+1; j++)
			std::cin >> y(j); 

		for(int j=0; j<N; j++){
			std::cin >> ans(j,0);
			std::cin >> ans(j,1);
			std::cin >> ans(j,2);
		}
		
		int flag;
		std::cin >> flag;

		auto start = std::chrono::steady_clock::now();	
		MatrixXd test;
		if(flag)
			test = quadraticSpline(x,y);
		else 
			test = quadraticSplineFast(x,y);
		auto end = std::chrono::steady_clock::now();
		double difftime = std::chrono::duration<double,std::milli>(end-start).count();

		if((test-ans).norm() > TOL * ans.norm()){
			std::cout << "Test " << i << " INCORRECT.\n";
			std::cerr << "Your answer:\n" << test << "\n\n" << "Correct answer:\n" << ans << "\n\n";
			break;
		}
		else
			std::cout << "Test " << i << " with N=" << N << " correct in " << difftime << "ms.\n";
	}
}
