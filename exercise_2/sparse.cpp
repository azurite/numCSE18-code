#include <chrono>
#include <functional>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

std::vector<Eigen::Triplet<double>> MakeTripletList(int n) {
	int nnz = 3 * n - 2;
	std::vector<Eigen::Triplet<double>> tripletList(nnz);

	for(int i = 0; i < n; i++) {
		tripletList.push_back(Eigen::Triplet<double>(i, i, 2));

		if(i > 0) {
			tripletList.push_back(Eigen::Triplet<double>(i, i - 1, -1));
		}

		if(i < n - 1) {
			tripletList.push_back(Eigen::Triplet<double>(i, i + 1, -1));
		}
	}

	return tripletList;
}

double Runtime(const std::function<void(void)> &f) {
	double runtime = std::numeric_limits<double>::max();
	std::chrono::duration<double, std::ratio<1> > duration;
	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

	for(int i = 0; i < 10; i++) {
		start = std::chrono::high_resolution_clock::now();
		f();
		end = std::chrono::high_resolution_clock::now();

		duration = end - start;
		runtime = std::min(runtime, duration.count());
	}

	return runtime;
}

template <class T>
std::ostream & operator<< (std::ostream &os, const std::vector<T> &v) {
	os << "[";
	if (!v.empty()) {
		os << v[0];
		for (int i = 1; i < v.size(); ++i) os << ", " << v[i];
	}
    os << "]";

    return os;
}

int main() {
	// print small example of the tridiagonal matrix
	int m = 4;
	std::vector<Eigen::Triplet<double>> tripletList = MakeTripletList(m);
	Eigen::SparseMatrix<double> S_(m, m);
	S_.setFromTriplets(tripletList.begin(), tripletList.end());
	std::cout << "If n = " << m << ", then T equals" << std::endl;
	std::cout << Eigen::MatrixXd(S_) << std::endl;

	// matrix sizes for benchmark
	std::vector<int> N = {64, 128, 256, 512};
	std::cout << "LU decomposition of T, where n = " << N << std::endl;

	// set up variables for runtime measurement
	std::vector<double> runtimeSparse;
	std::vector<double> runtimeDense;

	for (int n : N) {
		tripletList = MakeTripletList(n);

		// sparse LU decomposition
		Eigen::SparseMatrix<double> S(n, n);
		S.setFromTriplets(tripletList.begin(), tripletList.end());

		runtimeSparse.push_back(Runtime([&]() {
			Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
			solver.analyzePattern(S);
			solver.factorize(S);
		}));

		// dense LU decomposition
		Eigen::MatrixXd D(S);

		runtimeDense.push_back(Runtime([&]() {
			D.fullPivLu();
		}));
	}

	std::cout << "Runtime in seconds using storage format..." << std::endl;
	std::cout << "...sparse: " << runtimeSparse << std::endl;
	std::cout << "...dense:  " << runtimeDense << std::endl;

	return 0;
}
