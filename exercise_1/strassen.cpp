#include<bits/stdc++.h>
#include<eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

MatrixXf mult(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C = MatrixXf::Zero(N, N);

	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			for(int k = 0; k < N; k++) {
				C(i, j) += A(i, k) * B(k, j);
			}
		}
	}

	return C;
}

MatrixXf mult_rec(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C(N, N);

	if(N == 1) { return A * B; }
	int M = N / 2;

	MatrixXf A11 = A.topLeftCorner(M, M);
	MatrixXf A12 = A.topRightCorner(M, M);
	MatrixXf A21 = A.bottomLeftCorner(M, M);
	MatrixXf A22 = A.bottomRightCorner(M, M);

	MatrixXf B11 = B.topLeftCorner(M, M);
	MatrixXf B12 = B.topRightCorner(M, M);
	MatrixXf B21 = B.bottomLeftCorner(M, M);
	MatrixXf B22 = B.bottomRightCorner(M, M);

	C.topLeftCorner(M, M) = mult_rec(A11, B11) + mult_rec(A12, B21);
	C.topRightCorner(M, M) = mult_rec(A11, B12) + mult_rec(A12, B22);
	C.bottomLeftCorner(M, M) = mult_rec(A21, B11) + mult_rec(A22, B21);
	C.bottomRightCorner(M, M) = mult_rec(A21, B12) + mult_rec(A22, B22);

	return C;
}


MatrixXf strassen(const MatrixXf &A, const MatrixXf &B) {
	int N = A.rows();
	MatrixXf C(N, N);

	if(N == 1) { return A * B; }
	int M = N / 2;

	MatrixXf A11 = A.topLeftCorner(M, M);
	MatrixXf A12 = A.topRightCorner(M, M);
	MatrixXf A21 = A.bottomLeftCorner(M, M);
	MatrixXf A22 = A.bottomRightCorner(M, M);

	MatrixXf B11 = B.topLeftCorner(M, M);
	MatrixXf B12 = B.topRightCorner(M, M);
	MatrixXf B21 = B.bottomLeftCorner(M, M);
	MatrixXf B22 = B.bottomRightCorner(M, M);

	MatrixXf M1 = strassen(A11 + A22, B11 + B22);
	MatrixXf M2 = strassen(A21 + A22, B11);
	MatrixXf M3 = strassen(A11, B12 - B22);
	MatrixXf M4 = strassen(A22, B21 - B11);
	MatrixXf M5 = strassen(A11 + A12, B22);
	MatrixXf M6 = strassen(A21 - A11, B11 + B12);
	MatrixXf M7 = strassen(A12 - A22, B21 + B22);

	C.topLeftCorner(M, M) = M1 + M4 - M5 + M7;
	C.topRightCorner(M, M) = M3 + M5;
	C.bottomLeftCorner(M, M) = M2 + M4;
	C.bottomRightCorner(M, M) = M1 - M2 + M3 + M6;

	return C;
}

int main() {
	srand(time(0));
	cout << setprecision(6) << setfill(' ');

	for (int i = 1; i < 9; i++) {
		int N = 1 << i;
		cout << "Matrix size = " << N << endl;
		MatrixXd AA = MatrixXd::Random(N, N);
		MatrixXd BB = MatrixXd::Random(N, N);
		MatrixXd ans = AA*BB;
		MatrixXf A = AA.cast<float>();
		MatrixXf B = BB.cast<float>();

		auto start = std::chrono::steady_clock::now();
		MatrixXf W = mult(A, B);
		auto finish = std::chrono::steady_clock::now();
		cout << setw(24) << " " <<  setw(15)
			<< "Time (s)" << setw(20) << "Error (l2-norm)"  << endl;
		cout << setw(24) << "Naive iterative "<< setw(15)
			<< std::chrono::duration_cast<std::chrono::duration<double> >
			(finish - start).count() << setw(20) << (W.cast<double>() - ans).norm() << endl;

		start = std::chrono::steady_clock::now();
		MatrixXf X = mult_rec(A, B);
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Naive recursive "<< setw(15)
			<< std::chrono::duration_cast<std::chrono::duration<double> >
			(finish - start).count() << setw(20) << (X.cast<double>() - ans).norm() << endl;
		start = std::chrono::steady_clock::now();
		MatrixXf Y = strassen(A, B);
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Strassen recursive " << setw(15)
			<< std::chrono::duration_cast<std::chrono::duration<double> >
			(finish - start).count() << setw(20) << (Y.cast<double>() - ans).norm() << endl;
		start = std::chrono::steady_clock::now();
		MatrixXf Z = A*B;
		finish = std::chrono::steady_clock::now();
		cout << setw(24) << "Eigen built-in "<< setw(15)
			<< std::chrono::duration_cast<std::chrono::duration<double> >
			(finish - start).count() << setw(20) << (Z.cast<double>() - ans).norm() << "\n\n\n";
	}
}
