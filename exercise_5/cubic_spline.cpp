// to compile run: g++ -std=gnu++11 cubic_spline.cpp -lmgl
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <mgl2/mgl.h>

using namespace Eigen;

MatrixXd cubicSpline(const VectorXd &T, const VectorXd &Y) {
	// returns the matrix representing the spline interpolating the data
	// with abscissae T and ordinatae Y. Each column represents the coefficients
	// of the cubic polynomial on a subinterval.
	// Assumes T is sorted, has no repeated elements and T.size() == Y.size().

	int n = T.size() - 1; // T and Y have length n+1
	MatrixXd spline(4, n);
	VectorXd h(n + 1);
	VectorXd r(n - 1);

	// h(0) is never used by convention
	for(int i = 1; i <=n; i++) {
		h(i) = T(i) - T(i - 1);
	}

	for(int i = 1; i < n; i++) {
		r(i - 1) = (Y(i + 1) - Y(i)) / h(i + 1) - (Y(i) - Y(i - 1)) / h(i);
	}

	MatrixXd A = MatrixXd::Zero(n - 1, n - 1);

	for(int i = 1; i <= n - 2; i++) {
		A(i - 1, i - 1) = (h(i) + h(i + 1)) / 3;
		A(i - 1, i) = h(i + 1) / 6;
		A(i, i - 1) = h(i + 1) / 6;
	}

	A(n - 2, n - 2) = (h(n - 1) + h(n)) / 3;
	VectorXd rho = A.fullPivLu().solve(r);

	// eigen rho(i - 1) contains formula rho(i)
	for(int i = 1; i <= n; i++) {
		double rho_im1 = i == 1 ? 0 : rho(i - 2);
		double rho_i = i == n ? 0 : rho(i - 1);

		spline(0, i - 1) = Y(i - 1);
		spline(1, i - 1) = (Y(i) - Y(i - 1)) / h(i) - (h(i) * (2 * rho_im1 + rho_i)) / 6;
		spline(2, i - 1) = rho_im1 / 2;
		spline(3, i - 1) = (rho_i - rho_im1) / (6 * h(i));
	}

	return spline;
}

VectorXd evalCubicSpline(const MatrixXd &S, const VectorXd &T, const VectorXd &evalT) {
	// Returns the values of the spline S calculated in the points X.
	// Assumes T is sorted, with no repetetions.

	VectorXd out(evalT.size());

	for(int i = 0; i < evalT.size(); i++) {
		double t = evalT(i);

		for(int j = 1; j < T.size(); j++) {
				if((j == 1 && (t >= T(j - 1) && t <= T(j))) || (t > T(j - 1) && t <= T(j))) {
					double delta_t = t - T(j - 1);
					out(i) = ((S(3, j - 1) * delta_t + S(2, j - 1)) * delta_t + S(1, j - 1)) * delta_t + S(0, j - 1);
				}
		}
	}

	return out;
}

int main() {
	// tests
	VectorXd T(9);
	VectorXd Y(9);
	T << 0, 0.4802, 0.7634, 1, 1.232, 1.407, 1.585, 1.879, 2;
	Y << 0., 0.338, 0.7456, 0, -1.234, 0 , 1.62, -2.123, 0;

	int len = 1 << 9;
	VectorXd evalT = VectorXd::LinSpaced(len, T(0), T(T.size()-1));

	VectorXd evalSpline = evalCubicSpline(cubicSpline(T, Y), T, evalT);

 	mglData datx, daty;
	datx.Link(evalT.data(), len);
	daty.Link(evalSpline.data(), len);
	mglGraph gr;
	gr.SetRanges(0, 2, -3, 3);
	gr.Plot(datx, daty, "0");
	gr.WriteFrame("spline.eps");
}
