#include <iostream>
#include <vector>

#include <Eigen/Dense>

double f1(const Eigen::Vector2d &x) {
	return x(0) * x(0) + 2. * x(1) * x(1) - 1.;
}

double f2(const Eigen::Vector2d &x) {
	return x(1) - x(0) * x(0);
}

Eigen::Vector2d f(const Eigen::Vector2d &x) {
	return Eigen::Vector2d(f1(x), f2(x));
}

Eigen::Matrix2d Jacobian(const Eigen::Vector2d &x) {
	Eigen::Matrix2d J;
	J <<  2*x(0), 4*x(1),
	     -2*x(0), 	  1;
	return J;
}

Eigen::Vector2d Newton(Eigen::Vector2d x, int n) {
	for(int i = 0; i < n; i++) {
		x += Jacobian(x).fullPivLu().solve(-f(x));
	}	

	return x;
}

Eigen::IOFormat ShortFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[", "]");

int main() {
	int n = 100;
	std::vector<Eigen::Vector2d> startingPoints(3);

	startingPoints[0] = Eigen::Vector2d(-1., 1.);
	startingPoints[1] = Eigen::Vector2d(1., 1.);
	startingPoints[2] = Eigen::Vector2d(-2., -2.);

	std::cout << "---------------------------" << std::endl;
	for (const Eigen::Vector2d &x : startingPoints) {
		Eigen::Vector2d y = Newton(x, n);
		std::cout << "   x = " << x.format(ShortFmt) << std::endl;
		std::cout << "   y = " << y.format(ShortFmt) << std::endl;
		std::cout << "f(y) = " << f(y).format(ShortFmt) << std::endl;
		std::cout << "---------------------------" << std::endl;
	}

	return 0;
}
