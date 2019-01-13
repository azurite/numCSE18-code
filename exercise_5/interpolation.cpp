#include <iostream>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

const std::complex<double> I(0,1); // imaginary unit
const double PI = 3.14159265359;

struct Newton {
	Newton(const Eigen::VectorXd &x) : _x(x), _a(x.size()) { }
	void Interpolate(const Eigen::VectorXd &y);
	double operator()(double x) const;

private:
	Eigen::VectorXd _x;	// nodes
	Eigen::VectorXd _a;	// coefficients
};

// Compute the coefficients in the Newton basis.
void Newton::Interpolate(const Eigen::VectorXd &y) {
	int size = _x.size();

	Eigen::MatrixXd M = Eigen::MatrixXd::Zero(size, size);
	M.col(0) = Eigen::VectorXd::Ones(size);

	for(int col = 1; col < size; col++) {
		for(int row = col; row < size; row++) {
			M(row, col) = M(row, col - 1) * (_x(row) - _x(col - 1));
		}
	}

	_a = M.triangularView<Eigen::Lower>().solve(y);
}

// Evaluate the interpolant at x.
double Newton::operator()(double x) const {
	double p = 0;

	for(int i = _x.size() - 1; i >= 0; i--) {
		p = p * (x - _x(i)) + _a(i);
	}

	return p;
}

struct Lagrange {
	Lagrange(const Eigen::VectorXd &x);
	void Interpolate(const Eigen::VectorXd &y) { _y = y; }
	double operator()(double x) const;

private:
	Eigen::VectorXd _x;	// nodes
	Eigen::VectorXd _l;	// weights
	Eigen::VectorXd _y;	// coefficients
};

// Compute the weights l for given nodes x.
Lagrange::Lagrange(const Eigen::VectorXd &x) : _x(x), _l(x.size()), _y(x.size()) {
	for(int i = 0; i < x.size(); i++) {
		double lambda_i = 1;

		for(int j = 0; j < x.size(); j++) {
			if(i != j) {
				lambda_i *= (x(i) - x(j));
			}
		}

		_l(i) = 1 / lambda_i;
	}
}

// Evaluate the interpolant at x.
double Lagrange::operator()(double x) const {
	double p = 0;
	double omega = 1;

	for(int i = 0; i < _x.size(); i++) {
		omega *= (x - _x(i));
	}

	for(int i = 0; i < _x.size(); i++) {
		p += _y(i) * omega * _l(i) / (x - _x(i));
	}

	return p;
}

struct Chebychev {
	Chebychev(const Eigen::VectorXd &x) : _x(x), _a(x.size()) {};
	void Interpolate(const Eigen::VectorXd &y);
	double operator()(double x) const;

private:
	Eigen::VectorXd _x; // nodes
	Eigen::VectorXd _a; // coefficients
};

void Chebychev::Interpolate(const Eigen::VectorXd &y) {
	int n = y.size() - 1;
	Eigen::VectorXcd b(2*n + 2);
	std::complex<double> constant = -PI * n / (double)(n + 1) * I;

	// build [b, b^R] of length 2(n + 1) = 2 * y.size()
	for(int k = 0; k <= n; k++) {
		b(k) = y(k) * std::exp(constant * (double)k);
		b(2 * n + 1 - k) = y(k) * std::exp(constant * (double)(2 * n + 1 - k));
	}

	Eigen::FFT<double> fft;
	Eigen::VectorXcd c = fft.inv(b);
	Eigen::VectorXd beta(c.size());

	for(int k = 0; k < c.size(); k++) {
		beta(k) = (c(k) * std::exp(PI * (k - n) / (2 * (n + 1)) * I)).real();
	}

	_a(0) = beta(n);

	for(int j = 1; j <= n; j++) {
		_a(j) = 2 * beta(j + n);
	}
}

double Chebychev::operator()(double x) const {
	int n = _x.size() - 1;
	Eigen::VectorXd b(n + 3);
	b(n + 2) = 0;
	b(n + 1) = 0;

	for(int k = n; k >= 0; k--) {
		double lead_const = (k == 0 ? 2 : 1);
		b(k) = lead_const * _a(k) + 2 * x * b(k + 1) - b(k + 2);
	}

	return 0.5 * (b(0) - b(2));
}

// Runge function
Eigen::VectorXd r(const Eigen::VectorXd &x) {
	return (1.0 / (1.0 + 25.0 * x.array() * x.array())).matrix();
}

int main() {
	int n = 5;
	Eigen::VectorXd x;
	x.setLinSpaced(5, -1.0, 1.0);
	Eigen::VectorXd y = r(x);

	Newton p(x);
	p.Interpolate(y); // correct result: p._a = [0.0384615, 0.198939, 1.5252, -3.31565, 3.31565]

	Lagrange q(x);    // correct result: p._l = [0.666667, -2.66667, 4, -2.66667, 0.666667]
	q.Interpolate(y);

	// Use Chebychev nodes instead of linearly spaced nodes
	for(int i = 0; i < x.size(); i++) {
		x(i) = std::cos((2 * i + 1) * PI / (2 * (x.size() + 1)));
	}
	y = r(x);

	Chebychev c(x);
	c.Interpolate(y);

	// Compute difference of p and q.
	int m = 22;
	double offset = 0.08333333333;
	x.setLinSpaced(m, -1.0 + offset, 1.0 - offset);
	double norm2 = .0;
	for (int i = 0; i < m; ++i) {
		double d = p(x(i)) - q(x(i));
		norm2 += d * d;
	}

	// By uniquenss of the interpolation polynomial, we expect p = q.
	std::cout << "This number should be close to zero: " << norm2 << std::endl;

	norm2 = 0;
	for(int i = 0; i < m; ++i) {
		double d = c(x(i)) - p(x(i));
		norm2 += d * d;
	}

	// By uniquenss of the interpolation polynomial, we expect c = p.
	std::cout << "This number should be close to zero: " << norm2 << std::endl;

	norm2 = 0;
	for(int i = 0; i < m; ++i) {
		double d = c(x(i)) - q(x(i));
		norm2 += d * d;
	}

	// By uniquenss of the interpolation polynomial, we expect c = q.
	std::cout << "This number should be close to zero: " << norm2 << std::endl;

	return 0;
}
