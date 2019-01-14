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
	// the lecture notes model the data output points as [y_0, y_1, ..., y_n]^T
	int n = y.size() - 1;
	Eigen::VectorXcd b(2*n + 2);

	for(int k = 0; k < b.size(); k++) {
		// z_k from lecture notes (5.40)
		double z_k = k > n ? y(2 * n + 1 - k) : y(k);

		// b(k) from lecture notes (5.41)
		b(k) = z_k * std::exp(-PI * n * k / (n + 1) * I);
	}

	Eigen::FFT<double> fft;
	Eigen::VectorXcd c = fft.inv(b);
	Eigen::VectorXd beta(c.size());

	// extract beta's according to lecture notes (5.41) and index substitution k = j + n
	// to only deal with non-negative indices
	for(int k = 0; k < c.size(); k++) {
		// possible error in lecture notes should be j = -n, ..., n + 1 instead of j = 0, ..., 2n + 1
		beta(k) = (c(k) * std::exp(PI * (k - n) / (2 * (n + 1)) * I)).real();
	}

	// from lecture notes (5.38)
	// we have to acces beta(j + n) because of the substitution we did above
	// at the end we are only interested in alpha_j's from j = 0,...,n like in (5.32)
	for(int j = 0; j <= n; j++) {
		if(j == 0) {
			_a(j) = beta(j + n);
		}
		else {
			_a(j) = 2 * beta(j + n);
		}
	}
}

// Evaluate the interpolant at x.
double Chebychev::operator()(double x) const {
	// Clenshaw algorithm from lecture notes (5.35) and (5.36)
	int n = _a.size() - 1;
	Eigen::VectorXd beta(n + 3);
	beta(n + 2) = 0;
	beta(n + 1) = 0;

	for(int k = n; k >= 0; k--) {
		beta(k) = (k == 0 ? 2 : 1) * _a(k) + 2 * x * beta(k + 1) - beta(k + 2);
	}

	return 0.5 * (beta(0) - beta(2));
}

// Runge function
Eigen::VectorXd r(const Eigen::VectorXd &x) {
	return (1.0 / (1.0 + 25.0 * x.array() * x.array())).matrix();
}

int main() {
	int n = 5;

	/*
	Eigen::VectorXd x;
	x.setLinSpaced(5, -1.0, 1.0);
	Eigen::VectorXd y = r(x);

	Newton p(x);
	p.Interpolate(y); // correct result: p._a = [0.0384615, 0.198939, 1.5252, -3.31565, 3.31565]

	Lagrange q(x);    // correct result: p._l = [0.666667, -2.66667, 4, -2.66667, 0.666667]
	q.Interpolate(y);
	*/

	// Use Chebychev nodes instead of linearly spaced nodes
	Eigen::VectorXd x(n + 1);

	for(int i = 0; i < x.size(); i++) {
		x(i) = std::cos((2 * i + 1) * PI / (2 * (n + 1)));
	}

	Eigen::VectorXd y = r(x);

	Newton p(x);
	p.Interpolate(y);

	Lagrange q(x);
	q.Interpolate(y);

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
