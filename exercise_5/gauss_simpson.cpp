#include <cmath>
#include <iostream>
#include <functional>
#include <vector>
#include <algorithm>

double GaussLegendre5(const std::function<double(double)> &f, double a, double b) {
	std::vector<double> c = { -0.90617984593, -0.53846931010, 0.0, 0.53846931010, 0.90617984593 };
	std::vector<double> w = { 0.23692688505, 0.47862867049, 0.56888888888, 0.47862867049, 0.23692688505 };

	std::transform(c.begin(), c.end(), c.begin(), [&](double cj) {
		return 0.5 * ((1 - cj) * a + (1 + cj) * b);
	});

	std::transform(w.begin(), w.end(), w.begin(), [&](double wj) {
		return 0.5 * (b - a) * wj;
	});

	double integral = 0;

	for(int i = 0; i < c.size(); i++) {
		integral += f(c[i]) * w[i];
	}

	return integral;
}

double CompositeSimpson(const std::function<double(double)> &f, const std::vector<double> &x) {
	int m = x.size() - 1;
	double integral = 0;
	double a = x[0];
	double b = x[m];

	double sixth = 1 / (double)6;
	double third = 2 / (double)3;

	integral += sixth * (x[1] - x[0]) * f(a) + sixth * (x[m] - x[m - 1]) * f(b);

	for(int j = 1; j < m; j++) {
		integral += sixth * (x[j + 1] - x[j - 1]) * f(x[j]);
	}

	for(int j = 1; j < m + 1; j++) {
		integral += third * (x[j] - x[j - 1]) * f(0.5 * (x[j] + x[j - 1]));
	}

	return integral;
}

std::vector<double> LinSpace(int n, double a, double b) {
	std::vector<double> x(n);
	double d = (b - a) / (n - 1);
	for (int i = 0; i < n; ++i) {
		x[i] = i * d + a;
	}
	return x;
}

double f(double x) {
	return std::sqrt(x);
}

double F(double x) {
	double y = std::sqrt(x);
	return 2.0 / 3.0 * y * y * y;
}

int main() {
	int n = 5;
	double a = 0.0;
	double b = 1.0;

	std::cout << "Gauss-Legendre: " << GaussLegendre5(f, a, b) << std::endl;
	std::cout << "Simpson: " << CompositeSimpson(f, LinSpace(n, a, b)) << std::endl;
	std::cout << "Exact value: " << F(b) - F(a) << std::endl;

	return 0;
}
