#include <iostream>
#include <Eigen/Dense>
#include "gaussQuad.hpp"


using namespace Eigen;

const double PI = 3.141592653589793;

// small helper function to evaluate f_alpha_(x)
double f(VectorXd alpha, double x) {
	return std::exp(alpha(0) + alpha(1) * x + alpha(2) * x * x);
}

// rho(x) = (1, x, x^2)^T
VectorXd rho(double x) {
	VectorXd r(3);
	r << 1, x, x * x;
	return r;
}

// helper function to compute the jacobian of F
MatrixXd m_rho(double x) {
	MatrixXd m(3, 3);
	double x2 = x * x;
	double x3 = x * x * x;
	double x4 = x * x * x * x;

	m << 1,  x,  x2,
	     x,  x2, x3,
	     x2, x3, x4;

	return m;
}

// maps quadrature weights and nodes to arbitrary intervals [a, b]
struct QuadRule map_to_interval(const QuadRule& qr, const int a, const int b) {
	struct QuadRule m_qr;
	m_qr.nodes = VectorXd::Zero(qr.nodes.size());
	m_qr.weights = VectorXd::Zero(qr.weights.size());

	// map nodes from [-1, 1] to [a, b]
	for(int i = 0; i < qr.nodes.size(); i++) {
		m_qr.nodes(i) = (0.5 * (1 - qr.nodes(i)) * a) + (0.5 * (1 + qr.nodes(i)) * b);
	}

	// scale the weights with |[a, b]| / |[-1, 1]|
	for(int i = 0; i < qr.weights.size(); i++) {
		m_qr.weights(i) = 0.5 * (b - a) * qr.weights(i);
	}

	return m_qr;
}

// routine: evalF
// " Computes F(\alpha) "
// (in)  qr: quad nodes and weights
// (in)  \alpha: model parameters
// (in)  u: vector of moments
// (out) F(\alpha)
VectorXd evalF(const QuadRule& qr, const VectorXd& alpha, const VectorXd& u) {

	int m = alpha.size();
	VectorXd F = VectorXd::Zero(m);

	const QuadRule m_qr = map_to_interval(qr, -PI / 2, PI / 2);

	// calculate the integral with the quadrature
	for(int i = 0; i < m_qr.nodes.size(); i++) {
		double s = m_qr.nodes(i);
		double x = std::tan(s);
		F += (m_qr.weights(i) * f(alpha, x) * (1 + x * x)) * rho(x);
	}

	F -= u;

	return F;
}


// routine: evalJ
// " Computes the Jacobian of $F(\alpha)$ "
// (in)  qr: quad nodes and weights
// (in)  alpha: model parameters
// (out) jacobian: Jacobian
MatrixXd evalJ(const QuadRule& qr, const VectorXd& alpha, const VectorXd& u) {

	int m = alpha.size();
	MatrixXd jacobian = MatrixXd::Zero(m, m);

	const QuadRule m_qr = map_to_interval(qr, -PI / 2, PI / 2);

	for(int i = 0; i < m_qr.nodes.size(); i++) {
		double s = m_qr.nodes(i);
		double x = std::tan(s);
		jacobian += (m_qr.weights(i) * f(alpha, x) * (1 + x * x)) * m_rho(x);
	}

	return jacobian;
}


// routine: newtonMethod
// " Solves the non-linear system using Newton method "
// (in)  qr: quad nodes and weights
// (in)  u: vector of moments
// (in)  atol: absolute tolerance
// (in)  rtol: relative tolerance
// (in)  maxItr: maximum iterations
// (in/out) alpha\_: model parameters
void newtonMethod(const QuadRule& qr, const VectorXd& u, const double atol, const double rtol, const int maxItr, VectorXd& alpha_) {

	for(int i = 0; i < maxItr; i++) {
		MatrixXd J = evalJ(qr, alpha_, u);
		VectorXd F = evalF(qr, alpha_, u);

		VectorXd correction_term = J.fullPivLu().solve(F);

		// what kind of stopping criterion is this ?
		// souldn't it just be: correction_term.norm() <= rtol * alpha_ like in the lecture notes p.249 ?
		if(correction_term.norm() <= atol || F.norm() <= rtol * u.norm()) {
			return;
		}
		
		alpha_ -= correction_term;
	}
}


// -- DO NOT CHANGE THIS ROUTINE --
// routine: testEvalF
// " test evalF computation "
// (in) alpha: model parameters
// (in) u: vector of moments
// (in) nQuad: number of quad nodes and weights
void testEvalF(const VectorXd& alpha, const VectorXd& u, const int nQuad) {

	QuadRule qr;
	gaussQuad(nQuad, qr);

	VectorXd F = evalF(qr, alpha, u);
	std::cout << "Norm of F: " << F.norm() << ", should be small" << std::endl;  
}


// routine: main
int main() {

	int nQ = 40; // number of Gauss points

	// -- DO NOT CHANGE THIS TEST SNIPPET --
	// test evalF
	// optimal model parameters corresponding to a moment vector are given
	VectorXd u_test(3), alpha_test(3);
	u_test << 1, 0, 1;
	alpha_test << -0.9189385332, 0, -0.5;
	testEvalF (alpha_test, u_test, nQ);

	// subproblem (h)
	VectorXd u(3);
	u << 1, -1, 2;

	VectorXd alpha(3);
	alpha << 0, 0, -0.2;

	double atol = 1e-8;
	double rtol = 1e-6;
	int maxItr = 100;

	QuadRule qr;
	gaussQuad(nQ, qr);
	newtonMethod(qr, u, atol, rtol, maxItr, alpha);

	std::cout << alpha.transpose() << std::endl;  
}


// END OF FILE
