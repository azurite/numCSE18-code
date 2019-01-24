#include <vector>
#include <iostream>
#include <cmath>
#include "writer.hpp"

double f(double t, double u) {
    return std::exp(-2*t) - 2*u;
}

double exact(double t) {
    return t * std::exp(-2*t);
}


/// Uses the SSP RK3 method to compute u from time 0 to time T 
/// for the ODE $u'=e^{-2t}-2u$
///
/// @param[out] u at the end of the call (solution at all time steps)
/// @param[out] time contains the time levels
/// @param[in] u0 the initial data
/// @param[in] dt the step size
/// @param[in] T the final time up to which to compute the solution.
///

//----------------SSPRK3Begin----------------
void SSPRK3(std::vector<double> & u, std::vector<double> & time,
          const double & u0, double dt, double T) {
	const unsigned int nsteps = std::round(T/dt);
	u.resize(nsteps+1);
	time.resize(nsteps+1);

	u[0] = u0;
	time[0] = 0;
	
	for(int k = 0; k < nsteps - 1; k++) {
		time[k + 1] = time[k] + dt;
		double tk = k * dt;		
		double k1 = f(tk, u[k]);
		double k2 = f(tk + dt, u[k] + dt * k1);
		double k3 = f(tk + 0.5 * dt, u[k] + dt * 0.25 * (k1 + k2));
		u[k + 1] = u[k] + dt * ((1 / (double)6) * k1 + (1 / (double)6) * k2 + (2 / (double)3) * k3);
	}


}
//----------------SSPRK3End----------------

int main(int argc, char** argv) {

	double T = 10.0;

	double dt = 0.5;

	// To make some plotting easier, we take the dt parameter in as an optional
	// parameter.
	if (argc == 2) {
	dt = atof(argv[1]);
	}

	const double u0 = 0.;
	std::vector<double> error;
	std::vector<double> h = { 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625 };

	for(int k = 0; k < 8; k++) {
		std::vector<double> time;
		std::vector<double> u;
		SSPRK3(u, time, u0, h[k], T);
		error.push_back(std::abs(u.back() - exact(T)));
	}

	writeToFile("error.txt", error);
	writeToFile("dt.txt", h);

	return 0;
}
