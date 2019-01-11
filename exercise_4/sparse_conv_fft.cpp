#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <string>
#include <complex>

using namespace std;

const complex<double> I(0,1); // imaginary unit
const double PI = 3.14159265359;

template<class T> struct duplet {
	int ind;
	T val;

	duplet(int p, T v) {
		ind = p;
		val = v;
	}
};

template<class T> struct sparse_vec {
	double tol = 1e-6;
	vector<duplet<T> > duplets;
	int len=0;

	sparse_vec(int l) {
		len = l;
	}

	void append(int ind, T val) {
		if(std::abs(val) >= tol) {
			duplets.push_back(duplet<T>(ind, val));
		}
	}

	void cleanup() {
		std::sort(duplets.begin(), duplets.end(),
		[&](const struct duplet<T> &a, const struct duplet<T> &b) {
			return a.ind < b.ind;
		});

		std::remove_if(duplets.begin(), duplets.end(),
		[=](const struct duplet<T> &a) {
			return (std::abs(a.val) < tol) || a.ind > (len - 1);
		});

		if(duplets.size() == 0) { return; }

		vector<duplet<T> > uniq;
		uniq.push_back(*duplets.begin());

		for(auto it = duplets.begin() + 1; it != duplets.end(); it++) {
			struct duplet<T> *last = &uniq.back();

			if(it->ind == last->ind) {
				last->val += it->val;
			}
			else {
				uniq.push_back(*it);
			}
		}

		duplets = uniq;
	}

	T get_val(int ind) const {
		int lo = 0;
		int hi = duplets.size() - 1;

		while(lo <= hi) {
			int mid = lo + (hi - lo) / 2;

			if(duplets[mid].ind > ind) {
				hi = mid - 1;
			}
			else if(duplets[mid].ind < ind) {
				lo = mid + 1;
			}
			else {
				return duplets[mid].val;
			}
		}

		return 0;
	}

	static sparse_vec cwise_mult(const sparse_vec &a, const sparse_vec &b) {
		sparse_vec out(max(a.len,b.len));

		auto a_it = a.duplets.begin();
		auto b_it = b.duplets.begin();

		while(a_it != a.duplets.end() && b_it != b.duplets.end()) {
			if(a_it->ind < b_it->ind) {
				a_it++;
			}

			if(a_it->ind > b_it->ind) {
				b_it++;
			}

			if(a_it->ind == b_it->ind) {
				out.append(a_it->ind, a_it->val * b_it->val);
				a_it++;
				b_it++;
			}
		}

		return out;
	}

	static sparse_vec conv(const sparse_vec &a, const sparse_vec &b) {
		sparse_vec out(a.len + b.len - 1);

		for(auto a_it = a.duplets.begin(); a_it != a.duplets.end(); a_it++) {
			for(auto b_it = b.duplets.begin(); b_it != b.duplets.end(); b_it++) {
				out.append(a_it->ind + b_it->ind, a_it->val * b_it->val);
			}
		}

		out.cleanup();
		return out;
	}

	static sparse_vec fft(const sparse_vec &x) {
		int n = x.len;
		sparse_vec tot(n);

		if(n == 1) { return x; }

		int m = n / 2;
		sparse_vec x_even(m);
		sparse_vec x_odd(m);

		for(auto it = x.duplets.begin(); it != x.duplets.end(); it++) {
			if(it->ind % 2 == 0) {
				x_even.append(it->ind / 2, x.get_val(it->ind));
			}
			else {
				x_odd.append((it->ind - 1) / 2, x.get_val(it->ind));
			}
		}

		sparse_vec c_even = fft(x_even);
		sparse_vec c_odd = fft(x_odd);

		std::complex<double> omega = std::exp(-2 * PI / n * I);

		for(int k = 0; k < n; k++) {
			tot.append(k, c_even.get_val(k % m) + std::pow(omega, k) * c_odd.get_val(k % m));
		}

		return tot;
	}

	static sparse_vec ifft(const sparse_vec &x) {
		double n = x.len;
		sparse_vec out(n);

		sparse_vec x_conj(n);
		for(auto x_it = x.duplets.begin(); x_it != x.duplets.end(); x_it++) {
			x_conj.append(x_it->ind, std::conj(x_it->val));
		}

		sparse_vec fft_conj = fft(x_conj);
		for(auto fft_it = fft_conj.duplets.begin(); fft_it != fft_conj.duplets.end(); fft_it++) {
			out.append(fft_it->ind, std::conj(fft_it->val) / n);
		}

		return out;
	}

	static sparse_vec conv_fft(sparse_vec a, sparse_vec b) {
		int n = a.len + b.len - 1;
		a.len = n;
		b.len = n;
		return ifft(cwise_mult(fft(a), fft(b)));
	}

	std::string to_string() const {
		std::stringstream ss;
		for (auto p : this->duplets) {
			ss << "(" << p.ind << "," << p.val << "),";
		}
		ss << "\n";
		std::string out = ss.str();
		return out;
	}


};




/***** TESTING ******/

int main() {

	sparse_vec<complex<double> > x(5);
	x.append(0,complex<double>(8.2,0));
	x.append(1,complex<double>(1,-2));
	x.append(3,complex<double>(-3,4.66));
	x.append(4,complex<double>(0,4));
	x.cleanup();

	sparse_vec<complex<double> > y(4);
	y.append(1,complex<double>(5,0));
	y.append(2,complex<double>(1.21,-4));
	y.append(3,complex<double>(4,2.4));
	y.cleanup();

	auto m = sparse_vec<complex<double> >::cwise_mult(x,y);
	m.cleanup();
	cout << "TESTS. Correct componentwise multiplication between x and y: (1,(5,-10)),(3,(-23.184,11.44)),\n";
	cout << "cwise_mult(x,y) = " << m.to_string();

	auto c = sparse_vec<complex<double> >::conv(x,y);
	c.cleanup();
	cout << "Correct exact discrete convolution between x and y: (1,(41,0)),(2,(14.922,-42.8)),(3,(26.01,13.26)),(4,(-6.2,17.7)),(5,(15.01,37.6386)),(6,(-7.184,16.28)),(7,(-9.6,16)),\n";
	cout << "conv(x,y) = " << c.to_string();
	auto cf = sparse_vec<complex<double> >::conv_fft(x,y);
	cf.cleanup();
	cout << "conv_fft(x,y) = " << cf.to_string();
}
