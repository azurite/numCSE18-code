#include<bits/stdc++.h>

using namespace std;

float naive_sum(vector<float> &v) {
	float sum = 0;
	for(float x : v){
		sum += x;
	}
	return sum;
}

float sum(vector<float> v) {
	float sum = 0;

	std::sort(v.begin(), v.end(), [](float a, float b) {
		return std::abs(a) < std::abs(b);
	});

	// we add numers with small magnitude first
	for(float x : v) {
		sum += x;
	}

	return sum;
}

float acc_sum(vector<float> &v) {
	float sum = 0;
	float err = 0;

	for(float x : v) {
		/*
		 * case 1: x and err do not lead to cancellation
		 * new_error = x - (((sum + x + err) - sum) - err) = x - (x + err - err) = x - x = 0
		 *
		 * case 2: only err get's cancelled
		 * new_error = x - (((sum + x + err) - sum) - err) = x - ((sum + x - sum) - err) = x - (x - err) = err
		 *
		 * case 3: only x get's cancelled
		 * new_error = x - (((sum + x + err) - sum) - err) = x - ((sum + err - sum) - err) = x - (err - err) = x
		 *
		 * case 4: x and err lead to cancellation
		 * new_error = x - (((sum + x + err) - sum) - err) = x - ((sum - sum) - err) = x - (0 - err) = x + err
		 *
		 * Therefore in each case the error conatinas exactly the terms that were cancelled and we add that new_error
		 * back into the sum
		 */
		err = x - (((sum + x + err) - sum) - err);
		sum = sum + x + err;
	}

	return sum;
}

int main() {
	srand(time(0));
	cout << setprecision(15);
	int N = 1e7;
	double corr_sum = 0;
	vector<float> v(N);
	for (int i = 0; i < N; i++) {
		double x = 1e-8 + (rand() % 10) * 1e-9;
		corr_sum += x;
		v[i] = (float) x;
	}
	corr_sum += 1;
	v[1]=1;

	cout << "naive_sum(v) = " << naive_sum(v) << endl
		 << "      sum(v) = " << sum(v) << endl
		 << "  acc_sum(v) = " << acc_sum(v) << endl
		 << " correct sum = " << corr_sum << endl;
}
