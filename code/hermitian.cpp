#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>
#include <random>
#include "../../eigen/Eigen/Dense"
#include "../../eigen/Eigen/Eigenvalues"

using namespace Eigen;
using namespace std;

typedef long double ld;
typedef complex<double> co;
typedef MatrixXcd Mt;

int sort_k;

const ld eps = 1e-6;
const co I = {0.0, 1.0};

default_random_engine generator;
normal_distribution<double> distribution(0.0,1.0);

double rand_n() {
	return distribution(generator);
}

Mt abs_m(Mt A) {
	Mt B = A.adjoint()*A;
	SelfAdjointEigenSolver<Mt> es(B);
	return es.operatorSqrt();
}

double p_absx(double x, int k) {
	if (k == -1) {
		if (x > 0) {
			return 1.0;
		} else {
			return -1.0;
		}
	}
	double ans = abs(x);
	for (int i = 0; i < k; i++) {
		ans *= x;
	}
	return ans;
}

Mt xp_absx(Mt A, int k) {
	Mt B = abs_m(A);
	for (int i = 0; i < k; i++) {
		B = A*B;
	}
	return B;
}

double xp_absx_trace(Mt A, int k) {
	SelfAdjointEigenSolver<Mt> es(A);
	double ans = 0;
	Mt eigs = es.eigenvalues();
	for (int i = 0; i < (eigs.rows()); i++) {
		ans += p_absx(eigs(i).real(), k);
	}
	return ans;
}

Mt randH(int n) {
	Mt X(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			X(i, j) = rand_n() + I*rand_n();
		}
	}
	return X + X.adjoint();
}

Mt randHP(int n) {
	Mt X = randH(n);
	return X*X;
}

bool is_positive(Mt A) {
	LLT<Mt> lltOfA(A);
	if(lltOfA.info() == NumericalIssue) {
		return false;
	}
	return true;
}

double k_coeff(int k, int i) {
	if (k == i) {
		return 1.0;
	}
	if (i == 0) {
		if ((k%2) == 0) {
			return 1.0;
		} else {
			return -1.0;
		}
	}
	return k_coeff(k - 1, i - 1) - k_coeff(k - 1, i);
}

int binom(int n, int m) {
	if ((m == 0)||(n == m)) {
		return 1;
	}
	if (n < m) {
		return 0;
	}
	return binom(n - 1, m) + binom(n - 1, m - 1);
}

double divided_difference(vector<double> xs, int ala, int yla, int k) {
	sort(xs.begin() + ala, xs.begin() + yla + 1);
	if (xs[ala] + eps > xs[yla]) {
		if (xs[ala] < 0) {
			return 0;
		} else {
			int der = yla - ala;
			return pow(xs[ala], k - 1 - der)*binom(k - 1, der);
		}
	}
	else {
		return (divided_difference(xs, ala, yla - 1, k) - divided_difference(xs, ala + 1, yla, k))/(xs[ala] - xs[yla]);
	}
}

void go(vector<double> xs, Mt B, int k, vector<int> cur, ArrayXXcd &D) {
	int n = B.rows();
	if ((int)cur.size() == k + 1) {
		co val = 1;
		for (int i = 0; i < k; i++) {
			val *= B(cur[i], cur[i + 1]);
		}
		vector<double> ys;
		for (int i = 0; i < k + 1; i++) {
			ys.push_back(xs[cur[i]]);
		}
		val *= divided_difference(ys, 0, k, k);
		D(cur[0], cur[k]) = D(cur[0], cur[k]) + val;
	} else {
		for (int i = 0; i < n; i++) {
			vector<int> ne = cur;
			ne.push_back(i);
			go(xs, B, k, ne, D);
		}
	}
}

ArrayXXcd matrix_derivative(vector<double> xs, Mt B, int k) {
	int n = B.rows();
	ArrayXXcd D(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			D(i, j) = 0;
		}
	}
	go(xs, B, k, {}, D);
	return D;
}

bool test_k_trace(int n, int k) {
	vector<double> xs;
	for (int i = 0; i < n; i++) {
		xs.push_back(rand_n());
	}
	Mt B = randHP(n);
	ArrayXXcd D = matrix_derivative(xs, B, k);
	double tes = 0;
	for (int i = 0; i < n; i++) {
		cout << D(i, i).real() << " ";
		tes += D(i, i).real();
	}
	cout << endl;
	return (tes > -eps);
}

double get_val(vector<double> xs, double t, int k) {
	double ans = 0;
	for (double x : xs) {
		ans += p_absx(x - t, k - 2);
	}
	return ans;
}

double get_trace(Mt A, Mt H, int k) {
	double tes = 0;
	for (int i = 0; i <= k; i++) {
		tes += k_coeff(k, i)*((xp_absx(A + i*H, k - 2).trace()).real());
	}
	return tes;
}

double get_trace2(Mt A, Mt H, int k) {
	double tes = 0;
	for (int i = 0; i <= k; i++) {
		tes += k_coeff(k, i)*xp_absx_trace(A + i*H, k - 2);
	}
	return tes;
}

bool try_ktone(int n, int k) {
	Mt A = randH(n);
	if (is_positive(A)) {
		return true;
	}
	Mt H = randHP(n);
	double tes = get_trace2(A, H, k);
	if (tes < -eps) {
		cout << "tes: " << tes << endl;
		cout << "A = " << endl;
		cout << A << endl;
		cout << "H = " << endl;
		cout << H << endl;
	}
	return tes > -eps;
}

double approx_largest_eigenvalue(Mt A) {
	double norm_1 = A.norm();
	int n = A.rows();
	int m = A.cols();
	Mt B = A + norm_1*Mt::Identity(n, m);
	return B.norm() - norm_1;
}

double approx_smallest_eigenvalue(Mt A) {
	return -approx_largest_eigenvalue(-A);
}

vector<double> get_k_tone_vec(int n, int k, double step) {
	Mt A = randH(n);
	double eig1 = approx_smallest_eigenvalue(A);
	Mt Id = Mt::Identity(n, n);
	A -= Id*eig1;
	Mt H = randHP(n);
	eig1 = approx_largest_eigenvalue(A + k*H);
	vector<double> traces;
	for (double cur = 0; cur < eig1; cur += step) {
		double nt = get_trace(A - Id*cur, H, k);
		traces.push_back(nt);
	}
	return traces;
}

void test_conjecture(int k, int n, int cap) {
	for (int i = 0; i < cap; i++) {
		if (!try_ktone(n, k)) {
			cout << "FAIL: " << i << "\n";
			return;
		}
	}
	cout << "OK\n";
	return;
}

void test_conjecture2(int n, int k, int cap) {
	for (int i = 0; i < cap; i++) {
		if (!test_k_trace(n, k)) {
			cout << "FAIL: " << i << "\n";
			return;
		}
	}
	cout << "OK\n";
	return;
}

void test_random(int k, int n, int cap) {
	for (int i = 0; i < cap; i++) {
		int kk = rand()%k + 1;
		int nn = rand()%n + 1;
		if (!try_ktone(nn, kk)) {
			cout << "FAIL: " << i << "\n";
			cout << "n = " << nn << ", k = " << kk << "\n";
			return;
		}
	}
	cout << "OK\n";
	return;
}

bool cmp_tuples(vector<double> xs, vector<double> ys) {
	bool less = false;
	bool greater = false;
	for (double t : xs) {
		double d1 = get_val(xs, t, sort_k);
		double d2 = get_val(ys, t, sort_k);
		if (d1 + eps < d2) {
			less = true;
		} else if (d1 > d2 + eps) {
			greater = true;
		}
	}
	for (double t : ys) {
		double d1 = get_val(xs, t, sort_k);
		double d2 = get_val(ys, t, sort_k);
		if (d1 + eps < d2) {
			less = true;
		} else if (d1 > d2 + eps) {
			greater = true;
		}
	}
	if (less&&greater) {
		cout << "FAIL\n";
		exit(0);
	}
	return less;
}

vector<double> get_3_tuple() {
	//Generate 3 random numbers with sum 0 and sum of squares 1, and return them ordered in increasing order
	double a = rand_n();
	double b = rand_n();
	double c = rand_n();
	double s3 = (a + b + c)/3.0;
	a -= s3;
	b -= s3;
	c -= s3;
	double qs = sqrt(a*a + b*b + c*c);
	if (qs < eps) {
		return get_3_tuple();
	}
	a /= qs;
	b /= qs;
	c /= qs;
	vector<double> ans = {a, b, c};
	sort(ans.begin(), ans.end());
	return ans;
}

void try_3_tuples(int cnt) {
	sort_k = 3;
	vector<vector<double>> tupls;
	for (int i = 0; i < cnt; i++) {
		tupls.push_back(get_3_tuple());
	}
	sort(tupls.begin(), tupls.end(), cmp_tuples);
	int ind = 0;
	for (vector<double> vs : tupls) {
		double prod = 1;
		for (double x : vs) {
			prod *= x;
		}
		cout << ind << " " << prod;
		ind++;
		cout << endl;
	}
}

void print_k_tone_vec(int n, int k, double step) {
	vector<double> v = get_k_tone_vec(n, k, step);
	cout << setprecision(15);
	for (int i = 0; i < (int) v.size(); i++) {
		cout << i*step << " " << v[i] << "\n";
	}
}

int main() {
	int cap = 100;
	int n, k;
	cin >> n >> k;
	test_conjecture2(n, k, cap);
	return 0;
}