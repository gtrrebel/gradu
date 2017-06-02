from numpy.random import randn
import matplotlib.pyplot as plt
import numpy as np
eps = 1e-3

def sum_coeff(n, k):
	if (n < 0):
		return 0
	if (n < k):
		return 0
	if (k == 0):
		return 1
	return n*sum_coeff(n - 1, k - 1)/k

def go(i, n, k, l, subs):
	if (i == n):
		if (len(l) == k):
			subs.append(l)
		return
	if (len(l) > k):
		return
	l1 = [x for x in l]
	l2 = [x for x in l]
	l2.append(i)
	go(i + 1, n, k, l1, subs)
	go(i + 1, n, k, l2, subs)

def subsets_of_size_k(n, k):
	subs = []
	go(0, n, k, [], subs)
	return subs

def get_abc(n = 3):
	ass = randn(n)
	return sorted(ass)

def is_zero(tt):
	return (abs(tt) < eps)

def is_pos(tt):
	return (tt > -eps)

def get_B(n = 3):
	B = np.zeros((n, n))
	for _ in xrange(n):
		for __ in xrange(n):
			B[_][__] = randn(1)[0]
	return np.dot(B,np.transpose(B))

def get_x():
	return randn()

def get_formula(ass, B):
	a, b, c = ass
	ans = 0
	all_prod = B[0][1]*B[1][2]*B[2][0]
	#Negative terms
	ans += 2*B[0][0]*all_prod*(c - b)*(c - a)
	###ans += 2*B[1][1]*all_prod*(c - a)**2*(b - a)
	ans += 2*B[2][2]*all_prod*(c - a)*(b - a)
	#Odd term
	ans += 2*B[0][0]*B[0][2]**2*B[2][2]*(c - b)*(b - a)
	#Split terms
	##ans += B[0][1]**2*B[0][2]**2*(c - b)*(c - a)*(b - a)
	ans += B[0][1]**2*B[1][2]**2*(c - a)**2
	###ans += B[0][2]**2*B[1][2]**2*(c - a)*(b - a)**2
	#Half terms
	ans += B[0][0]**2*B[0][2]**2*(c - b)**2
	ans += B[2][2]**2*B[0][2]**2*(b - a)**2
	###ans += B[1][1]**2*B[0][1]**2*(c - a)**3
	#Four term
	##ans += B[0][2]**4*(c - b)*(b - a)**2
	return ans

def random_test():
	ass = get_abc()
	B = get_B()
	ans = get_formula(ass, B)
	if (ans < -1):
		print "FAIL"
		print ans
		print
		print ass
		print B

def do_tests(N = 10):
	for _ in xrange(N):
		random_test()

def get_sub_term(ass, B, js):
	ass_2 = []
	n = len(ass)
	ass_2 = [ass[j] for j in js]
	B_2 = [[B[j1][j2] for j2 in js] for j1 in js]
	return ass_2, B_2


def get_general_term(ass, B, x, i, first = True, second = True):
	n = len(ass)
	term1, term2, term3 = 0, 0, 0
	for j in xrange(n):
		if (j != i):
			insum = 0
			for l in xrange(n):
				if (l != i):
					insum += (x - ass[i])/(ass[l] - ass[i])*B[i][l]*B[l][j]
			insum += (x - ass[j])/(ass[i] - ass[j])*B[i][i]*B[i][j]
			term1 += insum*insum/(ass[j] - ass[i])
	for j in xrange(n):
		if (j != i):
			term2 += (ass[j] - x)/((ass[j] - ass[i])**2)*B[i][j]**2
	for j in xrange(n):
		if (j != i):
			term3 += B[i][j]**2/(ass[j] - ass[i])
	ans = 0
	if first:
		ans += 4*term1
	if second:
		ans += 4*term2*term3*(x - ass[i])
	return ans

def try_general_sum_4(ass, B, x):
	n = len(ass)
	if (x < ass[0]):
		return 0
	if (x > ass[-1]):
		return 0
	p = 0
	while (ass[p] < x):
		p += 1
	ans = 0
	for i in xrange(p):
		for j in xrange(p, n):
			term = 0
			term1 = 0
			for i1 in xrange(p):
				term1 += B[i][i1]*B[i1][j]/(ass[j] - ass[i1])
			term2 = 0
			for j1 in xrange(p, n):
				term2 += B[i][j1]*B[j1][j]/(ass[i] - ass[j1])
			term = term1*(x - ass[j]) + term2*(x - ass[i])
			ans += term**2/(ass[j] - ass[i])
	for i1 in xrange(p):
		termij_0 = (x - ass[i1])
		for i2 in xrange(p):
			for j in xrange(p, n):
				termij_1 = termij_0*(x - ass[j])/((ass[j] - ass[i1])*(ass[j] - ass[i2]))
				termij_1 *= B[i1][j]*B[i2][j]
				for l in xrange(p, n):
					termij_2 = termij_1*B[i1][l]*B[i2][l]/(ass[l] - ass[i1])
					#ans -= termij_2
	for i1 in xrange(p):
		for i2 in xrange(p):
			interm = 0
			for j in xrange(p, n):
				interm += (x - ass[j])/((ass[j] - ass[i1])*(ass[j] - ass[i2]))*B[i1][j]*B[i2][j]
			ans += interm**2*(x - ass[i1])
	for j1 in xrange(p, n):
		for j2 in xrange(p, n):
			interm = 0
			for i in xrange(p):
				interm += (x - ass[i])/((ass[i] - ass[j1])*(ass[i] - ass[j2]))*B[i][j1]*B[i][j2]
			ans += interm**2*(ass[j1] - x)
	return 4*ans


def get_3_general_term(ass, B, x, i):
	n = len(ass)
	term1, term2 = 0, 0
	for j in xrange(n):
		if (j != i):
			term1 += (ass[j] - x)/((ass[j] - ass[i])**2)*B[i][i]*B[i][j]**2
	for j1 in xrange(n):
		if (j1 != i):
			for j2 in xrange(n):
				if (j2 != i):
					term2 += B[i][j1]*B[i][j2]*B[j1][j2]/((ass[j1] - ass[i])*(ass[j2] - ass[i]))
	ans = 0
	ans += 3*term1
	ans += 3*term2*(x - ass[i])
	return ans

def get_sum_sum(ass, B, x, kk):
	n = len(ass)
	subs = subsets_of_size_k(n, kk)
	ans = 0
	for sub in subs:
		ass_2, B_2 = get_sub_term(ass, B, sub)
		for i in xrange(len(ass_2)):
			if (ass_2[i] < x):
				ans += get_general_term(ass_2, B_2, x, i)
	return ans/sum_coeff(n, kk)

def get_sums(ass, B, x):
	n = len(ass)
	for kk in xrange(n + 1):
		print kk, ":",
		print get_sum_sum(ass, B, x, kk), "  ",
	print

def test_sum_sum_right(ass, B, cap = 1):
	print ass
	for _ in xrange(cap):
		x = get_x()
		print x
		get_sums(ass, B, x)

def test_sum_conj(n = 4, cap = 1):
	for _ in xrange(cap):
		ass = get_abc(n)
		B = get_B(n)
		test_sum_sum_right(ass, B)

def test_sum_right(ass, B, cap = 1000):
	n = len(ass)
	for _ in xrange(cap):
		x = get_x()
		ts = []
		for i in xrange(n):
			ts.append(get_general_term(ass, B, x, i))
		su = 0
		for i in xrange(n):
			su += ts[i]
		ok = is_zero(su)
		if (not ok):
			print "FAIL"
			print su
			print ts
			print ass
			print B

def try_sum_4_conjecture(ass, B, cap = 100):
	n = len(ass)
	for _ in xrange(cap):
		x = get_x()
		ts = []
		for i in xrange(n):
			ts.append(get_general_term(ass, B, x, i))
		su = 0
		for i in xrange(n):
			if (ass[i] < x):
				su += ts[i]
		su2 = try_general_sum_4(ass, B, x)
		if (abs(su2 - su) > eps):
			print "FAIL"
			print su
			print x
			print ass

def test_3_sum_right(ass, B, cap = 1000):
	n = len(ass)
	for _ in xrange(cap):
		x = get_x()
		ts = []
		for i in xrange(n):
			ts.append(get_3_general_term(ass, B, x, i))
		su = 0
		for i in xrange(n):
			su += ts[i]
		ok = is_zero(su)
		if (not ok):
			print "FAIL"
			print su
			print ts
			print ass
			print B

def test_pos(ass, B, cap = 1000):
	n = len(ass)
	for _ in xrange(cap):
		x = get_x()
		ts = []
		for i in xrange(n):
			ts.append(get_general_term(ass, B, x, i))
		su = 0
		for i in xrange(n):
			if (ass[i] < x):
				su += ts[i]
		ok = is_pos(su)
		if (not ok):
			print "FAIL"
			print su
			print ts
			print ass
			print B

def test_pos_3(ass, B, cap = 1000):
	n = len(ass)
	for _ in xrange(cap):
		x = get_x()
		ts = []
		for i in xrange(n):
			ts.append(get_3_general_term(ass, B, x, i))
		su = 0
		for i in xrange(n):
			if (ass[i] < x):
				su += ts[i]
		ok = is_pos(su)
		if (not ok):
			print "FAIL"
			print su
			print ass
			print B
			print x
			print ts

def get_pol_vals(ass, B, cap = 1000):
	n = len(ass)
	xs = np.linspace(ass[0] - 1, ass[-1] + 1, cap)
	vals = []
	for x in xs:
		ts = []
		for i in xrange(n):
			ts.append(get_general_term(ass, B, x, i))
		su = 0
		for i in xrange(n):
			if (ass[i] < x):
				su += ts[i]
		vals.append(su)
	return xs, vals

def print_random_vals(n, cap = 100):
	print n, cap
	ass = get_abc(n)
	B = get_B(n)
	for i in xrange(n):
		print ass[i],
	print
	for i in xrange(n):
		for j in xrange(n):
			print B[i][j],
		print
	xs, vals = get_pol_vals(ass, B, cap=cap)
	for xx, yy in zip(xs, vals):
		print xx, yy

def plot_random_vals(n, cap = 100):
	ass = get_abc(n)
	B = get_B(n)
	xs, vals = get_pol_vals(ass, B, cap=cap)
	plt.plot(xs, vals)
	max_val = max(vals)
	for i in xrange(n):
		plt.plot([ass[i], ass[i]], [0, max_val])
	plt.show()

def try_sums(n):
	cap = 100
	for _ in xrange(cap):
		ass = get_abc(n)
		B = get_B(n)
		test_3_sum_right(ass, B)

def try_pos(n):
	cap = 100
	for _ in xrange(cap):
		ass = get_abc(n)
		B = get_B(n)
		test_pos_3(ass, B)

def try_conj(n):
	cap = 100
	for _ in xrange(cap):
		ass = get_abc(n)
		B = get_B(n)
		try_sum_4_conjecture(ass, B)

n = 7

try_conj(n)
