from numpy.random import randn
from numpy.random import randint
import matplotlib.pyplot as plt
import numpy as np
eps = 1e-3

# x is always 0
ass = []
B = [[]]
res_list = [[]]
equal = []

def get_ass(n):
	ass = randn(n)
	return sorted(ass)

def get_B(n, k = -1):
	B = np.zeros((n, n))
	for _ in xrange(n):
		for __ in xrange(n):
			B[_][__] = randn(1)[0]
	B = B + B.transpose()
	if (k >= 0):
		if ((k%2) == 1):
			B = np.dot(B,np.transpose(B))
	return B

def get_random_l(n, k):
	l = []
	for _ in xrange(k):
		l.append(randint(0, n - 1))
	return l

def list_index(l):
	ind = 0
	for i in xrange(len(l)):
		x1 = ass[l[i - 1]]
		x2 = ass[l[i]]
		if (x1*x2 < 0):
			ind += 1
	ind //= 2
	return ind

def list_neg(l):
	p = get_p()
	return sum(1 for x in l if x < p)

def get_p():
	return sum(1 for x in ass if x < 0)

def get_B_coeff(l):
	coeff = 1
	for i in xrange(len(l)):
		coeff *= B[l[i - 1]][l[i]]
	return coeff

def binom(n, k):
	if (n < 0):
		return 0
	if (k > n):
		return 0
	if (k == 0):
		return 1
	return binom(n - 1, k) + binom(n - 1, k - 1)

def get_basic_divided(xx, nn, k):
	if (xx < 0):
		return 0
	if (nn == k):
		return 0
	return xx**(k - 1 - nn)*binom(k - 2, nn - 1)

def get_divided_difference(l, k):
	l = sorted(l)
	if (l[0] < l[-1]):
		ans1 = get_divided_difference(l[:-1], k)
		ans2 = get_divided_difference(l[1:], k)
		return (ans1 - ans2)/(ass[l[0]] - ass[l[-1]])
	else:
		return get_basic_divided(ass[l[0]], len(l), k)

def get_k5_i2_p1(n):
	res1 = 0
	for j1 in xrange(1, n):
		for j2 in xrange(1, n):
			div = 0
			x1 = (ass[j1])/(ass[j1] - ass[0])
			x2 = (ass[j2])/(ass[j2] - ass[0])
			div += x1*x1 + x1*x2 + x2*x2
			div *= -ass[0]*5
			div /= (ass[j1] - ass[0])*(ass[j2] - ass[0])
			term = div*B[0][0]*B[0][j1]**2*B[0][j2]**2
			res1 += term
	res2 = 0
	for j1 in xrange(1, n):
		for j2 in xrange(1, n):
			for j3 in xrange(1, n):
				div = 0
				div += (ass[j1])/(ass[j1] - ass[0])
				div += (ass[j2])/(ass[j2] - ass[0])
				div += (ass[j3])/(ass[j3] - ass[0])
				div *= 5*ass[0]**2
				div /= (ass[j1] - ass[0])*(ass[j2] - ass[0])*(ass[j3]-ass[0])
				term = div*B[0][j1]*B[j1][j2]*B[0][j2]*B[0][j3]**2
				res2 += term
	res = res1 + res2
	if (res2 < 0):
		print "FAIL"
		print res2
		print ass
		print B
		exit(0)
	return res

def process(l, n, k):
	ind1 = list_index(l)
	ind2 = list_neg(l)
	B_coeff = get_B_coeff(l)
	div_dif = get_divided_difference(l, k)
	global res_list
	if (ind1 > 0):
		res_list[ind1 - 1][ind2 - 1] += B_coeff*div_dif
	return

def process_diag(l, n, k):
	B_coeff = get_B_coeff(l)
	l.append(l[0])
	div_dif = get_divided_difference(l, k + 1)
	global diag
	diag[l[0]] += B_coeff*div_dif
	return

def process_equal(l, n, k):
	B_coeff = get_B_coeff(l)
	global equal
	ind = list_index(l)
	if (ind < k):
		equal[ind - 1] += B_coeff*binom(k - 2, ind - 1)
	return

def go(l, n, k):
	if (len(l) == k):
		process(l, n, k)
		return
	for i in xrange(n):
		l2 = [xx for xx in l]
		l2.append(i)
		go(l2, n, k)
	return

def go_diag(l, n, k):
	if (len(l) == k):
		process_diag(l, n, k)
		return
	for i in xrange(n):
		l2 = [xx for xx in l]
		l2.append(i)
		go_diag(l2, n, k)
	return

def go_equal(l, n, k):
	if (len(l) == k):
		process_equal(l, n, k)
		return
	for i in xrange(n):
		l2 = [xx for xx in l]
		l2.append(i)
		go_equal(l2, n, k)
	return

def set_vals(n, k):
	global res_list
	res_list = [[0 for _ in xrange(k)] for __ in xrange(k//2)]
	go([], n, k)
	return res_list

def set_diagonal_vals(n, k):
	global diag
	diag = [0 for _ in xrange(n)]
	go_diag([], n, k)
	return diag

def set_equal_val(n, k):
	global equal
	equal = [0 for _ in xrange(k - 1)]
	go_equal([], n, k)
	return equal

def get_random_vals(n, k):
	global ass, B
	ass = get_ass(n)
	B = get_B(n, k)
	vals = set_vals(n, k)
	return vals

def get_random_diagonal_values(n, k):
	global ass, B
	ass = get_ass(n)
	B = get_B(n, k)
	vals = set_diagonal_vals(n, k)
	return vals

def get_random_equal_value(n, k):
	global ass, B
	ass = get_ass(n)
	B = get_B(n, k)
	val = set_equal_val(n, k)
	print "ok"
	print ass
	print val
	return val

def test_div_dif(n, k):
	for i in xrange(10000):
		global ass
		ass = get_ass(n)
		l = get_random_l(n, k)
		yy = get_divided_difference(l, k)
		if (yy < -eps):
			print "FAIL"
			print yy
			print ass

def print_table(table):
	for ll in table:
		print ass
		print ll, "sum: ", sum(ll)
		#print sum(ll)
		if (sum(ll) < -eps):
			print "FAIL"
			exit(0)

def test_conjecture(n, k, cap = 100):
	for _ in xrange(cap):
		vals = get_random_vals(n, k)
		print _
		print_table(vals)
		#for ll in vals:
		#	for xx in ll:
		#		if (xx < -eps):
		#			print "FAIL"
		#			print xx
		#			print B
		#			print ass
		#			print vals
		#			return

def test_if_diagonal_values_are_non_negative(n, k, cap = 100):
	for _ in xrange(cap):
		vals = get_random_diagonal_values(n, k)
		for i in xrange(n):
			if (vals[i] < -eps):
				print "FAIL"
				print vals
				print B
				print ass
				return

def test_if_only_two_distinct_lambdas(n, k, cap = 100):
	for _ in xrange(cap):
		vals = get_random_equal_value(n, k)
		for i in xrange(k - 1):
			if (vals[i] < -eps):
				print "FAIL"
				print vals
				print B
				print ass
				return

n = 4
k = 5
test_conjecture(n, k)
