from numpy.random import randn
from numpy.random import randint
import matplotlib.pyplot as plt
import numpy as np
eps = 1e-3

# x is always 0

B = [[]]

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
	return sum(1 for x in l if x <= p)

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

def get_weird_term():
	As = [B[0][0]**2*B[0][a] for a in xrange(1, n)]
	Bs = [sum(B[0][0]*B[0][b]*B[b][a] for b in xrange(1, n)) for a in xrange(1, n)]
	Cs = [sum(sum(B[0][b]*B[b][c]*B[c][a] for c in xrange(1, n)) for b in xrange(1, n)) for a in xrange(1, n)]
	term1 = sum(A**2 for A in As)
	term2 = 4*sum(A*B for A, B in zip(As, Bs))
	term3_1 = 6*sum(A*C for A, C in zip(As, Cs))
	term3_2 = 6*sum(B*B for B in Bs)
	if (abs(term3_1 - term3_2) > eps):
		print "FAIL",
		print term3_1, term3_2
		exit(0)
	term3 = term3_1
	term4 = 4*sum(B*C for B, C, in zip(Bs, Cs))
	term5 = sum(C*C for C in Cs)
	terms = [term1, term2, term3, term4, term5]
	term = sum(terms)
	if (term < -eps):
		print "FAIL",
		print term
		exit(0)
	return term

def test_weird_terms(n, cap = 10000):
	for _ in xrange(cap):
		global B
		B = get_B(n, 6)
		t = get_weird_term()
		if (t < -eps):
			print "FAIL"
			print t
			return

n = 4
k = 6
test_weird_terms(n)



