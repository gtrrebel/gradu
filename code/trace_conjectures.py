from numpy.random import randn
from numpy.random import randint
import matplotlib.pyplot as plt
import numpy as np
eps = 1e-3
from numpy.linalg import matrix_power as mat_pow

def binom(n, k):
	if (n < 0):
		return 0
	if (k > n):
		return 0
	if (k == 0):
		return 1
	return binom(n - 1, k) + binom(n - 1, k - 1)

def get_H(n):
	B = np.zeros((n, n))
	for _ in xrange(n):
		for __ in xrange(n):
			B[_][__] = randn(1)[0]
	B = B + B.transpose()
	return B

def get_P(n):
	B = get_H(n)
	return np.dot(B,np.transpose(B))

def get_M(n, m):
	C = np.zeros((n, m))
	for _ in xrange(n):
		for __ in xrange(m):
			C[_][__] = randn(1)[0]
	return C

def try_conjecture_0(n, m, cap = 10000):
	for _ in xrange(cap):
		M_m = get_H(n)
		M_p = get_H(m)
		C = get_M(n, m)
		Cp = np.transpose(C)
		terms = []
		for i in xrange(3):
			M1 = mat_pow(M_m, 2 - i)
			M2 = mat_pow(M_p, i)
			MM = np.dot(M1, C)
			MM2 = np.dot(M2, Cp)
			MM3 = np.dot(MM, MM2)
			terms.append(np.trace(MM3)*binom(2, i))
		ans = sum(terms)
		if (ans < -eps):
			print "FAIL"
			print ans
			exit(0)

def try_conjecture_1(n, m, cap = 10000):
	for _ in xrange(cap):
		M_m = get_H(n)
		M_p = get_H(m)
		C = get_M(n, m)
		Cp = np.transpose(C)
		terms = []
		for i in xrange(5):
			M1 = mat_pow(M_m, 4 - i)
			M2 = mat_pow(M_p, i)
			MM = np.dot(M1, C)
			MM2 = np.dot(M2, Cp)
			MM3 = np.dot(MM, MM2)
			terms.append(np.trace(MM3)*binom(4, i))
		ans = sum(terms)
		if (ans < -eps):
			print "FAIL"
			print ans
			exit(0)


try_conjecture_0(3, 2)



