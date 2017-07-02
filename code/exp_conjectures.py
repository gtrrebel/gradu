from numpy.random import randn
from numpy.random import randint
from math import exp
import matplotlib.pyplot as plt
import numpy as np
eps = 1e-3
from numpy.linalg import matrix_power as mat_pow

ass = []
B = [[]]
res = 0

def get_p():
	return sum(1 for x in ass if x < 0)

def fact(n):
	if (n < 0):
		return 0
	if (n == 0):
		return 1
	return n*fact(n - 1)

def binom(n, k):
	if (n < 0):
		return 0
	if (k > n):
		return 0
	if (k == 0):
		return 1
	return binom(n - 1, k) + binom(n - 1, k - 1)

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

def get_B_coeff(l):
	coeff = 1
	for i in xrange(len(l)):
		coeff *= B[l[i - 1]][l[i]]
	return coeff


def get_basic_exp_divided(xx, nn):
	return exp(-abs(xx))/fact(nn - 1)

def get_exp_divided_difference(l):
	if (len(l) == 0):
		return 0
	l = sorted(l)
	if (l[0] < l[-1]):
		ans1 = get_exp_divided_difference(l[:-1])
		ans2 = get_exp_divided_difference(l[1:])
		return (ans1 - ans2)/(ass[l[0]] - ass[l[-1]])
	else:
		return get_basic_exp_divided(ass[l[0]], len(l))

def process(l, n, k):
	B_coeff = get_B_coeff(l)
	p = get_p()
	l1 = [x for x in l if x < p]
	l2 = [x for x in l if x >= p]
	div_dif1 = get_exp_divided_difference(l1)
	div_dif2 = get_exp_divided_difference(l2)
	global res
	res += B_coeff*div_dif1*div_dif2
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

def try_conjecture_0(n, k, cap = 10000):
	for _ in xrange(cap):
		global res
		res = 0
		global ass, B
		B = get_B(n, k)
		ass = get_ass(n)
		go([], n, k)
		#print _
		#print ass
		#print res
		if (res < -eps):
			print "FAIL"
			print res
			print ass
			exit(0)


try_conjecture_0(4, 3)



