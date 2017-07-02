from numpy.random import randn, randint, uniform
from math import cos, sin
from math import pi as PI
import matplotlib.pyplot as plt
import numpy as np
eps = 1e-4

def maj_1(v1, v2):
	v1 = sorted(v1)
	v2 = sorted(v2)
	return all(x < y + eps for (x, y) in zip(v1, v2))

def get_random_vec():
	di = uniform()*2*PI
	return [cos(di), sin(di)], [-sin(di), cos(di)]

def random_mat():
	A = np.zeros((2, 2))
	A[0][0] = randn()
	A[1][1] = randn()
	A[0][1] = randn()
	A[1][0] = A[0][1]
	return A

def gf(A, v):
	su = 0
	n = len(v)
	for i in xrange(n):
		for j in xrange(n):
			su += A[i][j]*v[i]*v[j]
	return su

def is_pos(A, print_det=False):
	if (A[0][0] < -eps):
		return False
	det = A[0][0]*A[1][1] - A[0][1]*A[1][0]
	if print_det:
		print det
	return A[0][0]*A[1][1] - A[0][1]*A[1][0] + eps > 0

def test_low_order(A, B):
	return is_pos(B - A)

def test_1_order(A, B, cap = 1000000):
	for _ in xrange(cap):
		v1, v2 = get_random_vec()
		apu1s = [gf(A, v1), gf(A, v2)]
		apu2s = [gf(B, v1), gf(B, v2)]
		if not maj_1(apu1s, apu2s):
			return False
	return True

def is_same_order(cap = 100):
	for _ in xrange(cap):
		A = [[1, 0], [0, 0]]
		B = random_mat()
		ok1 = bool(test_low_order(A, B))
		ok2 = bool(test_1_order(A, B))
		ok = (ok1 == ok2)
		if not ok:
			print _
			print "fail"
			print A
			print B
			print ok1, ok2
			is_pos(B - A, print_det=True)
			return False
	return True

print is_same_order()






