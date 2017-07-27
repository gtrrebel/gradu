from numpy.random import randn, randint, uniform
import matplotlib.pyplot as plt
import numpy as np
eps = 1e-6
import cmath

step_len = 5e-3
iter_count = 5000

def get_B(n):
	B = np.zeros((n, n), dtype=complex)
	for i in xrange(n):
		for j in xrange(n):
			B[i][j] = randn()
	B = np.dot(B, np.transpose(B))
	return B

def get_ass(n):
	return sorted(randn(n))

def get_A(ass, z):
	A = np.zeros((n, n), dtype=complex)
	for i in xrange(n):
		A[i][i] = 1.0/(1.0 - ass[i]*z)
	return A

def get_term(ass, B, z, k, extra = False):
	n = len(ass)
	A = get_A(ass, z)
	C = np.eye(n, dtype=complex)
	for i in xrange(k):
		C = np.dot(C, A)
		C = np.dot(C, B)
	if (extra):
		C = np.dot(C, A)
	T = np.trace(C)
	di = -np.conj(T)
	le = abs(di)
	di /= le
	return di

def get_term2(ass, B, z, k, extra = False):
	ans = (z**2 - 1)/z
	ans /= abs(ans)
	return ans

def get_trajectory(ass, B, k, step = step_len, iters = iter_count):
	z0 = -0.5 + 0.25j
	zs = [z0]
	for ite in xrange(iters):
		print ite*1.0/iters
		di = get_term2(ass, B, z0, k)
		z0 += di*step
		zs.append(z0)
	return zs

def get_random_curve(n, k):
	ass = get_ass(n)
	B = get_B(n)
	zs = get_trajectory(ass, B, k)
	return zs

def plot_curve(ass, zs):
	xs = [np.real(z) for z in zs]
	ys = [np.imag(z) for z in zs]
	for a in ass:
		plt.plot([a], [0], "r*")
	plt.plot(xs, ys)
	plt.show()

def plot_random_curve(n, k):
	print "jee"
	ass = get_ass(n)
	if (ass[0]*ass[1] > 0):
		return
	B = get_B(n)
	zs = get_trajectory(ass, B, k)
	plot_curve(ass, zs)

n = 3
k = 4

print plot_random_curve(n, k)

