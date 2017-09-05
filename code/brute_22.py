from numpy.random import randn, randint, uniform
from numpy import linspace
import matplotlib.pyplot as plt
import numpy as np
eps = 1e-1

lambs = [0, 0]
B = [[0, 0], [0, 0]]
cs = [0, 0]
x = 0

Coeffs = [[0, 0], [0, 0]]

def get_B():
	n = 2
	B = np.zeros((n, n))
	for _ in xrange(n):
		for __ in xrange(n):
			B[_][__] = randn(1)[0]
	B = B + B.transpose()
	return B

def get_B2():
	n = 2
	B = np.zeros((n, n))
	for _ in xrange(n):
		for __ in xrange(n):
			if (_ == __):
				B[_][__] = 1.0
			else:
				B[_][__] = 1.0 - eps
	return B

def apu_fs(i):
	if (i == 0):
		return (x - lambs[1])/(lambs[0] - lambs[1])
	else:
		return (x - lambs[0])/(lambs[1] - lambs[0])

def apu_prod(i, j):
	if (i < 0):
		return 0
	if (j < 0):
		return 0
	if (i > 0):
		return apu_prod(i - 1, j)*apu_fs(0)
	elif (j > 0):
		return apu_prod(i, j - 1)*apu_fs(1)
	return 1.0

def get_apu_coeff(l1):
	ap = 1.0
	ap2 = l1
	ap *= apu_prod(ap2, 1 - ap2)
	return ap

def get_apu_coeff2(m1, m2, i, j):
	ap = 1.0
	if (m1 == 1):
		ap *= (x - lambs[1 - i])
	if (m2 == 1):
		ap *= (x - lambs[1 - j])
	return ap

def make_coeffs():
	global Coeffs
	for i in xrange(2):
		for j in xrange(2):
			Coeffs[i][j] = 0
	for l1 in xrange(2):
		for i in xrange(2):
			for j in xrange(2):
				for m1 in xrange(2):
					for m2 in xrange(2):
						apu = 1.0
						apu *= B[i][l1]*B[j][l1]
						apu *= cs[i]*cs[j]
						apu *= get_apu_coeff(l1)
						apu *= get_apu_coeff2(m1, m2, i, j)
						Coeffs[m1][m2] += apu

def set_coeffs():
	global lambs, x, B, cs
	lambs[0] = 0
	lambs[1] = 1
	x = uniform(0, 1)
	cs = randn(2)
	B = get_B()
	if 0:
		print "B:"
		print B
		print "cs:"
		print cs

def make_plot():
	set_coeffs()
	xs = linspace(0, 1, 101)
	ys1, ys2, ys3, ys4 = [], [], [], []
	min_det = 0
	global x
	for xx in xs:
		x = xx
		make_coeffs()
		y1 = Coeffs[0][0]
		y2 = Coeffs[0][1]
		y3 = Coeffs[1][1]
		y4 = y1*y3-y2**2
		ys1.append(y1)
		ys2.append(y2)
		ys3.append(y3)
		ys4.append(y4)
		min_det = min(min_det, y4)
	plt.plot(xs, ys1, label='X')
	plt.plot(xs, ys2, label='Y')
	plt.plot(xs, ys3, label='Z')
	plt.plot(xs, ys4, label='XZ - Y^2')
	plt.plot(xs, [0]*101)
	plt.legend()
	print min_det
	plt.show()

def get_min():
	set_coeffs()
	xs = linspace(0, 1, 101)
	min_det = 0
	global x
	for xx in xs:
		x = xx
		make_coeffs()
		y1 = Coeffs[0][0]
		y2 = Coeffs[0][1]
		y3 = Coeffs[1][1]
		y4 = y1*y3-y2**2
		min_det = min(min_det, y4)
	return min_det

def check_mins():
	cap = 100
	minmin = 0
	for i in xrange(cap):
		minmin = min(minmin, get_min())
	return minmin

print check_mins()








