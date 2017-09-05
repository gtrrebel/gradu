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
	B = np.dot(B,np.transpose(B))
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

def get_B3():
	n = 2
	B = np.zeros((n, n))
	B[0][0] = 4.0
	B[0][1] = 3.0
	B[1][0] = 3.0
	B[1][1] = 4.0
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

def get_apu_coeff(l1, l2):
	ap = 1.0
	if (l1 == l2):
		ap *= 2.0
	else:
		ap *= 3.0
	ap2 = l1 + l2
	ap *= apu_prod(2 - ap2, ap2)
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
		for l2 in xrange(2):
			for i in xrange(2):
				for j in xrange(2):
					for m1 in xrange(2):
						for m2 in xrange(2):
							apu = 1.0
							apu *= B[l1][l2]
							apu *= B[i][l1]*B[j][l2]
							apu *= cs[i]*cs[j]
							apu *= get_apu_coeff(l1, l2)
							apu *= get_apu_coeff2(m1, m2, i, j)
							Coeffs[m1][m2] += apu

def get_LHS(z):
	res = 0
	for i in xrange(2):
		for j in xrange(2):
			for i1 in xrange(2):
				for i2 in xrange(2):
					ap1 = B[i][i1]*B[i1][i2]*B[i2][j]*cs[i]*cs[j]
					ap2 = ((z - lambs[i1])*(z - lambs[i2])*(z - lambs[i])*(z - lambs[j]))
					res += ap1/ap2
	return res

def set_coeffs():
	global lambs, x, B, cs
	lambs[0] = 0
	lambs[1] = 1
	x = uniform(0, 1)
	cs = [1.0, -1.0]
	B = get_B3()
	if 0:
		print "B:"
		print B
		print "B_det"
		print B[0][0]*B[1][1] - B[0][1]**2
		print "cs:"
		print cs

def get_extra_term(x, c = 0.0, d = 0.0, f = 0.0):
	v1 = [0.5*(5*x**2 - 5*x + 1)*(x**2 - x), (x**2 - x)**2*(2*x - 1), 0]
	v1 = [c*v_ for v_ in v1]
	v2 = [0, (x**2 - x)*(2*x - 1)*0.2, 2*(x**2 - x)**2]
	v2 = [d*v_ for v_ in v2]
	v3 = [0, (2*x - 1)*(x**2 - x)**2*0.3, 2*(x**2 - x)**3]
	v3 = [f*v_ for v_ in v3]
	v = [x + y + z for (x, y, z) in zip(v1, v2, v3)]
	return v[0], v[1], v[2]

def get_help(y1, y2, y3, x, z1 = -2.0, z2 = -1.5, ch = 0.5):
	if (x <= ch):
		return (z1**2)*(y1/(z1 - x)**4 + 2*y2/(z1 - x)**5 + y3/(z1 - x)**6)
	else:
		return (z2**2)*(y1/(z2 - x)**4 + 2*y2/(z2 - x)**5 + y3/(z2 - x)**6)

# For B = [[1, 1], [1, 2]]

# X is given by x |-> 4*x^2
#plt.plot(xs, [4*xi**2 for xi in xs], label='X_test')
# Y is given by x |-> 4*x^3 + x^2 + 3*x
#plt.plot(xs, [3*xi  + xi**2 + 4*xi**3 for xi in xs], label='Y_test')
# Z is given by x |-> 4*x^4 + 2*x^3 + 6*x^2 + 2*x + 2
#plt.plot(xs, [2 + 2*xi + 6*xi**2 + 2*xi**3 + 4*xi**4 for xi in xs], label='Z_test')
#plt.plot(xs, [yi/(xi**2*((1.0-xi)**2)) for (xi, yi) in zip(xs, ys4)], label='Det_test')

# For B = [[4, 3], [3, 4]]

# X is given by x |-> 34*x^2 - 34*x + 8
#plt.plot(xs, [34*xi**2 - 34*xi + 8 for xi in xs], label='X_test')
# Y is given by x |-> 34*x^3 - 51*x^2 + 81*x - 32
#plt.plot(xs, [34*xi**3 - 51*xi**2 + 81*xi - 32 for xi in xs], label='Y_test')
# Z is given by x |-> 34*x^4 - 68*x^3 + 138*x^2 - 104*x + 128
#plt.plot(xs, [34*xi**4 - 68*xi**3 + 138*xi**2 - 104*xi + 128 for xi in xs], label='Z_test')

# c should be on interval -2...-8 for x to be non-negative


def make_plot():
	set_coeffs()
	cap = 10001
	xs = linspace(0, 1, cap)
	ys1, ys2, ys3, ys4, ys5 = [], [], [], [], []
	min_det = 0
	global x
	z = -3.141
	su = 0
	for xx in xs:
		x = xx
		make_coeffs()
		x1, x2, x3 = get_extra_term(x, c= 10.0)
		apu_p = apu_prod(1, 1)
		y1 = Coeffs[0][0]*apu_p + x1
		y2 = Coeffs[0][1]*apu_p + x2
		y3 = Coeffs[1][1]*apu_p + x3
		y4 = y1*y3-y2**2
		app = (y1/(z - x)**4 + 2*y2/(z - x)**5 + y3/(z - x)**6)
		y5 = app
		su += y5/cap
		ys1.append(y1)
		ys2.append(y2)
		ys3.append(y3)
		ys4.append(y4)
		ys5.append(y5)
		min_det = min(min_det, y4)
	plt.plot(xs, ys1, label='X')
	#plt.plot(xs, ys2, label='Y')
	#plt.plot(xs, ys3, label='Z')
	#plt.plot(xs, ys4, label='XZ - Y^2') # seems to be always given by x |-> c*x^2*(1 - x)^2 for some c depending on B and cs, probably
	#plt.plot(xs, ys5, label='Lol')
	plt.plot(xs, [0]*cap)
	plt.legend()
	lol2 = get_LHS(z)
	print "min_det:", min_det
	print su, lol2
	print su/lol2
	plt.show()

def ce_3(x):
	return 1/((2 - x)*(2 + x)**3)

def ce_4(x):
	return (x - 1)/((2 - x)**2*(x + 2)**4)

def ce_5(x):
	return (5*x**2 - 10*x + 8)/(5*(2 - x)**3*(x + 2)**5)

def get_help_lol(y1, y2, y3, x):
	return (y1*ce_3(x) + 2*y2*ce_4(x) + y3*ce_5(x))

def test_counterexample():
	set_coeffs()
	cap = 1001
	xs = linspace(0, 1, cap)
	ys1, ys2, ys3, ys4, ys5 = [], [], [], [], []
	min_det = 0
	global x
	z = -5.0
	su = 0
	for xx in xs:
		x = xx
		make_coeffs()
		y1 = Coeffs[0][0]
		y2 = Coeffs[0][1]
		y3 = Coeffs[1][1]
		y4 = y1*y3-y2**2
		app = get_help_lol(y1, y2, y3, x)
		#app = (y1/(z - x)**4 + 2*y2/(z - x)**5 + y3/(z - x)**6)
		y5 = app*apu_prod(1, 1)
		su += y5/cap
		ys1.append(y1)
		ys2.append(y2)
		ys3.append(y3)
		ys4.append(y4)
		ys5.append(y5)
		min_det = min(min_det, y4)
	plt.plot(xs, ys1, label='X')
	#plt.plot(xs, ys2, label='Y')
	#plt.plot(xs, ys3, label='Z')
	plt.plot(xs, ys4, label='XZ - Y^2')
	plt.plot(xs, ys5, label='Lol')
	plt.plot(xs, [0]*cap)
	plt.legend()
	lol2 = get_LHS(z)
	print "min_det:", min_det
	print su, lol2
	print su/lol2
	plt.show()

def make_all_plots():
	bnd_c = 20.0
	bnd_d = 5.0
	bnd_f = 30.0
	set_coeffs()
	z = -3.141
	lol2 = get_LHS(z)
	for c in linspace(-bnd_c, bnd_c, 7):
		for d in linspace(-8.0, -2.0, 7):
			for f in linspace(-bnd_f, bnd_f, 7):
				set_coeffs()
				cap = 100001
				xs = linspace(0, 1, cap)
				ys1, ys2, ys3, ys4, ys5 = [], [], [], [], []
				min_det = 0
				global x
				z = -3.141
				su = 0
				for xx in xs:
					x = xx
					make_coeffs()
					x1, x2, x3 = get_extra_term(x, c, d, f)
					apu_p = apu_prod(1, 1)
					y1 = Coeffs[0][0]*apu_p + x1
					y2 = Coeffs[0][1]*apu_p + x2
					y3 = Coeffs[1][1]*apu_p + x3
					y4 = y1*y3-y2**2
					app = (y1/(z - x)**4 + 2*y2/(z - x)**5 + y3/(z - x)**6)
					y5 = app
					su += y5/cap
					ys1.append(y1)
					ys2.append(y2)
					ys3.append(y3)
					ys4.append(y4)
					ys5.append(y5)
					min_det = min(min_det, y4)
				plt.plot(xs, ys1, label='X')
				#plt.plot(xs, ys2, label='Y')
				#plt.plot(xs, ys3, label='Z')
				#plt.plot(xs, ys4, label='XZ - Y^2') # seems to be always given by x |-> c*x^2*(1 - x)^2 for some c depending on B and cs, probably
				#plt.plot(xs, ys5, label='Lol')
				plt.plot(xs, [0]*cap)
				print su/lol2
	plt.plot(xs, [0]*cap, color="black")
	plt.show()



def test1():
	set_coeffs()
	xs = linspace(0, 1, 101)
	ys1, ys2, ys3, ys4 = [], [], [], []
	min_det = 0
	global x
	x = 0.5
	make_coeffs()
	y1 = Coeffs[0][0]
	y2 = Coeffs[0][1]
	y3 = Coeffs[1][1]
	y4 = y1*y3-y2**2
	print y1, y2, y3
	print y4
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

print make_plot()



