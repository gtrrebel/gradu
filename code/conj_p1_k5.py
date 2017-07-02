from numpy.random import randn, randint
import matplotlib.pyplot as plt
import numpy as np
eps = 1e-4

# x is always 0
lambs = []
B = [[]]
res_list = [[]]
equal = []
form_Bs = [[]]
form_Bg = [[]]

def mat_prods(l):
	if (len(l) == 1):
		return l[0]
	if (len(l) == 0):
		return 0
	return np.dot(mat_prods(l[:-1]), l[-1])

def set_form(Bs, Bg):
	global form_Bs, form_Bg
	form_Bs = Bs
	form_Bg = Bg

def gf(*arg):
	return mat_prods([form_Bs] + list(arg) + [form_Bg])[0][0]

def set_lambs(n, proj=True):
	global lambs
	if (proj):
		ind = randint(0, n - 2)
		lambs = [0 for _ in xrange(n - 1)]
		lambs[ind] = abs(randn(1)[0])
	else:
		lambs = randn(n - 1)
		lambs = [abs(lamb) for lamb in lambs]
		lambs = sorted(lambs)

def set_B(n, k = 5, split=True, proj=True):
	global B
	B = np.zeros((n, n))
	if proj:
		if split:
			v = randn(n - 1)
			for _ in xrange(n - 1):
				for __ in xrange(n - 1):
					B[_][__] = v[_]*v[__]
			v2 = randn(n - 1)
			for _ in xrange(n - 1):
				B[n - 1][_] = v2[_]
				B[_][n - 1] = v2[_]
			B[n - 1][n - 1] = randn(1)[0]**2
		else:
			v = randn(n)
			for _ in xrange(n):
				for __ in xrange(n):
					B[_][__] = v[_]*v[__]
	else:
		if split:
			for _ in xrange(n - 1):
				for __ in xrange(n - 1):
					B[_][__] = randn(1)[0]
			B = B + B.transpose()
			B = np.dot(B,np.transpose(B))
			v2 = randn(n - 1)
			for _ in xrange(n - 1):
				B[n - 1][_] = v2[_]
				B[_][n - 1] = v2[_]
			B[n - 1][n - 1] = randn(1)[0]**2
		else:
			for _ in xrange(n):
				for __ in xrange(n):
					B[_][__] = randn(1)[0]
			B = B + B.transpose()
			B = np.dot(B,np.transpose(B))

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
		return (ans1 - ans2)/(lambs[l[0]] - lambs[l[-1]])
	else:
		return get_basic_divided(lambs[l[0]], len(l), k)

def term_00():
	res = 0
	for j1 in xrange(n - 1):
		for j2 in xrange(n - 1):
			for j3 in xrange(n - 1):
				for j4 in xrange(n - 1):
					l = [j1, j2, j3, j4]
					B_coeff = get_B_coeff(l + [n - 1 for _ in xrange(1)])
					div_dif = get_divided_difference(l, 5)
					res += B_coeff*div_dif
	return res

def term_01():
	res = 0
	for j1 in xrange(n - 1):
		for j2 in xrange(n - 1):
			for j3 in xrange(n - 1):
				l = [j1, j2, j3]
				B_coeff = get_B_coeff(l + [n - 1 for _ in xrange(2)])
				div_dif = get_divided_difference(l, 5)
				res += B_coeff*div_dif
	return res

def term_02():
	res = 0
	for j1 in xrange(n - 1):
		for j2 in xrange(n - 1):
			l = [j1, j2]
			B_coeff = get_B_coeff(l + [n - 1 for _ in xrange(3)])
			div_dif = get_divided_difference(l, 5)
			res += B_coeff*div_dif
	return res

def term_03():
	res = 0
	for j1 in xrange(n - 1):
		l = [j1]
		B_coeff = get_B_coeff(l + [n - 1 for _ in xrange(4)])
		div_dif = get_divided_difference(l, 5)
		res += B_coeff*div_dif
	return res

def term_10():
	res = 0
	for j1 in xrange(n - 1):
		for j2 in xrange(n - 1):
			for j3 in xrange(n - 1):
				B_coeff = get_B_coeff([j1, j2, n - 1, j3, n - 1])
				div_dif = get_divided_difference([j1, j2, j3], 5)
				res += B_coeff*div_dif
	return res

def term_11():
	res = 0
	for j1 in xrange(n - 1):
		for j2 in xrange(n - 1):
			B_coeff = get_B_coeff([j1, n - 1, j2, n - 1, n - 1])
			div_dif = get_divided_difference([j1, j2], 5)
			res += B_coeff*div_dif
	return res

def get_Bp():
	return B[:-1,:-1]

def get_Bs():
	return B[-1,:-1].reshape(1, n - 1)

def get_Bg():
	return B[:-1,-1].reshape(n - 1, 1)

def get_Bm():
	return B[-1,-1]

def get_D():
	return np.diag(lambs)

def get_terms():
	res = 0
	t00 = term_00()
	t01 = term_01()
	t02 = term_02()
	t03 = term_03()
	t10 = term_10()
	t11 = term_11()
	#print t00, t01, t02, t03, t10, t11
	res = t00 + t01 + t02 + t03 + t10 + t11
	return res

def get_terms2():
	Bp = get_Bp()
	Bs = get_Bs()
	Bg = get_Bg()
	Bm = get_Bm()
	set_form(Bs, Bg)
	D = get_D()
	t00 = gf(Bp, Bp, Bp)
	t01 = Bm*(gf(D, Bp, Bp) + gf(Bp, D, Bp) + gf(Bp, Bp, D))
	t02 = Bm**2*(gf(Bp, D, D) + gf(D, Bp, D) + gf(D, D, Bp))
	t03 = Bm**3*gf(D, D, D)
	t10 = gf(D, Bp)*gf() + gf(Bp, D)*gf() + gf(Bp)*gf(D)
	t11 = Bm*(gf(D, D)*gf() + gf(D)*gf(D) + gf()*gf(D, D))
	res = t00 + t01 + t02 + t03 + t10 + t11
	return res

def get_terms3():
	Bp = get_Bp()
	Bs = get_Bs()
	Bg = get_Bg()
	Bm = get_Bm()
	set_form(Bs, Bg)
	D = get_D()
	Bps = Bp - np.dot(Bg, Bs)/Bm
	t00 = 0
	t00 += gf(Bps, Bps, Bps)
	t00 += Bm*(gf(D, Bps, Bps) + gf(Bps, Bps, D))
	t00 += Bm**2*gf(D, Bps, D)
	t01 = 0
	t01 += Bm*gf(Bps, D, Bps)
	t01 += Bm**2*(gf(Bps, D, D) + gf(D, D, Bps))
	t01 += Bm**3*gf(D, D, D)
	t10 = (2*gf(Bps, Bps)/Bm + 3*gf(Bps, D) + 3*gf(D, Bps) + 4*Bm*gf(D, D))*gf()
	t11 = gf(Bps)*gf(Bps)/Bm + 3*gf(Bps)*gf(D) + 2*gf(D)*gf(D)*Bm
	t20 = (6*gf(D)/Bm + 3*gf(Bps)/Bm**2)*gf()*gf()
	t30 = gf()*gf()*gf()*gf()/Bm**3
	res = 0
	res += t00
	res += t01
	res += t10
	res += t11
	res += t20
	res += t30
	print t00, t01, t10, t11, t20, t30
	return res

def get_terms4():
	Bps = get_Bp()
	Bs = get_Bs()
	Bg = get_Bg()
	Bm = get_Bm()
	set_form(Bs, Bg)
	D = get_D()
	t00 = 0
	t00 += gf(Bps, Bps, Bps)
	t00 += Bm*(gf(D, Bps, Bps) + gf(Bps, Bps, D))
	t00 += Bm**2*gf(D, Bps, D)
	t01 = 0
	t01 += Bm*gf(Bps, D, Bps)
	t01 += Bm**2*(gf(Bps, D, D) + gf(D, D, Bps))
	t01 += Bm**3*gf(D, D, D)
	t10 = (2*gf(Bps, Bps)/Bm + 3*gf(Bps, D) + 3*gf(D, Bps) + 4*Bm*gf(D, D))*gf()
	t11 = gf(Bps)*gf(Bps)/Bm + 3*gf(Bps)*gf(D) + 2*gf(D)*gf(D)*Bm
	t20 = (6*gf(D)/Bm + 3*gf(Bps)/Bm**2)*gf()*gf()
	t30 = gf()*gf()*gf()*gf()/Bm**3
	res = 0
	res += t00
	res += t01
	res += t10 #Only bad term
	res += t11
	res += t20 #This with t10 fails with eps = 1e-4
	res += t30
	return res

def get_terms5():
	#FOUND FAILURE
	Bps = get_Bp()
	Bs = get_Bs()
	Bg = get_Bg()
	set_form(Bs, Bg)
	D = get_D()
	return 2*gf(Bps, Bps) + 3*gf(Bps, D) + 3*gf(D, Bps) + 4*gf(D, D)

def try_conj(n, cap = 10000000):
	for _ in xrange(cap):
		set_lambs(n)
		set_B(n)
		res4 = get_terms4()
		if (res4 < -eps):
			print _
			print "FAIL"
			print res4
			exit(0)


n = 3
try_conj(n)








