RR.<t, x, a, b, c, d, X11, X12, X22, X13, X23, X33, X14, X24, X34, X44, d_ab, d_ac, d_ad, d_bc, d_bd, d_cd> = QQ[]
K = RR.fraction_field()

X = [[X11, X12, X13, X14], [X12, X22, X23, X24], [X13, X23, X33, X34], [X14, X24, X34, X44]]

ass = [a, b, c, d]
ds = [[0, d_ab, d_ac, d_ad], [-d_ab, 0, d_bc, d_bd], [-d_ac, -d_bc, 0, d_cd], [-d_ad, -d_bd, -d_cd, 0]]

k = 3
n = 4
p = 0

vals = [[0 for __ in xrange(k - 1)] for _ in xrange(n - 1)]

def list_index(l):
	ind = 0
	for i in xrange(len(l)):
		if ((l[i - 1] <= p) and (l[i] > p)):
			ind += 1
	return ind

def set_divided_difference(l, coeff):
	if (l[0] == l[-1]):
		return
	i_count = 0
	while (l[i_count] == l[0]):
		i_count += 1
	vals[l[0] - 1][k - 1 - i_count] += coeff


def open_divided_difference(l, coeff = 1):
	l = sorted(l)
	if (l[0] + 1 >= l[-1]):
		set_divided_difference(l, coeff)
	else:
		i = l[0]
		j = l[0] + 1
		ii = l[-1]
		coeff1 = ds[j - 1][i - 1]/ds[ii - 1][i - 1]
		coeff2 = ds[ii - 1][j - 1]/ds[ii - 1][i - 1]
		open_divided_difference(l[:-1] + [j], coeff1*coeff)
		open_divided_difference(l[1:] + [j], coeff2*coeff)

def pprint(i, j):
	print "[",
	for ii in xrange(k - j - 1):
		print ass[i], ", ",
	for ii in xrange(j):
		print ass[i + 1], ", ",
	print ass[i + 1], "]"

def print_table():
	for i in xrange(n - 1):
		for j in xrange(k - 1):
			pprint(i, j)
			print vals[i][j]

def fact(n):
	if (n == 0):
		return 1
	return n*fact(n - 1)

def binom(n, k):
	if (n < k):
		return 0
	if (k < 0):
		return 0
	return fact(n)/(fact(k)*fact(n - k))

def sum_table(l, d):
	if (d == 0):
		return l
	ans = []
	ans.append(l[0])
	for i in xrange(len(l) - 1):
		ans.append(l[i] + l[i + 1])
	ans.append(l[-1])
	return sum_table(ans, d - 1)

def get_terms():
	return [vals[0][j]*binom(k - 2, j) for j in xrange(k - 1)]

def print_expanded_table_first(d = 1):
	l = sum_table(get_terms(), d)
	for i in xrange(len(l)):
		print i
		print l[i]

def get_polynomial_i(i):
	ans = 0
	for j in xrange(i):
		z = (t + ds[0][j + 1])/ds[j][j + 1]
		w = 1 - z
		for jj in xrange(k - 1):
			ans += z**(k - 2 - jj)*w**(jj)*binom(k - 2, jj)*vals[j][jj]
	return ans


def go(l = [], coeff = 1):
	if (len(l) == k):
		open_divided_difference(l, coeff*X[l[0] - 1][l[-1] - 1])
		return
	for i in xrange(1, n + 1):
		coeff2 = 1
		if (len(l) > 0):
			coeff2 *= X[l[-1] - 1][i - 1]
		go([x for x in l] + [i], coeff2*coeff)

def set_all_terms():
	go()

polys = [0 for _ in xrange(k/2 + 1)]

def basic_div(i, m):
	if (m == k):
		return 0
	if (i <= p):
		return 0
	return (ass[i] - x)**(k - 1 - m)*binom(k - 2, m - 1)

def get_X_coeff(l):
	coeff = 1
	for i in xrange(k):
		coeff *= X[i - 1][i]
	return coeff

def get_divided_difference_polynomial(l):
	l = sorted(l)
	if (l[0] < l[-1]):
		return (get_divided_difference_polynomial(l[:-1]) - get_divided_difference_polynomial(l[1:]))/(ass[l[0]] - ass[l[-1]]) #ds[l[0]][l[-1]]
	return basic_div(l[0], len(l))

def set_div(l):
	ind = list_index(l)
	pol = get_divided_difference_polynomial(l)
	coeff = get_X_coeff(l)
	global polys
	polys[ind] += coeff*pol
	return

def go_poly(l=[]):
	if (len(l) == k):
		set_div(l)
		return
	for i in xrange(n):
		l2 = [xx for xx in l]
		l2.append(i)
		go_poly(l2)
	return

go_poly()
for i in xrange(k/2 + 1):
	print i, polys[i]
