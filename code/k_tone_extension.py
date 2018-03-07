from numpy.random import randn

def get_numbers(k):
	l = randn(2*k + 1)
	l = sorted(l)
	xs = l[:k]
	ys = l[k + 1:]
	x = l[k]
	return xs, ys, x

def get_divs(k):
	l = randn(k)
	l = [abs(x) for x in l]
	return l

def get_val(k, xs, ys, divs, i, j):
	if ((i < 0) and (j < 0)):
		return 0
	if ((j - i) == k):
		return divs[j]
	diff = ys[j]
	if (i - 1 >= 0):
		diff -= ys[i - 1]
	else:
		diff -= xs[i - 1 + k]
	return diff*get_val(k, xs, ys, divs, i - 1, j) + get_val(k, xs, ys, divs, i - 1, j - 1)

def get_zs(k, xs, ys, divs):
	return [get_val(k, xs, ys, divs, i, i) for i in xrange(k)]

def get_interpolating(k, xs, ys, divs, x, i):
	i -= k
	coeff = 1
	su = 0
	for j in xrange(k):
		su += get_val(k, xs, ys, divs, i, i + j)*coeff
		if (i + j < 0):
			coeff *= (x - xs[i + j + k])
		else:
			coeff *= (x - ys[i + j])
	return su


def get_boundary_values(k, xs, ys, divs, x):
	return [get_interpolating(k, xs, ys, divs, x, i) for i in xrange(k + 1)]

def do_random_check(k):
	xs, ys, x = get_numbers(k)
	divs = get_divs(k)
	bvs = get_boundary_values(k, xs, ys, divs, x)
	evens = bvs[0::2]
	odds = bvs[1::2]
	mi = max(evens)
	ma = min(odds)
	if (mi > ma):
		print "FAIL"
		print xs, ys, x
		print divs
		print bvs
		print ma - mi

def brute_test(k, N = 100):
	for _ in xrange(N):
		do_random_check(k)

k = 4
N = 100

brute_test(k, N)





