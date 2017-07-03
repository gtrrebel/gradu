from numpy.random import randn, randint
import numpy as np
eps = 1e-4

def get_random_pos(n):
	A = np.zeros((n, n))
	for i in xrange(n):
		for j in xrange(n):
			A[i][j] = randn()
	return np.dot(A, np.transpose(A))

def is_pos(x):
	return x + eps > 0

def get_Id(n):
	return np.eye(n)

def get_random_product(l, n):
	A = get_random_pos(n)
	B = get_random_pos(n)
	C = get_Id(n)
	for x in l:
		if (x == 0):
			C = np.dot(C, A)
		else:
			C = np.dot(C, B)
	return C

def test_random_products(l, n, cap = 10000000):
	for _ in xrange(cap):
		C = get_random_product(l, n)
		tr = np.trace(C)
		if not is_pos(tr):
			print "FAIL"
			print "tr: ", tr
			print "l: ", l
			print "C: "
			print C
			return False
	return True


n = 3
l = [0, 1, 0, 1, 1, 0]

print test_random_products(l, n)