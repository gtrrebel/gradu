from numpy.random import randn, randint, uniform
import matplotlib.pyplot as plt
import numpy as np
eps = 1e-6
brute_check_gap = 100

def is_same(x, y):
	return abs(x - y) < eps

def is_at_most(x, y):
	return x < y + eps

def get_random_tuple(n):
	return randn(n)

def normalize_tuple(l):
	l = sorted(l)
	su = sum(l)*1.0
	n = len(l)*1.0
	l = [x - su/n for x in l]
	su2 = sum(x**2 for x in l)
	ap = su2**(0.5)
	l = [x/ap for x in l]
	return l

def get_random_normalized_tuple(n):
	return normalize_tuple(get_random_tuple(n))

def get_f(x):
	if (x < 0):
		return 0
	else:
		return x**2

def get_val(x, l):
	return sum(get_f(y - x) for y in l)

def brute_check_majorization(l1, l2):
	su11 = sum(l1)
	su12 = sum(x**2 for x in l1)
	su21 = sum(l2)
	su22 = sum(x**2 for x in l2)
	if not is_same(su11, su21):
		print "fail1"
		return False
	if not is_same(su12, su22):
		print "fail2"
		return False
	min_val = min(min(l1), min(l2))
	max_val = max(max(l1), max(l2))
	for _ in xrange(brute_check_gap):
		x = uniform(min_val, max_val)
		v1 = get_val(x, l1)
		v2 = get_val(x, l2)
		if not is_at_most(v1, v2):
			#print _
			#print "lol"
			#print v1, v2
			#print l1
			#print l2
			return False
	return True

def get_cut_point(cur_speed, new_speed, le):
	t = abs(cur_speed)
	s = abs(new_speed)
	return le*t/(t + s)

def get_total_area(cur_speed, new_speed, le):
	return (cur_speed + new_speed)*0.5*le

def get_lowest_total(cur_speed, new_speed, le):
	tot = get_total_area(cur_speed, new_speed, le)
	if is_at_most(cur_speed, 0) and is_at_most(new_speed, 0):
		return tot, tot
	elif is_at_most(0, cur_speed) and is_at_most(0, new_speed):
		return 0, tot
	elif is_at_most(0, cur_speed) and is_at_most(new_speed, 0):
		return min(0, tot), tot
	else:
		lep = get_cut_point(cur_speed, new_speed, le)
		return get_total_area(cur_speed, 0, lep), tot

def clever_check_majorization(l1, l2):
	n = len(l1)
	started = False
	prev_x = 0
	cur_acc = 0
	cur_speed = 0
	cur_area = 0
	cur_x = 0
	new_speed = 0
	i1 = 0
	i2 = 0
	while (i1 + i2 < 2*n):
		if (i1 == n):
			sgn = -1
		elif i2 == n:
			sgn = 1
		elif (l1[i1] < l2[i2]):
			sgn = 1
		else:
			sgn = -1
		if (sgn == 1):
			cur_x = l1[i1]
			i1 += 1
		else:
			cur_x = l2[i2]
			i2 += 1
		if started:
			new_speed = cur_acc*(cur_x - prev_x) + cur_speed
			low, tot = get_lowest_total(cur_speed, new_speed, cur_x - prev_x)
			if not is_at_most(0, low + cur_area):
				return False
			cur_area += tot
		else:
			started = True
		cur_acc += sgn
		prev_x = cur_x
		cur_speed = new_speed
	if not is_same(cur_area, 0):
		return False
	if not is_same(cur_speed, 0):
		return False
	return True

def lol_check(l1, l2):
	As1 = [sum(l1[i]**2 for i in xrange(n) if i <= j) for j in xrange(n)]
	Bs1 = [sum(l2[i]**2 for i in xrange(n) if i <= j) for j in xrange(n)]
	for i in xrange(n):
		apu = As1[i] - Bs1[i]
		if (apu + eps < 0):
			print apu,
	print

def random_check(l1, l2):
	as1 = [sum(l1[i] for i in xrange(n) if i <= j) for j in xrange(n)]
	bs1 = [sum(l2[i] for i in xrange(n) if i <= j) for j in xrange(n)]
	As1 = [sum(l1[i]**2 for i in xrange(n) if i <= j) for j in xrange(n)]
	Bs1 = [sum(l2[i]**2 for i in xrange(n) if i <= j) for j in xrange(n)]
	if not is_same(as1[-1], bs1[-1]):
		return False
	if not is_same(As1[-1], Bs1[-1]):
		return False
	for i in xrange(1, n - 1, 2):
		apu = (Bs1[i] - As1[i - 1] - (bs1[i] - as1[i - 1])**2)
		if not is_at_most(apu, 0):
			return False
	for i in xrange(2, n - 1, 2):
		apu = (As1[i] - Bs1[i - 1] - (as1[i] - bs1[i - 1])**2)*(-1)
		if not is_at_most(apu, 0):
			return False
	return True

def get_order(l1, l2):
	l1p = [(x, 0) for x in l1]
	l2p = [(x, 1) for x in l2]
	l3 = l1p + l2p
	l3 = sorted(l3)
	return [a for (x, a) in l3]

def check_if_works(n = 3, cap = 10000):
	for _ in xrange(cap):
		l1 = get_random_normalized_tuple(n)
		l2 = get_random_normalized_tuple(n)
		ok1 = clever_check_majorization(l1, l2)
		if (ok1):
			lol_check(l1, l2)
		#ok2 = brute_check_majorization(l1, l2)
		#ok3 = random_check(l1, l2)
		if 0:
			print _
			print ok1, ok3
			print "ls:"
			print l1
			print l2
			print "FAIL"
			print get_order(l1, l2)
			return False
	return True

l1 = normalize_tuple([1, 2, 3])
l2 = normalize_tuple([1, 2, 4])
print l1
print l2

n = 7
print check_if_works(n)

