import matplotlib.pyplot as plt
import numpy as np

eps = 1e-9

default_t_min = -1.0
default_t_max = 1.0
default_N = 1000
default_k = 3

def get_random_ts(k = default_k, t_min = default_t_min, t_max = default_t_max):
	return sorted(np.random.uniform(t_min, t_max, k + 1))

def get_spline_value(t, xs):
	if (t < xs[0] + eps):
		return 0
	if (t > xs[-1] - eps):
		return 0
	if (len(xs) == 2):
		return 1.0/(xs[1] - xs[0])
	val1 = get_spline_value(t, xs[:-1])
	val2 = get_spline_value(t, xs[1:])
	return ((t - xs[0])/(xs[-1] - xs[0])*val1 + (t - xs[-1])/(xs[0] - xs[-1])*val2)/(len(xs) - 2)

def add_y_to_xs(xs, y):
	xs.append(y)
	xs = sorted(xs)
	return xs

def get_spline(ts, xs):
	return [get_spline_value(t, xs) for t in ts]

def plot_spline(xs, t_min = default_t_min, t_max = default_t_max, N = default_N):
	ts = np.linspace(t_min, t_max, N)
	ys = get_spline(ts, xs)
	plt.plot(ts, ys)
	plt.show()

def recursive_spline_split(xs, coeff, a, b, c, d, level, spline_list):
	if (len(spline_list) <= level):
		spline_list.append([])
	if (xs[-1] < b) or (xs[0] > c):
		spline_list[level].append((xs, coeff, True))
		return
	else:
		spline_list[level].append((xs, coeff, False))
	new_x = np.random.uniform(c, b)
	xs_1 = add_y_to_xs(xs[:-1], new_x)
	co_1 = coeff*(new_x - xs[0])/(xs[-1] - xs[0])
	xs_2 = add_y_to_xs(xs[1:], new_x)
	co_2 = coeff*(xs[-1] - new_x)/(xs[-1] - xs[0])
	recursive_spline_split(xs_1, co_1, a, b, c, d, level + 1, spline_list)
	recursive_spline_split(xs_2, co_2, a, b, c, d, level + 1, spline_list)

def get_split_maximum(xs):
	t_min = xs[0]
	t_max = xs[-1]
	N = default_N
	ts = np.linspace(t_min, t_max, N)
	return max(get_spline(ts, xs))

def plot_beauty_split(spline_list, a = default_t_min, d = default_t_max, N = default_N, do_copies = False):
	xs = spline_list[0][0][0]
	ts = np.linspace(a, d, N)
	ma = 1.5*get_split_maximum(xs)
	le = len(spline_list)
	for l in xrange(le):
		for spline in spline_list[l]:
			gap = l + 1
			if ((do_copies) and spline[2]):
				gap = le
			for ll in xrange(l, gap):
				ys = get_spline(ts, spline[0])
				ys = [y*spline[1] - ll*ma for y in ys]
				plt.plot(ts, ys)
	plt.show()


def beauty_plot_split(xs, a, b, c, d):
	spline_list = []
	recursive_spline_split(xs, 1.0, a, b, c, d, 0, spline_list)
	plot_beauty_split(spline_list, a, d)

xs = get_random_ts(k = 4)
a = -1.0
b = 0.25
c = -0.25
d = 1.0

beauty_plot_split(xs, a, b, c, d)
