import numpy as np
import cmath
from cmath import exp
import matplotlib.pyplot as plt

def get_val(x, theta, n):
	return (x - exp(1j*theta*n))

def get_vals(x, theta, n):
	l = [get_val(x, theta, i) for i in xrange(1, n + 1)]
	for i in xrange(n - 1):
		l[i + 1] *= l[i]
	return l

def make_plot(x, theta, n):
	seq = get_vals(x, theta, n)
	if 0:
		seq = np.real(seq)
	else:
		seq = np.abs(seq)
	seq = sorted(seq)
	plt.plot(range(n), seq)
	plt.show()

def get_x_min_max(x, theta, n):
	seq = np.abs(get_vals(x, theta, n))
	return min(seq), max(seq)

def plot_min_max_curves(xs, theta, n):
	mins, maxs = [], []
	for x in xs:
		mi, ma = get_x_min_max(x, theta, n)
		mins.append(mi)
		maxs.append(ma*(1.0 - x))
	plt.plot(xs, mins)
	plt.plot(xs, maxs)
	plt.show()

xs = np.linspace(0, 0.9, 50)
theta = 1.0
n = 5000000

plot_min_max_curves(xs, theta, n)
	
