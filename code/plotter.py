import matplotlib.pyplot as plt
import sys

inp = sys.argv
name = inp[1]

f = open(name, 'r')

xs = []
ys = []

params = f.readline().split()

n = int(params[0])
cap = int(params[1])

ass = [float(xx) for xx in f.readline().split()]
B = [[float(xx) for xx in f.readline().split()] for i in xrange(n)]

max_val = 0

for _ in xrange(cap):
	line = f.readline()
	print line
	zs = [float(xx) for xx in line.split()]
	xs.append(zs[0])
	ys.append(zs[1])
	max_val = max(max_val, ys[-1])

f.close()

plt.plot(xs, ys)
for i in xrange(n):
	plt.plot([ass[i], ass[i]], [0, max_val])

plt.show()
