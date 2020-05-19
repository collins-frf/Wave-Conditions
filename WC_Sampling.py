import lhsmdu
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy

l = lhsmdu.sample(3,10) # Latin Hypercube Sampling of two variables, and 10 samples each.
k = lhsmdu.createRandomStandardUniformMatrix(3,10) # Monte Carlo Sampling
pdir=numpy.matrix('55 55 115 55 115 55 115 80 70 60 100 84 63 108 73 94 85 63 98 112')
H=numpy.matrix('1.2 1.2 1.2 2.3 2.3 2.3 2.3 .7 1.2 1.7 1.4 1.1 2.1 .9 2.0 1.5 1.9 1.6 2.1 1.8')
F=numpy.matrix('.091 .125 .125 .167 .167 .091 .091 .091 .091 .140 .11 .131 .152 .116 .097 .149 .162 .144 .138 .169')
k[0]=60*k[0]+55
l[0]=60*l[0]+55
k[1]=2.5*k[1]
l[1]=2.5*l[1]
k[2]=.1*k[2]+.085
l[2]=.1*l[2]+.085

fig = plt.figure()
ax = fig.add_subplot(131, projection='3d')
ax.scatter([k[0]], [k[1]], [k[2]], color="b", label="MC")
ax.scatter([l[0]], [l[1]], [l[2]], color="r", label="LHS")
ax.scatter([pdir], [H], numpy.transpose(F), color="g", label="OURS")
ax.set_xlabel('Peak Dir')
ax.set_ylabel('Wave Height')
ax.set_zlabel('Frequency')
ax.legend()
ax = fig.add_subplot(132)
ax.scatter([pdir], [H], color="g", label="OURS")
ax.set_xlabel('Peak Dir')
ax.set_ylabel('Wave Height')
plt.grid()
ax = fig.add_subplot(133)
ax.scatter([pdir], [F], color="g", label="OURS")
ax.set_xlabel('Peak Dir')
ax.set_ylabel('Frequency')
plt.grid()
plt.show()