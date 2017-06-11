from numpy import *

NE = 10  # 100
dX = 0.1  # 0.05
AA = 0.1
EE = 1000
BC = [[1, 0.]]
LOAD = [[7, 1, 2.], [10, 1, -3.5]]  # [[70,2.],[100,-3.5]]

fileName = "input3.txt"
ff = open(fileName, 'w')
NN = NE + 1

ff.write('*node\n')
for i in xrange(NN):
    ff.write('%4i,     %8.3f,  0.,  0.\n' % (i + 1, i * dX))
ff.write('*element, eltype=T1D2\n')
for i in xrange(NE):
    ff.write('%4i,   %4i,%4i,   %8.3f\n' % ((i + 1), i + 1, i + 2, AA))
ff.write('*material\n   %8.2e\n' % (EE))
ff.write('*loading\n')
for i in LOAD:
    ff.write('%4i, 1, %8.3f\n' % (i[0], i[1]))
ff.write('*boundary\n')
for i in BC:
    ff.write('%4i, 1, %8.3f\n' % (i[0], i[1]))

ff.close()
print 'finish'