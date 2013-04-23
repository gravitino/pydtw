import libdtw as dtw
import math as m

# define query and subject as shifted cosine waves
query = dtw.TimeSeries([m.cos(i/128.0) * m.exp(-0.5*(i-512)**2/256**2) for i in range(1024)])
subject = dtw.TimeSeries([m.cos(i/128.0+1) * m.exp(-0.5*(i-512)**2/256**2) for i in range(1024)])

print "EUCLIDEAN\n"

print "Euclidean DTW:", dtw.dtw(query, subject, True)
print "Euclidean constrained DTW:", dtw.cdtw(query, subject, 32, True)
print "Euclidean L_2-norm:", dtw.euclidean(query, subject)

print "\nMANHATTEN\n"

print "Manhatten DTW:", dtw.dtw(query, subject)
print "Manhatten constrained DTW:", dtw.cdtw(query, subject, 32)
print "Manhatten L_1-norm:", dtw.manhatten(query, subject)

try:
    import pylab as pl
    pl.plot(query)
    pl.plot(subject)
    pl.show()
except:
    print("install matplotlib for graphical output!")
