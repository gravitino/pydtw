import numpy as np
import pylab as pl
import libdtw as dtw

# generate a cosine wave as query
T = np.linspace(0, 10, 1024)
Q = dtw.TimeSeries(np.cos(T))
L = dtw.TimeSeries()
U = dtw.TimeSeries()

# calculate envelope with window length w
w = 50
dtw.lb_envelope(Q, w, L, U)

# plot the envelope for the query
pl.plot(Q)
pl.plot(L)
pl.plot(U)
pl.show()

# generate a sine wave as subject
S = dtw.TimeSeries(np.sin(T))

# mode (Euclidean flavoured: True, Manhatten flavoured: False)
mode = True

# calculate lower bounds and associated constrained DTW measure
print "lb_Keogh on query", \
      dtw.lb_keogh_onQuery(Q, S, w, mode), \
      dtw.lb_keogh_onEnvelope(S, L, U, w, mode)

# now on subject
dtw.lb_envelope(S, w, L, U)
print "lb_Keogh on subject", \
      dtw.lb_keogh_onSubject(Q, S, w, mode), \
      dtw.lb_keogh_onEnvelope(Q, L, U, w, mode)
      
# the constrained DTW measure
print "constrained DTW", \
      dtw.dist_cdtw(Q, S, w, mode)
