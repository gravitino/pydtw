import libdtw as dtw
import pylab as pl
import math as m

# define query and subject as shifted cosine waves
query = dtw.TimeSeries([m.cos(i/64.0) * m.exp(-0.5*(i-512)**2/256**2) for i in range(1024)])
subject = dtw.TimeSeries([m.cos(i/64.0+1) * m.exp(-0.5*(i-512)**2/256**2) for i in range(1024)])

# plot query and subject for a quick overview
pl.title("query and subject")
pl.plot(query)
pl.plot(subject)
pl.show()

# do for Euclidean and Manhatten mode
for mode, name in [(True, 'Euclidean'), (False, 'Manhatten')]:

    # calculate naive L_p-norm
    print  dtw.dist_euclidean(query, subject) \
           if mode == True else \
           dtw.dist_manhatten(query, subject)
    
    # calculate a bunch of windowed DTWs (error should decrease monotonically)
    for window in range(0, max(len(subject), len(query)), 32):
        gamma = dtw.WarpingPath()
        print dtw.dist_cdtw_backtrace(query, subject, window, gamma, mode)
        pl.plot(*zip(*[node[::-1] for node in gamma]))

    # calculate full DTW
    print dtw.dist_dtw(query, subject, mode)

    # plot objective function H(i,j) := |query(i)-subject(j)|
    pl.imshow([[abs(query[i]-subject[j]) for j in range(len(subject))] 
                                         for i in range(len(query))], aspect="auto")
    pl.title(name)
    pl.show()

    # draw explicit alignment for full dtw
    pl.plot(query)
    pl.plot(subject)
    
    for i, j in gamma[0:len(gamma):4]:
        pl.plot([i, j], [query[i], subject[j]], c="grey")
    
    pl.show()
