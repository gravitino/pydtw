import libdtw as dtw
import pylab as pl
import math as m

# define query and subject as shifted cosine waves
query = dtw.TimeSeries([m.cos(i/128.0) * m.exp(-0.5*(i-512)**2/256**2) for i in range(1024)])
subject = dtw.TimeSeries([m.cos(i/128.0+1) * m.exp(-0.5*(i-512)**2/256**2) for i in range(1024)])

pl.title("query and subject")
pl.plot(query)
pl.plot(subject)
pl.show()

for mode, name in [(True, 'Euclidean'), (False, 'Manhatten')]:

    print  dtw.dist_euclidean(query, subject) \
           if mode == True else \
           dtw.dist_manhatten(query, subject)
    
    for window in range(0, max(len(subject), len(query)), 32):
        gamma = dtw.WarpingPath()
        print dtw.dist_cdtw_backtrace(query, subject, window, gamma, mode)
        pl.plot(*zip(*[node[::-1] for node in gamma]))

    print dtw.dist_dtw(query, subject, mode)

    pl.imshow([[abs(query[i]-subject[j]) for j in range(len(subject))] 
                                         for i in range(len(query))], aspect="auto")
    pl.title(name)
    pl.show()

