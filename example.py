import libdtw as dtw

query = dtw.TimeSeries([1,2,1,2])
subject = dtw.TimeSeries([0,0,0,0])

print "Euclidian:", dtw.dtw(query, subject, True)
print "Manhatten:", dtw.dtw(query, subject, False)
