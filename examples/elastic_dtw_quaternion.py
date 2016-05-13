# -*- coding: utf-8 -*-

# this is usually not needed, just make sure pydtw is in your path
import os.path, sys
parentdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)
sys.path.append(parentdir)

# pydtw integrates with numpy
import numpy as np
import pydtw as pd
import pylab as pl

#reads csv-File and returns an np.array (n-dims)
def read(filename, columns, splitter=None, verbose=False):
    print "reading file '", filename, "'..."
    result = []
    with open(filename, "r") as f:
        for line_str in f:
            try:
               line = line_str.split(splitter)
               row = []
               for c in columns:
                   row.append(line[c])
               result.append(map(float, row))
            except:
                if verbose:
                    print "ignored line: " + line;
    print "...done."
    return np.array(result)

def plot_range(subject, start, length, stride):
    end = start + length
    w = subject[start*stride+0:end*stride+0:stride]
    x = subject[start*stride+1:end*stride+1:stride]
    y = subject[start*stride+2:end*stride+2:stride]
    z = subject[start*stride+3:end*stride+3:stride]
    pl.plot(w)
    pl.plot(x)
    pl.plot(y)
    pl.plot(z)


## Main-------------------------------------------------------------------------

stride = 4

shape_indices_open_fridge =  [[1546, 1631], [4315, 4388], [6721, 6794], [9244, 9318], [12937, 13004], [15487, 15555], [18047, 18120], [20623, 20694], [23407, 23475], [26104, 26169], [28372, 28443], [30712, 30780], [33058, 33129], [35482, 35559], [38218, 38283], [40606, 40680], [43054, 43143], [45484, 45562], [48004, 48084], [50557, 50637]]
shape_indices_close_fridge = [[1636, 1718], [4398, 4464], [6795, 6867], [9320, 9402], [13009, 13089], [15565, 15635], [18124, 18198], [20698, 20780], [23485, 23560], [26190, 26249], [28444, 28515], [30787, 30861], [33133, 33210], [35560, 35640], [38284, 38360], [40687, 40761], [43150, 43233], [45565, 45645], [48085, 48159], [50638, 50712]]
## =========
## read data
## =========

subject_raw = read("Y:/MTMRepository/__BIG_DATA__/Opportunity/S1-Drill.dat", [72,73,74,75]) # RLA quaternion

# map data to [-1,1]

norms = np.sqrt(np.sum(subject_raw**2, axis=1))
subject_raw[:, 0] /= norms
subject_raw[:, 1] /= norms
subject_raw[:, 2] /= norms
subject_raw[:, 3] /= norms

S = subject_raw.flatten()

QUERY_NUMBER = 8

start = shape_indices_open_fridge[QUERY_NUMBER][0]
end =   shape_indices_open_fridge[QUERY_NUMBER][1]
Q = np.copy(S[start*stride:end*stride]) # COPY (!!!) query

L = len(Q)/stride
N = len(S)/stride
w = L/10

# alienate query
p, sigma = 1.5, 0.05
T = np.linspace(0, 1, L)
step_size = 1 # -1 => invert shape
Q[0::stride] = np.interp(T**p, T, Q[0::stride])[::step_size]
Q[1::stride] = np.interp(T**p, T, Q[1::stride])[::step_size]
Q[2::stride] = np.interp(T**p, T, Q[2::stride])[::step_size]
Q[3::stride] = np.interp(T**p, T, Q[3::stride])[::step_size]
# restore normalization
Q += np.random.normal(0, sigma, len(Q))
for i in range(L):
    norm = np.sqrt(np.sum(Q[i*stride:(i+1)*stride]**2))
    Q[i*stride + 0] /= norm
    Q[i*stride + 1] /= norm
    Q[i*stride + 2] /= norm
    Q[i*stride + 3] /= norm


DTW = pd.host.elasticQuaternionCDTWd
LB_Kim = pd.host.elasticLBKim4d
offset = 0

counter_Kim = 0
counter_DTW = 0
counter_total = 0
best_value = float("infinity")
best_index = -1
for i in range(offset, N-L):     
    if (i % 10000 == 0):
        print i, "/", N-L, "(", 100.0*i/(N-L), "%)"
        print best_value, "@ [", best_index, ",", best_index + L, "]"
    start = i
    end = i + L
    
    C = S[start*stride:end*stride]

    value = LB_Kim(Q,C)
    if (value < best_value):
        value = DTW(Q,C,w)
        counter_DTW += 1
        if (value < best_value):
            best_value = value
            best_index = i
    else:
        counter_Kim += 1
    
    counter_total += 1

print best_value, "@ [", best_index, ",", best_index + L, "]"
print "DTW: ", counter_DTW
print "Kim: ", counter_Kim
print "total: ", counter_total

## ==========
## upper plot
## ==========
ax1 = pl.subplot(211)
pl.title("Subject")
plot_range(S, best_index, L, stride)


## ==========
## lower plot
## ==========
ax2 = pl.subplot(212, sharex=ax1, sharey=ax1)
pl.title("Query")
plot_range(Q, 0, L, stride)

pl.show()