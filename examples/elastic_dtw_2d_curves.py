###############################################################################
#    This file is part of pydtw.
#
#    pydtw is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    pydtw is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with pydtw.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

# this is usually not needed, just make sure pydtw is in your path
import os.path, sys
parentdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)
sys.path.append(parentdir)

# pydtw integrates with numpy
import numpy as np
import pydtw as pd
import pylab as pl

# number of samples and warping window
length = 1 << 10
window = int(0.2*length) # 20% relative warping window

# create two spirals (array of 2D struct)
T = np.linspace(0, 1, length)
A = np.zeros(2*length)
B = np.zeros(2*length)

A[0::2] = T*np.cos(6*np.pi*T)
A[1::2] = T*np.sin(6*np.pi*T)

B[0::2] = T*np.cos(6*np.pi*T**0.5)
B[1::2] = T*np.sin(6*np.pi*T**0.5)

# compute best DTW alignment in 2D
DTWback = pd.host.elasticEuclideanCDTW2dBacktrace
DTWpath = pd.host.WarpingPath()
distance = DTWback(A, B, window, DTWpath)
print distance

# compute warping envelopes
LA = np.zeros(A.shape)
UA = np.zeros(A.shape)
LB = np.zeros(B.shape)
UB = np.zeros(B.shape)

# compute warping envelopes
pd.host.elasticWarpingEnvelopeNd(A, LA, UA, window, 2)
pd.host.elasticWarpingEnvelopeNd(B, LB, UB, window, 2)

# compute the matrix of residues
D = np.zeros((len(A)/2, len(B)/2))
pd.host.residuesMatrixEuclideanNd(A, B, D, 2);

# plot the alignment over spatial domain
pl.figure(1)
for i, j in DTWpath:
    pl.plot([A[2*i], B[2*j]], [A[2*i+1], B[2*j+1]], color="grey")
pl.plot(A[0::2], A[1::2])
pl.plot(B[0::2], B[1::2])

# plot alignment over time domain
pl.figure(2)

# the matrix of residues
pl.imshow(D)

# please not J: x-axis, I: y-axis
I, J = zip(*DTWpath)
pl.plot(J, I, color="red")

# plot warping envelopes
pl.figure(3)
pl.subplot(221)
pl.plot( A[0::2])
pl.plot(LA[0::2])
pl.plot(UA[0::2])

pl.subplot(222)
pl.plot( A[1::2])
pl.plot(LA[1::2])
pl.plot(UA[1::2])

pl.subplot(223)
pl.plot( B[0::2])
pl.plot(LB[0::2])
pl.plot(UB[0::2])

pl.subplot(224)
pl.plot( B[1::2])
pl.plot(LB[1::2])
pl.plot(UB[1::2])

pl.show()
