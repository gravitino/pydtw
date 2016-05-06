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
length = 1024
window = int(0.2*length) # 20% relative warping window

# create two spirals (array of 2D struct)
T = np.linspace(0, 1, length)
X = np.zeros(2*length)
Y = np.zeros(2*length)

X[0::2] = T*np.cos(6*np.pi*T)
X[1::2] = T*np.sin(6*np.pi*T)

Y[0::2] = T*np.cos(6*np.pi*T**0.5)
Y[1::2] = T*np.sin(6*np.pi*T**0.5)

# compute best DTW alignment in 2D
DTWback = pd.host.elasticEuclideanCDTW2dBacktrace
DTWpath = pd.host.WarpingPath()
distance = DTWback(X, Y, window, DTWpath)
print distance

# compute warping envelopes
LX = np.zeros(X.shape)
UX = np.zeros(X.shape)
LY = np.zeros(Y.shape)
UY = np.zeros(Y.shape)

# compute warping envelopes
pd.host.elasticWarpingEnvelopeNd(X, LX, UX, window, 2)
pd.host.elasticWarpingEnvelopeNd(Y, LY, UY, window, 2)

# plot the alignment
pl.figure(1)
for i, j in DTWpath:
    pl.plot([X[2*i], Y[2*j]], [X[2*i+1], Y[2*j+1]], color="grey")
pl.plot(X[0::2], X[1::2])
pl.plot(Y[0::2], Y[1::2])

# plot warping envelopes
pl.figure(2)
pl.subplot(221)
pl.plot( X[0::2])
pl.plot(LX[0::2])
pl.plot(UX[0::2])

pl.subplot(222)
pl.plot( X[1::2])
pl.plot(LX[1::2])
pl.plot(UX[1::2])

pl.subplot(223)
pl.plot( Y[0::2])
pl.plot(LY[0::2])
pl.plot(UY[0::2])

pl.subplot(224)
pl.plot( Y[1::2])
pl.plot(LY[1::2])
pl.plot(UY[1::2])

pl.show()
