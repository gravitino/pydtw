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

T = np.linspace(0, 1, 1024)
X = np.cos(4*np.pi*T)
Y = np.cos(4*np.pi*T**0.5)

pl.plot(T, X)
pl.plot(T, Y)
pl.show()

DTWpure = pd.host.elasticEuclideanDTW1d
DTWback = pd.host.elasticEuclideanDTW1dBacktrace
Path = pd.host.WarpingPath()

print DTWpure(X, Y)
print DTWback(X, Y, Path)

I, J = zip(*list(Path))

pl.plot(X[np.array(I)])
pl.plot(Y[np.array(J)])
pl.show()
