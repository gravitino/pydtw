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

# length of the compared arrays
length = (2*3*4*5)*(1 << 20)

# try different fixed-dimension distance/similarity measures
for dtype, measure, name in [(np.float64, pd.host.lockstepEuclidean1d, "E1d"),
                             (np.float64, pd.host.lockstepEuclidean2d, "E2d"),
                             (np.float64, pd.host.lockstepEuclidean3d, "E3d"),
                             (np.float64, pd.host.lockstepEuclidean4d, "E4d"),
                             (np.float32, pd.host.lockstepEuclidean1f, "E1f"),
                             (np.float32, pd.host.lockstepEuclidean2f, "E2f"),
                             (np.float32, pd.host.lockstepEuclidean3f, "E3f"),
                             (np.float32, pd.host.lockstepEuclidean4f, "E4f"),
                             (np.float64, pd.host.lockstepManhattan2d, "M1d"),
                             (np.float64, pd.host.lockstepManhattan1d, "M2d"),
                             (np.float64, pd.host.lockstepManhattan3d, "M3d"),
                             (np.float64, pd.host.lockstepManhattan4d, "M4d"),
                             (np.float32, pd.host.lockstepManhattan1f, "M1f"),
                             (np.float32, pd.host.lockstepManhattan2f, "M2f"),
                             (np.float32, pd.host.lockstepManhattan3f, "M3f"),
                             (np.float32, pd.host.lockstepManhattan4f, "M4f")]:

    X = np.zeros(length, dtype=dtype)
    Y = np.ones (length, dtype=dtype)*2

    print name, measure(X, Y)

# try different variable dimension distance/similarity measures
for dtype, measure, name in [(np.float64, pd.host.lockstepEuclideanNd, "ENd"),
                             (np.float32, pd.host.lockstepEuclideanNf, "ENf"),
                             (np.float64, pd.host.lockstepManhattanNd, "MNd"),
                             (np.float32, pd.host.lockstepManhattanNf, "MNf")]:

    X = np.zeros(length, dtype=dtype)
    Y = np.ones (length, dtype=dtype)*2

    # test variable dimension
    for dimension in [1, 2, 3, 4]:

        print name, dimension, measure(X, Y, dimension)
