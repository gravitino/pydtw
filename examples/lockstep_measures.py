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

# length of the processed arrays
length = (2*3*4*5)*(1 << 20)
print "length=", length

# short cuts for the measures
Euc, Man = pd.lockstepEuclidean, pd.lockstepManhattan

# try different fixed-dimension distance/similarity measures
for dtype, measure, name, shape in [(np.float64, Euc, "E1d", (length/1, 1)),
                                    (np.float64, Euc, "E2d", (length/2, 2)),
                                    (np.float64, Euc, "E3d", (length/3, 3)),
                                    (np.float64, Euc, "E4d", (length/4, 4)),
                                    (np.float64, Euc, "E5d", (length/5, 5)),
                                    (np.float32, Euc, "E1f", (length/1, 1)),
                                    (np.float32, Euc, "E2f", (length/2, 2)),
                                    (np.float32, Euc, "E3f", (length/3, 3)),
                                    (np.float32, Euc, "E4f", (length/4, 4)),
                                    (np.float32, Euc, "E5f", (length/5, 5)),
                                    (np.float64, Man, "M1d", (length/1, 1)),
                                    (np.float64, Man, "M2d", (length/2, 2)),
                                    (np.float64, Man, "M3d", (length/3, 3)),
                                    (np.float64, Man, "M4d", (length/4, 4)),
                                    (np.float64, Man, "M5d", (length/5, 5)),
                                    (np.float32, Man, "M1f", (length/1, 1)),
                                    (np.float32, Man, "M2f", (length/2, 2)),
                                    (np.float32, Man, "M3f", (length/3, 3)),
                                    (np.float32, Man, "M4f", (length/4, 4)),
                                    (np.float32, Man, "M5f", (length/5, 5))]:

    X = np.zeros(length, dtype=dtype)
    Y = np.ones (length, dtype=dtype)*2

    print name, measure(X.reshape(shape), Y.reshape(shape))
