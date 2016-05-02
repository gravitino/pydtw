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

import numpy as np
import pydtw as pd

# try different distance/similarity measures
for dtype, measure, name in [(np.float32, pd.host.lockstepEuclidean1f, "E1f"),
                             (np.float32, pd.host.lockstepManhattan1f, "M1f"),
                             (np.float64, pd.host.lockstepEuclidean1d, "E1D"),
                             (np.float64, pd.host.lockstepManhattan1d, "M1D")]:

    length = 1 << 27
    X = np.zeros(length, dtype=dtype)
    Y = np.ones (length, dtype=dtype)

    print name, measure(X, Y)
