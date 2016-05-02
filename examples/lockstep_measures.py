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

import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

import numpy as np
import pydtw as pd

for dtype, measure, name in [(np.float32, pd.host.lockstepEuclidean1f, "E1f"),
                             (np.float32, pd.host.lockstepManhattan1f, "M1f"),
                             (np.float64, pd.host.lockstepEuclidean1d, "E1D"),
                             (np.float64, pd.host.lockstepManhattan1d, "M1D")]:

    X = np.zeros(1024, dtype=dtype)
    Y = np.linspace(0, 1, 1024, dtype=dtype)

    print name, measure(X, Y)
