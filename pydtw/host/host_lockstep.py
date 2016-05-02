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

import numpy as np
import pydtw as pd
import ctypes

def lockstepEuclidean(series0, series1):
    """This method computes the Euclidean distance for two time series of same
       shape and data type. Use low-level API for faster calls."""

    # make sure both time series have the same shape and datatype
    assert(series0.shape == series1.shape)
    assert(series0.dtype == series1.dtype)
    assert(0 < len(series0.shape) < 3)

    # remember the shape and flatten
    shape = series0.shape
    dtype = series0.dtype
    series0 = series0.reshape((np.prod(shape),))
    series1 = series1.reshape((np.prod(shape),))

    # arbitrary shape double-precision calls
    if dtype == np.float64:
        if len(shape) ==1 or shape[1] == 1:
            return pd.host.lockstepEuclidean1d(series0, series1)
        if shape[1] == 2:
            return pd.host.lockstepEuclidean2d(series0, series1)
        if shape[1] == 3:
            return pd.host.lockstepEuclidean3d(series0, series1)
        if shape[1] == 4:
            return pd.host.lockstepEuclidean4d(series0, series1)
        return pd.host.lockstepEuclideanNd(series0, series1, shape[1])

    # arbitrary shape single-precision calls
    if dtype == np.float32:
        if len(shape) ==1 or shape[1] == 1:
            return pd.host.lockstepEuclidean1f(series0, series1)
        if shape[1] == 2:
            return pd.host.lockstepEuclidean2f(series0, series1)
        if shape[1] == 3:
            return pd.host.lockstepEuclidean3f(series0, series1)
        if shape[1] == 4:
            return pd.host.lockstepEuclidean4f(series0, series1)
        return pd.host.lockstepEuclideanNf(series0, series1, shape[1])

    raise TypeError("ERROR: Type %s not supported." % series0.dtype)

def lockstepManhattan(series0, series1):
    """This method computes the Manhattan distance for two time series of same
       shape and data type. Use low-level API for faster calls."""

    # make sure both time series have the same shape and datatype
    assert(series0.shape == series1.shape)
    assert(series0.dtype == series1.dtype)
    assert(0 < len(series0.shape) < 3)

    # remember the shape and flatten
    shape = series0.shape
    series0 = series0.reshape((np.prod(shape),))
    series1 = series1.reshape((np.prod(shape),))

    # arbitrary shape double-precision calls
    if series0.dtype == np.float64:
        if len(shape) ==1 or shape[1] == 1:
            return pd.host.lockstepManhattan1d(series0, series1)
        if shape[1] == 2:
            return pd.host.lockstepManhattan2d(series0, series1)
        if shape[1] == 3:
            return pd.host.lockstepManhattan3d(series0, series1)
        if shape[1] == 4:
            return pd.host.lockstepManhattan4d(series0, series1)
        return pd.host.lockstepManhattanNd(series0, series1, shape[1])

    # arbitrary shape single-precision calls
    if series0.dtype == np.float32:
        if len(shape) ==1 or shape[1] == 1:
            return pd.host.lockstepManhattan1f(series0, series1)
        if shape[1] == 2:
            return pd.host.lockstepManhattan2f(series0, series1)
        if shape[1] == 3:
            return pd.host.lockstepManhattan3f(series0, series1)
        if shape[1] == 4:
            return pd.host.lockstepManhattan4f(series0, series1)
        return pd.host.lockstepManhattanNf(series0, series1, shape[1])

    raise TypeError("ERROR: Type %s not supported." % series0.dtype)
