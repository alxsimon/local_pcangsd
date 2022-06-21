# Wrapper functions to use the corners function from the R lostruct package

from rpy2.robjects.packages import STAP
import rpy2.robjects as robj
import numpy as np
import local_pcangsd as lp
import os

rsource = f"{os.path.dirname(lp.__file__)}/corners.R"

with open(rsource, "r") as f:
    string = f.read()
lo_corners = STAP(string, "lo_corners")


def corners(xy, prop, k=3):
    # Convert xy to R matrix
    nr, nc = xy.shape
    xvec = robj.FloatVector(xy.reshape(xy.size))
    xyr = robj.r.matrix(xvec, nrow=nr, ncol=nc, byrow=True)
    # Get the corners from R function
    corners_r = lo_corners.corners(xyr, prop=prop, k=k)
    # Convert back result to numpy.array
    corners_py = np.array(list(corners_r)).reshape(corners_r.dim, order="F") - 1
    # transforms index from R to Python by subtracting 1

    return corners_py
