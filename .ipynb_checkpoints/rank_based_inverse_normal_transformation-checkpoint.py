# Rank-based inverse normal transformation
# This function taks in a np.arrary or list, returns INT values in a np.array
# Remove any missing values from input first!!!
# Converted from this in R: y <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
# Parameter: input_x is a np.array or list without missing values

import scipy.stats as ss
import numpy as np
def inverse_normal_transformation(input_x):
    """
    This function taks in a np.arrary or list, returns INT values in a np.array
    Remove any missing values from input first!!!
    Converted from R code: y <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
    Parameter: input_x is a np.array or list without missing values

    Example:
    >>> x = np.array([0.1, 0.3, 0.2, 0.5])
    >>> inv.inverse_normal_transformation(x)
    array([-1.15034938,  0.31863936, -0.31863936,  1.15034938])
    """
    input_x = np.array(input_x)
    return ss.norm.ppf((ss.rankdata(input_x)-0.5)/len(input_x))
