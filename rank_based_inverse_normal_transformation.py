# Rank-based inverse normal transformation
# This function taks in a np.arrary or list, returns INT values in a np.array
# Remove any missing values from input first!!!
# Converted from this in R: y <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
# Parameter: input_x is a np.array or list without missing values

import scipy.stats as ss
import numpy as np
def inverse_normal_transformation(input_x):
    """
    This function takes in a np.arrary or list, returns INT values in a np.array
    Remove any missing values from input first!!!
    Converted from R code: y <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
    Parameter: input_x is a np.array or list without missing values

    Note:
    May still have slight skewness after transformation.
    Check inverse_normal_transformation_blom() for better solutions

    Example:
    >>> x = np.array([0.1, 0.3, 0.2, 0.5])
    >>> inv.inverse_normal_transformation(x)
    array([-1.15034938,  0.31863936, -0.31863936,  1.15034938])
    """
    input_x = np.array(input_x)
    return ss.norm.ppf((ss.rankdata(input_x)-0.5)/len(input_x))

# Rank-based inverse normal transformed using the Blom transformation method
def inverse_normal_transformation_blom(input_x):
    '''
    This function takes in a np.arrary or list, returns INT values in a np.array
    Remove any missing values from input first!!!

    Notes:
    1 .This version use specific setting of some parameters,
        but overall is equivalent to inverse_normal_transformation().
    2. This version may create better normal dristributed values.

    
    Rank-based inverse normal transformed using the Blom transformation method with the recommended parameters:
        Y_i = Φ^{-1}((r_i - 3/8) / (N + 1/4))
    where r_i are ranks among the N non-missing observations.
    
    
    Reference:
    Beasley TM, Erickson S, Allison DB. Rank-based inverse normal transformations are increasingly used, but are they merited? Behav Genet. 2009 Sep;39(5):580-95. doi: 10.1007/s10519-009-9281-0. Epub 2009 Jun 14. PMID: 19526352; PMCID: PMC2921808.

    Param:
    - input_x: array of values to be transformed
    
    '''
    input_x = np.array(input_x)
    return ss.norm.ppf((ss.rankdata(input_x)-3/8)/(len(input_x)+1/4))