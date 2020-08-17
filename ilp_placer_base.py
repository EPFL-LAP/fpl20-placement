"""Module holding basic ILP constraints that are often reused.
"""

##########################################################################
def linearize(a, b, y):
    """Generates linearization constraints for the product y = ab.

    Parameters
    ----------
    a : str
        First factor.
    b : str
        Second factor.
    y : Product.

    Returns
    -------
    List[str]
        A list holding the three linearization constraints.
    """
    
    force_zero_a = "- %s + %s >= 0" % (y, a)
    force_zero_b = "- %s + %s >= 0" % (y, b)
    force_one = "- %s + %s + %s <= 1" % (y, a, b)
   
    return [force_zero_a, force_zero_b, force_one] 
##########################################################################

##########################################################################
def binarize(a, a_bin, M):
    """Generates binarization constraints for variable a.

    Parameters
    ----------
    a : str
        Variable to be binarized.
    a_bin : str
        Identifier of the binarized variable.
    M : int
        An upper bound on a.

    Returns
    -------
    List[str]
        A list holding the two binarization constraints.
    """
    
    lb = "- %s + %s <= 0" % (a, a_bin)
    ub = "- %s + %d %s >= 0" % (a, M, a_bin)

    return [lb, ub]
##########################################################################

##########################################################################
def geq_indicate(var, indicator, var_max, thr):
    """Generates constraints that make indicator 1 iff var >= thr, else 0.

    Parameters
    ----------
    var : str
        Variable on which thresholding is performed.
    indicator : str
        Identifier of the indicator variable.
    var_max : int
        An upper bound on var.
    the : int
        Comparison threshold.

    Returns
    -------
    List[str]
        A list holding the two constraints.
    """
    
    lb = "- %s + %d %s <= 0" % (var, thr, indicator)
    ub = "- %s + %d %s >= -%d" % (var, var_max - thr + 1, indicator, thr - 1)

    return [lb, ub]
##########################################################################

##########################################################################
def negate(a, a_neg):
    """Generates a negation constraint for variable a.

    Parameters
    ----------
    a : str
        Variable to be negated.
    a_neg : str
        Identifier of the negated variable.

    Returns
    -------
    str
        Text of the constraint.
    """
    
    return "%s + %s = 1" % (a, a_neg)
##########################################################################
