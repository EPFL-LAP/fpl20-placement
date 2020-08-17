"""Module holding type declarations used by ilp_placer.

Types
-----
BasicEdge
    A pair representing an edge.
NodePos
    Fully specifies a node position on the FPGA grid,
    by listing the x and y coordinates of its cluster
    and the LUT-index within it.
NodeCls
    Only specifies the cluster of a node on the FPGA grid.
PhysDirCon
    Specifies a concrete direct connection of the FPGA grid.
FreeDirCon
    Specifies a direct connection in the form of an edge of
    the pattern-representing static graph.
CovPair
    Specifies a concrete direct connection of the FPGA grid
    together with its delay.
PhysProgCon
    Specifies a concrete programmable connection of the FPGA grid.
FreeProgCon
    Specifies a programmable connection in the form of a free vector
    as in the VPR's delay lookup table.
StretchPair
    Specifies a concrete programmable connection of the FPGA grid
    together with its delay.

Methods
-------
A legality checker for each of the types.
"""
 
from collections import namedtuple
from config_vars import N

#Type declarations. Each one is follow by a legality checker.

BasicEdge = namedtuple("BasicEdge", ['u', 'v'])
##########################################################################
def check_basic_edge(a):
    """Checks if a conforms to the specification of the BasicEdge type.
    
    Parameters
    ----------
    a : Any
        The variable to be checked.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If a is not BasicEdge.
        If a.u is not a nonnegative integer.
        If a.v is not a nonnegative integer.
    """

    assert (isinstance(a, BasicEdge)), "a not BasicEdge."
    check_nonneg_int(a.u)
    check_nonneg_int(a.v)
##########################################################################

NodePos = namedtuple("NodePos", ['x', 'y', 'i'])
##########################################################################
def check_node_pos(a):
    """Checks if a conforms to the specification of the NodePos type.
    
    Parameters
    ----------
    a : Any
        The variable to be checked.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If a is not NodePos.
        If a.x is not a nonnegative integer.
        If a.y is not a nonnegative integer.
        If a.i is not a nonnegative integer.
        If a.i > N.
    """

    assert (isinstance(a, NodePos)), "a not NodePos."
    check_nonneg_int(a.x)
    check_nonneg_int(a.y)
    check_nonneg_int(a.i)
    assert (a.i < N), "LUT-index exceeds cluster capacity."
##########################################################################

NodeCls = namedtuple("NodeCls", ['x', 'y'])
##########################################################################
def check_node_cls(a):
    """Checks if a conforms to the specification of the NodeCls type.
    
    Parameters
    ----------
    a : Any
        The variable to be checked.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If a is not NodePos.
        If a.x is not a nonnegative integer.
        If a.y is not a nonnegative integer.
    """

    assert (isinstance(a, NodeCls)), "a not NodeCls."
    check_nonneg_int(a.x)
    check_nonneg_int(a.y)
##########################################################################

PhysDirCon = namedtuple("PhysDirCon", ['u', 'v'])
##########################################################################
def check_phys_dir_con(a):
    """Checks if a conforms to the specification of the PhysDirCon type.
    
    Parameters
    ----------
    a : Any
        The variable to be checked.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If a is not PhysDirCon.
        If a.u is not a legal NodePos.
        If a.v is not a legal NodePos.
    """

    assert (isinstance(a, PhysDirCon)), "a not PhysDirCon."
    check_node_pos(a.u)
    check_node_pos(a.v)
##########################################################################

FreeDirCon = namedtuple("FreeDirCon", ['u', 'v', "x_offset", "y_offset"])
##########################################################################
def check_free_dir_con(a):
    """Checks if a conforms to the specification of the FreeDirCon type.
    
    Parameters
    ----------
    a : Any
        The variable to be checked.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If a is not FreeDirCon.
        If a.u is not a nonnegative integer.
        If a.u > N.
        If a.v is not a nonnegative integer.
        If a.v > N.
        If a.x_offset is not an integer.
        If a.y_offset is not an integer.
    """

    assert (isinstance(a, FreeDirCon)), "a not FreeDirCon."
    check_nonneg_int(a.u)
    assert (a.u < N), "Tail LUT-index exceeds cluster capacity."
    check_nonneg_int(a.v)
    assert (a.v < N), "Head LUT-index exceeds cluster capacity."

    assert (isinstance(a.x_offset, int)), "x-offset noninteger."
    assert (isinstance(a.y_offset, int)), "y-offset noninteger."
##########################################################################

CovPair = namedtuple("CovPair", ['u', 'v', "td"])
##########################################################################
def check_cov_pair(a):
    """Checks if a conforms to the specification of the CovPair type.
    
    Parameters
    ----------
    a : Any
        The variable to be checked.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If a is not CovPair.
        If a.u is not a legal NodePos.
        If a.v is not a legal NodePos.
        If a.td is not a nonnegative float.
    """

    assert (isinstance(a, CovPair)), "a not CovPair."
    check_node_pos(a.u)
    check_node_pos(a.v)
    check_nonneg_float(a.td)
##########################################################################

PhysProgCon = namedtuple("PhysProgCon", ['u', 'v'])
##########################################################################
def check_phys_prog_con(a):
    """Checks if a conforms to the specification of the PhysProgCon type.
    
    Parameters
    ----------
    a : Any
        The variable to be checked.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If a is not PhysProgCon.
        If a.u is not a legal NodeCls.
        If a.v is not a legal NodeCls.
    """

    assert (isinstance(a, PhysProgCon)), "a not PhysProgCon."
    check_node_cls(a.u)
    check_node_cls(a.v)
##########################################################################

FreeProgCon = namedtuple("FreeProgCon", ["x_offset", "y_offset"])
##########################################################################
def check_free_prog_con(a):
    """Checks if a conforms to the specification of the FreeProgCon type.
    
    Parameters
    ----------
    a : Any
        The variable to be checked.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If a is not FreeProgCon.
        If a.x_offset is not an integer.
        If a.y_offset is not an integer.
    """

    assert (isinstance(a, FreeProgCon)), "a not FreeProgCon."
    assert (isinstance(a.x_offset, int)), "x-offset noninteger."
    assert (isinstance(a.y_offset, int)), "y-offset noninteger."
##########################################################################

StretchPair = namedtuple("StretchPair", ['u', 'v', "td"])
##########################################################################
def check_stretch_pair(a):
    """Checks if a conforms to the specification of the StretchPair type.
    
    Parameters
    ----------
    a : Any
        The variable to be checked.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If a is not StretchPair.
        If a.u is not a legal NodeCls.
        If a.v is not a legal NodeCls.
        If a.td is not a nonnegative float.
    """

    assert (isinstance(a, StretchPair)), "a not StretchPair."
    check_node_cls(a.u)
    check_node_cls(a.v)
    check_nonneg_float(a.td)
##########################################################################

##########################################################################
def check_nonneg_int(a):
    """Checks if a is a nonnegative integer.

    Parameters
    ----------
    a : Any
        The variable to be checked.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If a is not a nonnegative float.
    """

    assert (isinstance(a, int)), "a not an integer."
    assert (a >= 0), "a negative."
##########################################################################

##########################################################################
def check_nonneg_float(a):
    """Checks if a is a nonnegative float.

    Integers are also accepted.

    Parameters
    ----------
    a : Any
        The variable to be checked.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If a is not a nonnegative float.
    """

    assert (isinstance(a, int) or isinstance(a, float)), "a not a float."
    assert (a >= 0), "a negative."
##########################################################################

#Conversions:

##########################################################################
def node_pos_to_node_cls(a):
    """Converts a from NodePos to NodeCls.

    Parameters
    ----------
    a : NodePos
        The variable to be converted.

    Returns
    -------
    NodeCls
        Converted value.
    """

    return NodeCls(a.x, a.y)
##########################################################################

##########################################################################
def phys_to_free_dir_con(a):
    """Converts a from PhysDirCon to FreeDirCon.

    Parameters
    ----------
    a : PhysDirCon
        The variable to be converted.

    Returns
    -------
    FreeDirCon
        Converted value.
    """

    return FreeDirCon(a.u.i, a.v.i, a.v.x - a.u.x, a.v.y - a.u.y)
##########################################################################

##########################################################################
def phys_to_free_prog_con(a):
    """Converts a from PhysProgCon to FreeProgCon.

    Parameters
    ----------
    a : PhysProgCon
        The variable to be converted.

    Returns
    -------
    FreeProgCon
        Converted value.
    """

    return FreeProgCon(a.v.x - a.u.x, a.v.y - a.u.y)
##########################################################################
