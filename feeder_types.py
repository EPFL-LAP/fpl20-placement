"""Module holding type declarations used by feeder.

Types
-----
Ble
    A pair representing a BLE.

Methods
-------
A legality checker for each of the types.

ble_empty(ble : Ble)
    Checks if the BLE is empty.
"""

from collections import namedtuple

Ble = namedtuple("Ble", ["lut", "ff"])
##########################################################################
def check_ble(a):
    """Checks if a conforms to the specification of the Ble type.
    
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
        If a is not Ble.
        If a.lut is not a string.
        If a.ff is not a string.
    """

    assert (isinstance(a, Ble)), "a not Ble."

    assert (isinstance(a.lut, str)), "LUT not a string."
    assert (isinstance(a.ff, str)), "FF not a string."
##########################################################################


##########################################################################
def ble_empty(ble):
    """Checks if the BLE is empty.

    Parameters
    ----------
    ble : Ble

    Returns
    -------
    bool
        True if ble is empty, False otherwise.

    Raises
    ------
    AssertionError
        If ble is not a legal Ble.
    """

    check_ble(ble)

    return (not ble.lut) and (not ble.ff)
##########################################################################
