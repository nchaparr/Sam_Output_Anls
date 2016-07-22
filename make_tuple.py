from collections import namedtuple

def make_tuple(in_dict,tupname='values'):
    """
    make a named tuple from a dictionary

    Parameters
    ==========

    in_dict: dictionary
         Any python object with key/value pairs

    tupname: string
         optional name for the new namedtuple type

    Returns
    =======

    the_tup: namedtuple
          named tuple with keys as attributes
    """
    the_tup = namedtuple(tupname, in_dict.keys())
    the_tup = the_tup(**in_dict)
    return the_tup
