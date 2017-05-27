"""C++ wrapper for discretize geometry header."""
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp cimport bool

cdef extern from "moab/Types.hpp" namespace "moab":

    ctypedef size_t EntityHandle
    ctypedef enum ErrorCode:
        pass

cdef extern from "discretize.h" namespace "pyne":

    vector[vector[double]] discretize_geom(vector[vector[double]],
                                        map[EntityHandle, int],
                                        int,
                                        bool) except +
