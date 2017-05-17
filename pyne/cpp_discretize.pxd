"""C++ wrapper for discretize geometry header."""
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "discretize.h" namespace "pyne":
    ctypedef struct disc_result:
        int id
        int cell
        double vol_frac;
        double rel_error;

    cdef cppclass MeshRow:
        int num_rays
        bool grid
        vector[disc_result] sums
        double start_point_x, start_point_y
        double d1div1, d1div2, d2div1, d2div2
        vector[double] d3divs
        # ~MeshRow() except +
        MeshRow(int, bool) # except +
        void setDimension1(double, double) except +
        void setDimension2(double, double) except +
        void fireRays(int, vector[double]) except +
        vector[disc_result] getSums() except +
    vector[disc_result] discretize_geom(vector[vector[double]],
                                        int, bool) except +
