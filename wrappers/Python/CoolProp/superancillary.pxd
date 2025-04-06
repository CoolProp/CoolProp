from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.string cimport string


cdef extern from "superancillary/superancillary.h" namespace "CoolProp::superancillary":
    cdef struct MonotonicExpansionMatch:
        const size_t idx # The index of the expansion that has been matched
        const double ymin, ymax, xmin, xmax 
        # Check if a value of the dependent variable is within this match
        # bool contains_y (double) const 

    # Data associated with a monotonic interval
    cdef struct IntervalMatch:
        const vector[MonotonicExpansionMatch] expansioninfo # The information about the expansions for this interval
        const double xmin, xmax, ymin, ymax
        # Check if a value of the dependent variable is within this interval
        # bool contains_y (double) const

    cdef cppclass ChebyshevExpansion[ArrayType]:
        ChebyshevExpansion() except +ValueError
        ChebyshevExpansion(double, double, ArrayType) except +ValueError
        const double xmin() const
        const double xmax() const
        const ArrayType& coeff() const
        U eval[U](const U& x) const
        void eval_manyC[U](const U x[], U v[], size_t) const
        double solve_for_x(double, double, double, unsigned int, size_t, double) const
        void solve_for_x_manyC[U](const U[], size_t, double, double, unsigned int, size_t, double, U[], U[]) const

    cdef cppclass ChebyshevApproximation1D[ArrayType]:
        ChebyshevApproximation1D(vector[ChebyshevExpansion[ArrayType]]&&) except +ValueError
        bool is_monotonic() const
        double xmin()
        double xmax()
        void eval_manyC[U](const U x[], U v[], size_t) const
        const vector[ChebyshevExpansion[ArrayType]] m_expansions # The collection of expansions forming the approximation
        const vector[double] m_x_at_extrema # The values of the independent variable at the extrema of the expansions
        const vector[IntervalMatch] get_monotonic_intervals() # The intervals that are monotonic
        vector[pair[double, int]] get_x_for_y(double, unsigned int, size_t, double) const
        void count_x_for_y_manyC[U](const U[], size_t, unsigned int, size_t, double, U[]) const

    cdef cppclass SuperAncillary[ArrayType]:
        SuperAncillary(const string&) except +ValueError
        double eval_sat(double, char, short) except +ValueError
        void eval_sat_manyC[U](const U[], size_t, char, short, U[]) except +ValueError
    