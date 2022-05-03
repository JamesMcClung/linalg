#pragma once

#include <ostream>

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DESCRIPTION                             //
////////////////////////////////////////////////////////////////////////

// A Matrix is an abstract base class for several implementations of
//   matrices. Matrix only enforces a getter for its elements.

// Classes that inherit from Matrix include:
//   - FullMatrix: stores every element on the stack
//   - ZeroMatrix: immutable, and all elements are 0
//   - RecursiveMatrix: allocates blocks when needed
//   - BandedMatrix: only stores elements within bandwidth
//   - PermutationMatrix: for efficient storage and multiplication
//   - FuzzyMatrix?: coarsens elements
//   - LazyMatrix?: lazily evaluated
//   - MatrixView?: a view of another matrix

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

template <int nr, int nc, typename Real>
class Matrix {
   public:
    using real_t = Real;
    static constexpr int nrows = nr;
    static constexpr int ncols = nc;

    virtual Real operator()(int i, int j) const = 0;
    Real get(int i, int j) const;  // an alias for operator()(int, int) const
};

template <int nr, int nc, typename Real>
bool operator==(const Matrix<nr, nc, Real> &m1, const Matrix<nr, nc, Real> &m2);
template <int nr, int nc, typename Real>
bool operator!=(const Matrix<nr, nc, Real> &m1, const Matrix<nr, nc, Real> &m2);

template <int nr, int nc, typename Real>
std::ostream &operator<<(std::ostream &os, const Matrix<nr, nc, Real> &m);

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

template <int nr, int nc, typename Real>
Real Matrix<nr, nc, Real>::get(int i, int j) const {
    return this->operator()(i, j);
}

template <int nr, int nc, typename Real>
bool operator==(const Matrix<nr, nc, Real> &m1, const Matrix<nr, nc, Real> &m2) {
    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nc; j++) {
            if (m1.get(i, j) != m2.get(i, j)) return false;
        }
    }
    return true;
}

template <int nr, int nc, typename Real>
bool operator!=(const Matrix<nr, nc, Real> &m1, const Matrix<nr, nc, Real> &m2) {
    return !(m1 == m2);
}

template <int nr, int nc, typename Real>
std::ostream &operator<<(std::ostream &os, const Matrix<nr, nc, Real> &m) {
    os << nr << "x" << nc << " " << typeid(Real).name() << " Matrix:" << std::endl;
    for (int r = 0; r < nr; r++) {
        os << "[";
        for (int c = 0; c < nc; c++) {
            os << " " << m.get(r, c);
        }
        os << " ]" << std::endl;
    }
    return os;
}

}  // namespace linalg