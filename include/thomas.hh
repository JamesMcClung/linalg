#pragma once

#include <type_traits>

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DESCRIPTION                             //
////////////////////////////////////////////////////////////////////////

// An implementation of the Thomas algorithm for solving a tridiagonal
//   system.

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

template <class TridiagMatrix, class SomeMatrix>
SomeMatrix thomas(const TridiagMatrix &bmat, SomeMatrix rhs);

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

template <class TridiagMatrix, class SomeMatrix>
SomeMatrix thomas(const TridiagMatrix &bmat, SomeMatrix rhs) {
    static_assert(TridiagMatrix::nrows == TridiagMatrix::ncols, "Coefficient matrix must be square.");
    static_assert(TridiagMatrix::nrows == SomeMatrix::nrows, "Coefficient matrix and rhs must have same number of rows.");
    static_assert(std::is_same_v<typename TridiagMatrix::real_t, typename SomeMatrix::real_t>, "Coefficient matrix and rhs must be over the same field.");

    using Real = typename TridiagMatrix::real_t;
    Real altered_diag[bmat.nrows];
    for (int c = 0; c < rhs.ncols; c++) {
        // forward sweep
        altered_diag[0] = bmat.get(0, 0);
        for (int r = 1; r < bmat.nrows; r++) {
            Real w = bmat.get(r, r - 1) / altered_diag[r - 1];
            altered_diag[r] = bmat.get(r, r) - w * bmat.get(r - 1, r);
            rhs(r, c) -= w * rhs.get(r - 1, c);
        }
        // back substitution
        rhs(rhs.nrows - 1, c) /= altered_diag[bmat.nrows - 1];
        for (int r = rhs.nrows - 2; r >= 0; r--) {
            (rhs(r, c) -= bmat.get(r, r + 1) * rhs.get(r + 1, c)) /= altered_diag[r];
        }
    }
    return rhs;
}

}  // namespace linalg