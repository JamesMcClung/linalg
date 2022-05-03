#pragma once

#include <iostream>
#include <string>
#include <type_traits>
#include <utility>

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DESCRIPTION                             //
////////////////////////////////////////////////////////////////////////

// An implementation of LU decomposition.
// A matrix is decomposed in-place or not depending on whether an
//   r-value or l-value is passed to the constructor, respectively.

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

template <class SquareMatrix>
class LU_Decomp {
    static_assert(SquareMatrix::nrows == SquareMatrix::ncols, "Cannot LU-decompose a non-square matrix.");

   private:
    using Real = typename SquareMatrix::real_t;
    static constexpr int dim = SquareMatrix::nrows;

    SquareMatrix l_and_u;

    Real &l_ref(int r, int c);
    Real &u_ref(int r, int c);

    void init();

   public:
    LU_Decomp(const SquareMatrix &mat);
    LU_Decomp(SquareMatrix &&mat);

    template <class SomeMatrix>
    SomeMatrix solve(SomeMatrix rhs);

    Real l(int r, int c) const;
    Real u(int r, int c) const;
};

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

template <class SquareMatrix>
void LU_Decomp<SquareMatrix>::init() {
    for (int i = 0; i < dim; i++) {
        if (l_and_u(i, i) == Real(0)) {
            throw std::domain_error("Cannot LU-decompose a matrix with a 0 on the diagonal");
        }
        for (int r = i + 1; r < dim; r++) {
            l_ref(r, i) /= u(i, i);
            for (int c = i + 1; c < dim; c++) {
                l_and_u(r, c) -= l(r, i) * u(i, c);
            }
        }
    }
}

template <class SquareMatrix>
LU_Decomp<SquareMatrix>::LU_Decomp(const SquareMatrix &mat) : l_and_u(mat) {
    init();
}

template <class SquareMatrix>
LU_Decomp<SquareMatrix>::LU_Decomp(SquareMatrix &&mat) : l_and_u(std::move(mat)) {
    init();
}

template <class SquareMatrix>
template <class SomeMatrix>
SomeMatrix LU_Decomp<SquareMatrix>::solve(SomeMatrix rhs) {
    static_assert(SomeMatrix::nrows == SquareMatrix::nrows, "RHS must have same number of rows as decomposed matrix.");
    static_assert(std::is_same_v<typename SomeMatrix::real_t, typename SquareMatrix::real_t>, "RHS must have same element type as decomposed matrix.");

    // forward substitution
    for (int c = 0; c < rhs.ncols; c++) {
        for (int r = 0; r < rhs.nrows; r++) {
            for (int i = 0; i < r; i++) {
                rhs(r, c) -= l(r, i) * rhs(i, c);
            }
        }
    }

    // backward substitution
    for (int c = 0; c < rhs.ncols; c++) {
        for (int r = dim - 1; r >= 0; r--) {
            for (int i = dim - 1; i > r; i--) {
                rhs(r, c) -= u(r, i) * rhs(i, c);
            }
            rhs(r, c) /= u(r, r);
        }
    }

    return rhs;
}

template <class SquareMatrix>
typename LU_Decomp<SquareMatrix>::Real LU_Decomp<SquareMatrix>::l(int r, int c) const {
    if (r == c)
        return Real(1);
    if (r < c)
        return Real(0);
    return l_and_u(r, c);
}

template <class SquareMatrix>
typename LU_Decomp<SquareMatrix>::Real LU_Decomp<SquareMatrix>::u(int r, int c) const {
    if (r > c)
        return Real(0);
    return l_and_u(r, c);
}

template <class SquareMatrix>
typename LU_Decomp<SquareMatrix>::Real &LU_Decomp<SquareMatrix>::l_ref(int r, int c) {
    if (r <= c)
        throw std::out_of_range("Diagonal and upper triangle of L are fixed at index " + std::to_string(r) + "," + std::to_string(c));
    return l_and_u(r, c);
}

template <class SquareMatrix>
typename LU_Decomp<SquareMatrix>::Real &LU_Decomp<SquareMatrix>::u_ref(int r, int c) {
    if (r > c)
        throw std::out_of_range("Lower triangle of U is fixed at index " + std::to_string(r) + "," + std::to_string(c));
    return l_and_u(r, c);
}

}  // namespace linalg