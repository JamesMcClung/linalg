#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>

#include "matrix.hh"
#include "permutation.hh"

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DESCRIPTION                             //
////////////////////////////////////////////////////////////////////////

// An implementation of LUP decomposition.
// A matrix is decomposed in-place or not depending on whether an
//   r-value or l-value is passed to the constructor, respectively.

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

template <class SquareMatrix>
class LUP_Decomp {
    static_assert(SquareMatrix::nrows == SquareMatrix::ncols, "Cannot LUP-decompose a non-square matrix.");

   private:
    using Real = typename SquareMatrix::real_t;
    static constexpr int dim = SquareMatrix::nrows;

    SquareMatrix l_and_u;
    PermutationMatrix<dim, Real> p;

    Real &l_ref(int r, int c);
    Real &u_ref(int r, int c);

    void init();

   public:
    LUP_Decomp(const SquareMatrix &mat);
    LUP_Decomp(SquareMatrix &&mat);

    template <class SomeMatrix>
    SomeMatrix solve(const SomeMatrix &rhs);

    Real l(int r, int c) const;
    Real u(int r, int c) const;

    SquareMatrix get_L() const {
        SquareMatrix L;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                L(i, j) = l(i, j);
        return L;
    }
    FullMatrix<dim, dim, Real> get_U() const {
        SquareMatrix U;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                U(i, j) = u(i, j);
        return U;
    }
    PermutationMatrix<dim, Real> get_P() const {
        return PermutationMatrix<dim, Real>(p);
    }
};

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

template <class SquareMatrix>
void LUP_Decomp<SquareMatrix>::init() {
    // adapted from https://lampx.tugraz.at/~hadley/num/ch2/2.3a.php
    for (int i = 0; i < dim; i++) {
        // pivot rows
        Real largest_u = 0;
        int row_with_largest_u;
        for (int j = i; j < dim; j++) {
            Real Uii = l_and_u(j, i);
            for (int q = 0; q < i; q++) {
                Uii -= l_and_u(j, q) * l_and_u(q, j);
            }
            if (std::abs(Uii) > largest_u) {
                largest_u = std::abs(Uii);
                row_with_largest_u = j;
            }
        }
        if (i != row_with_largest_u) {
            p.swapRows(i, row_with_largest_u);
            for (int c = 0; c < dim; c++) {
                std::swap(l_and_u(i, c), l_and_u(row_with_largest_u, c));
            }
        }

        // determine U across row i
        for (int c = i; c < dim; c++) {
            for (int j = 0; j < i; j++) {
                u_ref(i, c) -= l_and_u(i, j) * l_and_u(j, c);
            }
        }
        // determine L down column i
        for (int r = i + 1; r < dim; r++) {
            for (int j = 0; j < i; j++) {
                l_ref(r, i) -= l_and_u(r, j) * l_and_u(j, i);
            }
            l_ref(r, i) /= l_and_u(i, i);
        }
    }
}

template <class SquareMatrix>
LUP_Decomp<SquareMatrix>::LUP_Decomp(const SquareMatrix &mat) : l_and_u(mat) {
    init();
}

template <class SquareMatrix>
LUP_Decomp<SquareMatrix>::LUP_Decomp(SquareMatrix &&mat) : l_and_u(std::move(mat)) {
    init();
}

template <class SquareMatrix>
template <class SomeMatrix>
SomeMatrix LUP_Decomp<SquareMatrix>::solve(const SomeMatrix &rhs) {
    static_assert(SomeMatrix::nrows == SquareMatrix::nrows, "RHS must have same number of rows as decomposed matrix.");
    static_assert(std::is_same_v<typename SomeMatrix::real_t, typename SquareMatrix::real_t>, "RHS must have same element type as decomposed matrix.");

    auto prhs = p * rhs;

    // forward substitution
    for (int c = 0; c < prhs.ncols; c++) {
        for (int r = 0; r < prhs.nrows; r++) {
            for (int i = 0; i < r; i++) {
                prhs(r, c) -= l(r, i) * prhs(i, c);
            }
        }
    }

    // backward substitution
    for (int c = 0; c < prhs.ncols; c++) {
        for (int r = dim - 1; r >= 0; r--) {
            for (int i = dim - 1; i > r; i--) {
                prhs(r, c) -= u(r, i) * prhs(i, c);
            }
            prhs(r, c) /= u(r, r);
        }
    }

    return prhs;
}

template <class SquareMatrix>
typename LUP_Decomp<SquareMatrix>::Real LUP_Decomp<SquareMatrix>::l(int r, int c) const {
    if (r == c)
        return Real(1);
    if (r < c)
        return Real(0);
    return l_and_u(r, c);
}

template <class SquareMatrix>
typename LUP_Decomp<SquareMatrix>::Real LUP_Decomp<SquareMatrix>::u(int r, int c) const {
    if (r > c)
        return Real(0);
    return l_and_u(r, c);
}

template <class SquareMatrix>
typename LUP_Decomp<SquareMatrix>::Real &LUP_Decomp<SquareMatrix>::l_ref(int r, int c) {
    if (r <= c)
        throw std::out_of_range("Diagonal and upper triangle of L are fixed at index " + std::to_string(r) + "," + std::to_string(c));
    return l_and_u(r, c);
}

template <class SquareMatrix>
typename LUP_Decomp<SquareMatrix>::Real &LUP_Decomp<SquareMatrix>::u_ref(int r, int c) {
    if (r > c)
        throw std::out_of_range("Lower triangle of U is fixed at index " + std::to_string(r) + "," + std::to_string(c));
    return l_and_u(r, c);
}

}  // namespace linalg
