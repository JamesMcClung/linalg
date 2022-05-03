#pragma once

#include "matrix.hh"

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DESCRIPTION                             //
////////////////////////////////////////////////////////////////////////

// A PermutationMatrix is a permutation matrix.

////////////////////////////////////////////////////////////////////////
//                              UTILITY                               //
////////////////////////////////////////////////////////////////////////

namespace _permutation_util {

}  // namespace _permutation_util

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// PermutationMatrix

template <int dim, typename Real>
class PermutationMatrix : public Matrix<dim, dim, Real> {
   private:
    int row_lookups[dim];
    int col_lookups[dim];

   public:
    // default constructor yields identity matrix
    PermutationMatrix();

    void swapRows(int r1, int r2);
    void swapCols(int c1, int c2);

    Real operator()(int i, int j) const;

    template <int _dim, typename _Real, class _SomeMatrix, typename>
    friend _SomeMatrix operator*(const PermutationMatrix<_dim, _Real>& p, _SomeMatrix m);
    template <int _dim, typename _Real, class _SomeMatrix, typename>
    friend _SomeMatrix operator*(_SomeMatrix m, const PermutationMatrix<_dim, _Real>& p);
};

template <int dim, typename Real, class SomeMatrix, typename = std::enable_if_t<SomeMatrix::nrows == dim>>
SomeMatrix operator*(const PermutationMatrix<dim, Real>& p, SomeMatrix m);

template <int dim, typename Real, class SomeMatrix, typename = std::enable_if_t<SomeMatrix::ncols == dim>>
SomeMatrix operator*(SomeMatrix m, const PermutationMatrix<dim, Real>& p);

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// PermutationMatrix

template <int dim, typename Real>
PermutationMatrix<dim, Real>::PermutationMatrix() {
    for (int i = 0; i < dim; i++) {
        row_lookups[i] = i;
        col_lookups[i] = i;
    }
}

template <int dim, typename Real>
void PermutationMatrix<dim, Real>::swapRows(int r1, int r2) {
    std::swap(row_lookups[r1], row_lookups[r2]);
}

template <int dim, typename Real>
void PermutationMatrix<dim, Real>::swapCols(int c1, int c2) {
    std::swap(col_lookups[c1], col_lookups[c2]);
}

template <int dim, typename Real>
Real PermutationMatrix<dim, Real>::operator()(int i, int j) const {
    return Real(row_lookups[i] == col_lookups[j]);
}

template <int dim, typename Real, class SomeMatrix, typename>
SomeMatrix operator*(const PermutationMatrix<dim, Real>& p, SomeMatrix m) {
    bool wasRowHandled[dim];
    // permute rows of m
    for (int c = 0; c < SomeMatrix::ncols; c++) {
        for (int r = 0; r < dim; r++) {
            if (!wasRowHandled[r]) {
                Real tmp = m(r, c);
                int thisr = r;
                for (int nextr = p.row_lookups[r]; nextr != r; nextr = p.row_lookups[nextr]) {
                    m(thisr, c) = m(nextr, c);
                    thisr = nextr;
                    wasRowHandled[nextr] = true;
                }
                m(thisr, c) = tmp;
            }
        }
    }
    return m;
}

template <int dim, typename Real, class SomeMatrix, typename>
SomeMatrix operator*(SomeMatrix m, const PermutationMatrix<dim, Real>& p) {
    bool wasColHandled[dim];
    // permute cols of m
    for (int r = 0; r < SomeMatrix::nrows; r++) {
        for (int c = 0; c < dim; c++) {
            if (!wasColHandled[c]) {
                Real tmp = m(r, c);
                int thisc = c;
                for (int nextc = p.col_lookups[c]; nextc != c; nextc = p.col_lookups[nextc]) {
                    m(r, thisc) = m(r, nextc);
                    thisc = nextc;
                    wasColHandled[nextc] = true;
                }
                m(r, thisc) = tmp;
            }
        }
    }
    return m;
}

}  // namespace linalg
