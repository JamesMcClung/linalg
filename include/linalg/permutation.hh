#pragma once

#include "matrix.hh"

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DESCRIPTION                             //
////////////////////////////////////////////////////////////////////////

// A PermutedMatrix wraps another matrix and proves a permuted view
//   of that matrix's rows or columns.

// A PermutationMatrix is a matrix used to create PermutedMatrices via
//   multiplication.

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

    template <int _dim, typename _Real, class _SomeMatrix>
    friend auto operator*(const PermutationMatrix<_dim, _Real>& p, _SomeMatrix&& m);
    template <int _dim, typename _Real, class _SomeMatrix>
    friend auto operator*(_SomeMatrix&& m, const PermutationMatrix<_dim, _Real>& p);
};

////////////////////////////////////////////////////////////////////////
// PermutedMatrix

template <class SomeMatrix, typename = void>
class PermutedMatrix;

template <class SomeMatrix>
class PermutedMatrix<SomeMatrix, std::enable_if_t<!std::is_const_v<SomeMatrix>>> : public Matrix<SomeMatrix::nrows, SomeMatrix::ncols, typename SomeMatrix::real_t> {
   private:
    SomeMatrix& original_matrix;
    int row_lookups[SomeMatrix::nrows];
    int col_lookups[SomeMatrix::ncols];

    PermutedMatrix(SomeMatrix& m) : original_matrix(m) {}

   public:
    typename SomeMatrix::real_t operator()(int i, int j) const {
        return original_matrix(row_lookups[i], col_lookups[j]);
    }
    typename SomeMatrix::real_t& operator()(int i, int j) {
        return original_matrix(row_lookups[i], col_lookups[j]);
    }

    template <int _dim, typename _Real, class _SomeMatrix>
    friend auto operator*(const PermutationMatrix<_dim, _Real>& p, _SomeMatrix&& m);
    template <int _dim, typename _Real, class _SomeMatrix>
    friend auto operator*(_SomeMatrix&& m, const PermutationMatrix<_dim, _Real>& p);
};

template <class SomeMatrix>
class PermutedMatrix<SomeMatrix, std::enable_if_t<std::is_const_v<SomeMatrix>>> : public Matrix<SomeMatrix::nrows, SomeMatrix::ncols, typename SomeMatrix::real_t> {
   private:
    const SomeMatrix& original_matrix;
    int row_lookups[SomeMatrix::nrows];
    int col_lookups[SomeMatrix::ncols];

    PermutedMatrix(const SomeMatrix& m) : original_matrix(m) {}

   public:
    typename SomeMatrix::real_t operator()(int i, int j) const {
        return original_matrix(row_lookups[i], col_lookups[j]);
    }

    template <int _dim, typename _Real, class _SomeMatrix>
    friend auto operator*(const PermutationMatrix<_dim, _Real>& p, _SomeMatrix&& m);
    template <int _dim, typename _Real, class _SomeMatrix>
    friend auto operator*(_SomeMatrix&& m, const PermutationMatrix<_dim, _Real>& p);
};

template <int dim, typename Real, class SomeMatrix>
auto operator*(const PermutationMatrix<dim, Real>& p, SomeMatrix&& m);

template <int dim, typename Real, class SomeMatrix>
auto operator*(SomeMatrix&& m, const PermutationMatrix<dim, Real>& p);

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

////////////////////////////////////////////////////////////////////////
// PermutedMatrix

template <int dim, typename Real, class SomeMatrix>
auto operator*(const PermutationMatrix<dim, Real>& p, SomeMatrix&& m) {
    auto pm = PermutedMatrix<std::remove_reference_t<SomeMatrix>>(std::forward<SomeMatrix>(m));
    for (int i = 0; i < pm.nrows; i++) {
        pm.row_lookups[i] = p.row_lookups[i];
    }
    for (int i = 0; i < pm.ncols; i++) {
        pm.col_lookups[i] = i;
    }
    return pm;
}

template <int dim, typename Real, class SomeMatrix>
auto operator*(SomeMatrix&& m, const PermutationMatrix<dim, Real>& p) {
    auto pm = PermutedMatrix<std::remove_reference_t<SomeMatrix>>(std::forward<SomeMatrix>(m));
    for (int i = 0; i < pm.nrows; i++) {
        pm.row_lookups[i] = i;
    }
    for (int i = 0; i < pm.ncols; i++) {
        pm.col_lookups[i] = p.col_lookups[i];
    }
    return pm;
}

}  // namespace linalg
