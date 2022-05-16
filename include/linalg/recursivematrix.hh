#pragma once

#include <memory>
#include <type_traits>

#include "zeromatrix.hh"

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DESCRIPTION                             //
////////////////////////////////////////////////////////////////////////

// A RecursiveMatrix is essentially a matrix of matrices.
// Viewed as a normal matrix, its dimensions are nrows-by-ncols.
// As a matrix of matrices, its dimensions are nbrows-by-nbcols.
// For brevity, nrows, ncols, nbrows, and nbcols are abbreviated
//   nr, nc, nbr, and nbc, respectively.

////////////////////////////////////////////////////////////////////////
//                              UTILITY                               //
////////////////////////////////////////////////////////////////////////

namespace _recursivematrix_util {

constexpr int max(int a, int b) {
    return a < b ? b : a;
}

constexpr int min(int a, int b) {
    return a < b ? a : b;
}

}  // namespace _recursivematrix_util

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

// A submatrix with either dimension larger than this is not allowed to
//   be a FullMatrix.
#define MAX_FMAT_DIM 256

template <int nrows, int ncols, int nbrows, int nbcols, typename Real>
class RecursiveMatrix : public Matrix<nrows, ncols, Real> {
   private:
    static constexpr int submat_nrows = _recursivematrix_util::max(nrows / nbrows, 1);
    static constexpr int submat_ncols = _recursivematrix_util::max(ncols / nbcols, 1);
    static constexpr int submat_nbrows = submat_nrows <= MAX_FMAT_DIM ? 1 : _recursivematrix_util::min(nbrows, submat_nrows);
    static constexpr int submat_nbcols = submat_ncols <= MAX_FMAT_DIM ? 1 : _recursivematrix_util::min(nbcols, submat_ncols);
    static constexpr bool has_full_submats = submat_nrows <= MAX_FMAT_DIM && submat_ncols <= MAX_FMAT_DIM;

    using Submatrix = typename std::conditional<has_full_submats,
                                                FullMatrix<submat_nrows, submat_ncols, Real>,
                                                RecursiveMatrix<submat_nrows, submat_ncols, submat_nbrows, submat_nbcols, Real>>::type;

    std::unique_ptr<Submatrix> submat_ptrs[nbrows][nbcols];

    bool is_submat_zero(const std::unique_ptr<Submatrix> &p) const;

   public:
    RecursiveMatrix() = default;
    RecursiveMatrix(const Real &r);
    explicit RecursiveMatrix(const Matrix<nrows, ncols, Real> &m);
    explicit RecursiveMatrix(const ZeroMatrix<nrows, ncols, Real> &z);
    RecursiveMatrix(const RecursiveMatrix &r);
    RecursiveMatrix(RecursiveMatrix &&r);

    Real operator()(int i, int j) const override;
    Real &operator()(int i, int j);

    bool operator==(const RecursiveMatrix &m2) const;
    bool operator!=(const RecursiveMatrix &m2) const;

    bool operator==(const ZeroMatrix<nrows, ncols, Real> &m2) const;
    bool operator!=(const ZeroMatrix<nrows, ncols, Real> &m2) const;

    RecursiveMatrix operator+(const RecursiveMatrix &m) const;
    RecursiveMatrix operator-(const RecursiveMatrix &m) const;

    RecursiveMatrix operator-() const;

    RecursiveMatrix operator*(const Real &r) const;
    RecursiveMatrix operator/(const Real &r) const;

    RecursiveMatrix &operator=(const RecursiveMatrix &m);
    RecursiveMatrix &operator=(RecursiveMatrix &&m);

    RecursiveMatrix &operator+=(const RecursiveMatrix &m);
    RecursiveMatrix &operator-=(const RecursiveMatrix &m);

    RecursiveMatrix &operator*=(const Real &r);
    RecursiveMatrix &operator/=(const Real &r);

    template <int nr, int ni, int nc, int nbr, int nbi, int nbc, typename R>
    friend RecursiveMatrix<nr, nc, nbr, nbc, R> operator*(const RecursiveMatrix<nr, ni, nbr, nbi, R> &m1, const RecursiveMatrix<ni, nc, nbi, nbc, R> &m2);
};

template <int nr, int nc, int nbr, int nbc, typename Real>
bool operator==(const ZeroMatrix<nr, nc, Real> &m1, const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m2);
template <int nr, int nc, int nbr, int nbc, typename Real>
bool operator!=(const ZeroMatrix<nr, nc, Real> &m1, const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m2);

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> operator*(const Real &r, const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m);

template <int nr, int ni, int nc, int nbr, int nbi, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> operator*(const RecursiveMatrix<nr, ni, nbr, nbi, Real> &m1, const RecursiveMatrix<ni, nc, nbi, nbc, Real> &m2);

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

template <int nr, int nc, int nbr, int nbc, typename Real>
bool RecursiveMatrix<nr, nc, nbr, nbc, Real>::is_submat_zero(const std::unique_ptr<Submatrix> &p) const {
    return p == nullptr || *p == ZeroMatrix<submat_nrows, submat_ncols, Real>();
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real>::RecursiveMatrix(const Real &r) : RecursiveMatrix() {
    if (r != Real(0)) {
        for (int br = 0; br < nbr; br++) {
            for (int bc = 0; bc < nbc; bc++) {
                submat_ptrs[br][bc].reset(new Submatrix(r));
            }
        }
    }
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real>::RecursiveMatrix(const Matrix<nr, nc, Real> &m) : RecursiveMatrix() {
    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nc; j++) {
            (*this)(i, j) = m.get(i, j);
        }
    }
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real>::RecursiveMatrix(const ZeroMatrix<nr, nc, Real> &) : RecursiveMatrix() {}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real>::RecursiveMatrix(const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m) : RecursiveMatrix() {
    *this = m;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real>::RecursiveMatrix(RecursiveMatrix<nr, nc, nbr, nbc, Real> &&m) : RecursiveMatrix() {
    *this = std::move(m);
}

template <int nr, int nc, int nrb, int ncb, typename Real>
Real RecursiveMatrix<nr, nc, nrb, ncb, Real>::operator()(int i, int j) const {
    auto &submat_ptr = submat_ptrs[i / submat_nrows][j / submat_ncols];

    if (submat_ptr == nullptr) {
        return Real(0);
    }
    return submat_ptr->get(i % submat_nrows, j % submat_ncols);
}

template <int nr, int nc, int nbr, int nbc, typename Real>
Real &RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator()(int i, int j) {
    int ib = i / submat_nrows;
    int jb = j / submat_ncols;

    auto &submat_ptr = submat_ptrs[ib][jb];

    if (submat_ptr == nullptr) {
        submat_ptr.reset(new Submatrix(Real(0)));
    }
    return (*submat_ptr)(i % submat_nrows, j % submat_ncols);
}

template <int nr, int nc, int nbr, int nbc, typename Real>
bool RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator==(const ZeroMatrix<nr, nc, Real> &) const {
    ZeroMatrix<submat_nrows, submat_ncols, Real> zmat;
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &submat_ptr = submat_ptrs[i][j];
            if (!(submat_ptr == nullptr || *submat_ptr == zmat)) {
                return false;
            }
        }
    }
    return true;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
bool RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator!=(const ZeroMatrix<nr, nc, Real> &m) const {
    return !(*this == m);
}

template <int nr, int nc, int nbr, int nbc, typename Real>
bool operator==(const ZeroMatrix<nr, nc, Real> &m1, const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m2) {
    return m2 == m1;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
bool operator!=(const ZeroMatrix<nr, nc, Real> &m1, const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m2) {
    return !(m2 == m1);
}

template <int nr, int nc, int nbr, int nbc, typename Real>
bool RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator==(const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m) const {
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &p1 = submat_ptrs[i][j];
            auto &p2 = m.submat_ptrs[i][j];
            if (is_submat_zero(p1) != is_submat_zero(p2)) {
                return false;
            } else if (p1 != nullptr && p2 != nullptr && *p1 != *p2) {
                return false;
            }
        }
    }
    return true;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
bool RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator!=(const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m) const {
    return !(*this == m);
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator+(const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m) const {
    RecursiveMatrix<nr, nc, nbr, nbc, Real> ret;
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &p1 = submat_ptrs[i][j];
            auto &p2 = m.submat_ptrs[i][j];
            bool isZero1 = is_submat_zero(p1);
            bool isZero2 = is_submat_zero(p2);
            if (isZero1 && !isZero2) {
                ret.submat_ptrs[i][j].reset(new Submatrix(*p2));
            } else if (!isZero1 && isZero2) {
                ret.submat_ptrs[i][j].reset(new Submatrix(*p1));
            } else if (!isZero1 && !isZero2) {
                ret.submat_ptrs[i][j].reset(new Submatrix(*p1 + *p2));
            }
        }
    }
    return ret;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator-(const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m) const {
    RecursiveMatrix<nr, nc, nbr, nbc, Real> ret;
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &p1 = submat_ptrs[i][j];
            auto &p2 = m.submat_ptrs[i][j];
            bool isZero1 = is_submat_zero(p1);
            bool isZero2 = is_submat_zero(p2);
            if (isZero1 && !isZero2) {
                ret.submat_ptrs[i][j].reset(new Submatrix(-*p2));
            } else if (!isZero1 && isZero2) {
                ret.submat_ptrs[i][j].reset(new Submatrix(*p1));
            } else if (!isZero1 && !isZero2) {
                ret.submat_ptrs[i][j].reset(new Submatrix(*p1 - *p2));
            }
        }
    }
    return ret;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator-() const {
    RecursiveMatrix<nr, nc, nbr, nbc, Real> ret;
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &p = submat_ptrs[i][j];
            if (!is_submat_zero(p)) {
                ret.submat_ptrs[i][j].reset(new Submatrix(-*p));
            }
        }
    }
    return ret;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator*(const Real &r) const {
    RecursiveMatrix<nr, nc, nbr, nbc, Real> ret;
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &p = submat_ptrs[i][j];
            if (!is_submat_zero(p)) {
                ret.submat_ptrs[i][j].reset(new Submatrix((*p) * r));
            }
        }
    }
    return ret;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> operator*(const Real &r, const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m) {
    return m * r;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator/(const Real &r) const {
    RecursiveMatrix<nr, nc, nbr, nbc, Real> ret;
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &p = submat_ptrs[i][j];
            if (!is_submat_zero(p)) {
                ret.submat_ptrs[i][j].reset(new Submatrix((*p) / r));
            }
        }
    }
    return ret;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> &RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator=(const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m) {
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &p = m.submat_ptrs[i][j];
            if (!is_submat_zero(p)) {
                submat_ptrs[i][j].reset(new Submatrix(*p));
            }
        }
    }
    return *this;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> &RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator=(RecursiveMatrix<nr, nc, nbr, nbc, Real> &&m) {
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            submat_ptrs[i][j] = std::move(m.submat_ptrs[i][j]);
        }
    }
    return *this;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> &RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator+=(const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m) {
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &p1 = submat_ptrs[i][j];
            auto &p2 = m.submat_ptrs[i][j];
            if (!is_submat_zero(p2)) {
                if (is_submat_zero(p1)) {
                    submat_ptrs[i][j].reset(new Submatrix(*p2));
                } else {
                    *submat_ptrs[i][j] += *p2;
                }
            }
        }
    }
    return *this;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> &RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator-=(const RecursiveMatrix<nr, nc, nbr, nbc, Real> &m) {
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &p1 = submat_ptrs[i][j];
            auto &p2 = m.submat_ptrs[i][j];
            if (!is_submat_zero(p2)) {
                if (is_submat_zero(p1)) {
                    Submatrix *newp1 = new Submatrix(Real(0));
                    *newp1 -= *p2;
                    submat_ptrs[i][j].reset(newp1);
                } else {
                    *submat_ptrs[i][j] -= *p2;
                }
            }
        }
    }
    return *this;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> &RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator*=(const Real &r) {
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &p = submat_ptrs[i][j];
            if (!is_submat_zero(p)) {
                *p *= r;
            }
        }
    }
    return *this;
}

template <int nr, int nc, int nbr, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> &RecursiveMatrix<nr, nc, nbr, nbc, Real>::operator/=(const Real &r) {
    for (int i = 0; i < nbr; i++) {
        for (int j = 0; j < nbc; j++) {
            auto &p = submat_ptrs[i][j];
            if (!is_submat_zero(p)) {
                *p /= r;
            }
        }
    }
    return *this;
}

template <int nr, int ni, int nc, int nbr, int nbi, int nbc, typename Real>
RecursiveMatrix<nr, nc, nbr, nbc, Real> operator*(const RecursiveMatrix<nr, ni, nbr, nbi, Real> &m1, const RecursiveMatrix<ni, nc, nbi, nbc, Real> &m2) {
    RecursiveMatrix<nr, nc, nbr, nbc, Real> ret;
    for (int br = 0; br < nbr; br++) {
        for (int bc = 0; bc < nbc; bc++) {
            auto &p = ret.submat_ptrs[br][bc];
            for (int bi = 0; bi < nbi; bi++) {
                auto &p1 = m1.submat_ptrs[br][bi];
                auto &p2 = m2.submat_ptrs[bi][bc];
                if (!m1.is_submat_zero(p1) && !m2.is_submat_zero(p2)) {
                    if (p == nullptr) {
                        p.reset(new typename decltype(ret)::Submatrix(*p1 * *p2));
                    } else {
                        *p += *p1 * *p2;
                    }
                }
            }
        }
    }
    return ret;
}

}  // namespace linalg