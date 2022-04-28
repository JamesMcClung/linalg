#pragma once

#include "matrix.hh"

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

template <int nrows, int ncols, typename Real>
class FullMatrix : public Matrix<nrows, ncols, Real> {
   private:
    Real vals[nrows][ncols];

   public:
    FullMatrix() = default;
    FullMatrix(const Real &r);
    FullMatrix(const Real (&_vals)[nrows][ncols]);
    explicit FullMatrix(const Matrix<nrows, ncols, Real> &m);

    ////////////////////////////////////////////////////////////////////
    // Indexing

    Real &operator()(int i, int j) {
        return vals[i][j];
    }
    Real operator()(int i, int j) const override {
        return vals[i][j];
    }

    ////////////////////////////////////////////////////////////////////
    // Assignment operators

    FullMatrix &operator*=(const Real &r);
    FullMatrix &operator/=(const Real &r);

    FullMatrix &operator+=(const Matrix<nrows, ncols, Real> &m);
    FullMatrix &operator-=(const Matrix<nrows, ncols, Real> &m);

    ////////////////////////////////////////////////////////////////////
    // Arithmetic operators

    FullMatrix operator*(const Real &r) const {
        return FullMatrix(*this) *= r;
    }
    friend FullMatrix operator*(const Real &r, const FullMatrix &m) {
        return m * r;
    }

    FullMatrix operator/(const Real &r) const {
        return FullMatrix(*this) /= r;
    }

    FullMatrix operator-() const {
        return FullMatrix(Real(0)) -= *this;
    }

    FullMatrix operator-(const FullMatrix &m) const {
        return FullMatrix(*this) -= m;
    }

    FullMatrix operator+(FullMatrix m) const {
        return m += *this;
    }
};

template <int nrows, int ninner, int ncols, typename Real>
FullMatrix<nrows, ncols, Real> operator*(const FullMatrix<nrows, ninner, Real> &m1, const FullMatrix<ninner, ncols, Real> &m2);

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

template <int nrows, int ncols, typename Real>
FullMatrix<nrows, ncols, Real>::FullMatrix(const Real &r) {
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            vals[i][j] = r;
        }
    }
}

template <int nrows, int ncols, typename Real>
FullMatrix<nrows, ncols, Real>::FullMatrix(const Real (&_vals)[nrows][ncols]) {
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            vals[i][j] = _vals[i][j];
        }
    }
}

template <int nrows, int ncols, typename Real>
FullMatrix<nrows, ncols, Real>::FullMatrix(const Matrix<nrows, ncols, Real> &m) {
    for (int r = 0; r < nrows; r++) {
        for (int c = 0; c < ncols; c++) {
            vals[r][c] = m(r, c);
        }
    }
}

template <int nrows, int ncols, typename Real>
FullMatrix<nrows, ncols, Real> &FullMatrix<nrows, ncols, Real>::operator*=(const Real &r) {
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            vals[i][j] *= r;
        }
    }
    return *this;
}

template <int nrows, int ncols, typename Real>
FullMatrix<nrows, ncols, Real> &FullMatrix<nrows, ncols, Real>::operator/=(const Real &r) {
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            vals[i][j] /= r;
        }
    }
    return *this;
}

template <int nrows, int ncols, typename Real>
FullMatrix<nrows, ncols, Real> &FullMatrix<nrows, ncols, Real>::operator+=(const Matrix<nrows, ncols, Real> &m) {
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            vals[i][j] += m(i, j);
        }
    }
    return *this;
}

template <int nrows, int ncols, typename Real>
FullMatrix<nrows, ncols, Real> &FullMatrix<nrows, ncols, Real>::operator-=(const Matrix<nrows, ncols, Real> &m) {
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            vals[i][j] -= m(i, j);
        }
    }
    return *this;
}

template <int nrows, int ninner, int ncols, typename Real>
FullMatrix<nrows, ncols, Real> operator*(const FullMatrix<nrows, ninner, Real> &m1, const FullMatrix<ninner, ncols, Real> &m2) {
    FullMatrix<nrows, ncols, Real> m;
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            Real mij = m1(i, 0) * m2(0, j);
            for (int k = 1; k < ninner; k++) {
                mij += m1(i, k) * m2(k, j);
            }
            m(i, j) = mij;
        }
    }
    return m;
}

}  // namespace linalg