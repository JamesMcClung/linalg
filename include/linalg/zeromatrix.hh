#pragma once

#include <type_traits>

#include "fullmatrix.hh"

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

template <int nrows, int ncols, typename Real>
class ZeroMatrix : public Matrix<nrows, ncols, Real> {
   public:
    Real operator()(int i, int j) const override;

    ZeroMatrix<nrows, ncols, Real> &operator*=(const Real &r);
    ZeroMatrix<nrows, ncols, Real> &operator/=(const Real &r);
};

template <int nrows, int ninner, int ncols, typename Real>
bool operator==(const ZeroMatrix<nrows, ninner, Real> &m1, const ZeroMatrix<ninner, ncols, Real> &m2);
template <int nrows, int ninner, int ncols, typename Real>
bool operator!=(const ZeroMatrix<nrows, ninner, Real> &m1, const ZeroMatrix<ninner, ncols, Real> &m2);

template <int nrows, int ninner, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator*(const Matrix<nrows, ninner, Real> &m1, const ZeroMatrix<ninner, ncols, Real> &m2);
template <int nrows, int ninner, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator*(const ZeroMatrix<nrows, ninner, Real> &m1, const Matrix<ninner, ncols, Real> &m2);
template <int nrows, int ninner, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator*(const ZeroMatrix<nrows, ninner, Real> &m1, const ZeroMatrix<ninner, ncols, Real> &m2);

template <int nrows, int ncols, typename Real, class AnyMatrix>
std::enable_if_t<std::is_base_of_v<Matrix<nrows, ncols, Real>, AnyMatrix>, AnyMatrix> operator+(const ZeroMatrix<nrows, ncols, Real> &m1, const AnyMatrix &m2);
template <int nrows, int ncols, typename Real, class AnyMatrix>
std::enable_if_t<std::is_base_of_v<Matrix<nrows, ncols, Real>, AnyMatrix>, AnyMatrix> operator+(const AnyMatrix &m1, const ZeroMatrix<nrows, ncols, Real> &m2);
template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator+(const ZeroMatrix<nrows, ncols, Real> &m1, const ZeroMatrix<nrows, ncols, Real> &m2);

template <int nrows, int ncols, typename Real, class AnyMatrix>
std::enable_if_t<std::is_base_of_v<Matrix<nrows, ncols, Real>, AnyMatrix>, AnyMatrix> operator-(const ZeroMatrix<nrows, ncols, Real> &m1, const AnyMatrix &m2);
template <int nrows, int ncols, typename Real, class AnyMatrix>
std::enable_if_t<std::is_base_of_v<Matrix<nrows, ncols, Real>, AnyMatrix>, AnyMatrix> operator-(const AnyMatrix &m1, const ZeroMatrix<nrows, ncols, Real> &m2);
template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator-(const ZeroMatrix<nrows, ncols, Real> &m1, const ZeroMatrix<nrows, ncols, Real> &m2);

template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator-(const ZeroMatrix<nrows, ncols, Real> &m);

template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator*(const ZeroMatrix<nrows, ncols, Real> &m, const Real &r);
template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator*(const Real &r, const ZeroMatrix<nrows, ncols, Real> &m);

template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator/(const ZeroMatrix<nrows, ncols, Real> &m, const Real &r);

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

template <int nrows, int ncols, typename Real>
Real ZeroMatrix<nrows, ncols, Real>::operator()(int i, int j) const {
    return Real(0);
}

template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> &ZeroMatrix<nrows, ncols, Real>::operator*=(const Real &r) {
    return *this;
}
template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> &ZeroMatrix<nrows, ncols, Real>::operator/=(const Real &r) {
    return *this;
}

template <int nrows, int ninner, int ncols, typename Real>
bool operator==(const ZeroMatrix<nrows, ninner, Real> &m1, const ZeroMatrix<ninner, ncols, Real> &m2) {
    return true;
}
template <int nrows, int ninner, int ncols, typename Real>
bool operator!=(const ZeroMatrix<nrows, ninner, Real> &m1, const ZeroMatrix<ninner, ncols, Real> &m2) {
    return false;
}

template <int nrows, int ninner, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator*(const Matrix<nrows, ninner, Real> &m1, const ZeroMatrix<ninner, ncols, Real> &m2) {
    return ZeroMatrix<nrows, ncols, Real>();
}
template <int nrows, int ninner, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator*(const ZeroMatrix<nrows, ninner, Real> &m1, const Matrix<ninner, ncols, Real> &m2) {
    return ZeroMatrix<nrows, ncols, Real>();
}
template <int nrows, int ninner, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator*(const ZeroMatrix<nrows, ninner, Real> &m1, const ZeroMatrix<ninner, ncols, Real> &m2) {
    return ZeroMatrix<nrows, ncols, Real>();
}

template <int nrows, int ncols, typename Real, class AnyMatrix>
std::enable_if_t<std::is_base_of_v<Matrix<nrows, ncols, Real>, AnyMatrix>, AnyMatrix> operator+(const ZeroMatrix<nrows, ncols, Real> &m1, const AnyMatrix &m2) {
    return AnyMatrix(m2);
}
template <int nrows, int ncols, typename Real, class AnyMatrix>
std::enable_if_t<std::is_base_of_v<Matrix<nrows, ncols, Real>, AnyMatrix>, AnyMatrix> operator+(const AnyMatrix &m1, const ZeroMatrix<nrows, ncols, Real> &m2) {
    return AnyMatrix(m1);
}
template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator+(const ZeroMatrix<nrows, ncols, Real> &m1, const ZeroMatrix<nrows, ncols, Real> &m2) {
    return ZeroMatrix<nrows, ncols, Real>();
}

template <int nrows, int ncols, typename Real, class AnyMatrix>
std::enable_if_t<std::is_base_of_v<Matrix<nrows, ncols, Real>, AnyMatrix>, AnyMatrix> operator-(const ZeroMatrix<nrows, ncols, Real> &m1, const AnyMatrix &m2) {
    return -m2;
}
template <int nrows, int ncols, typename Real, class AnyMatrix>
std::enable_if_t<std::is_base_of_v<Matrix<nrows, ncols, Real>, AnyMatrix>, AnyMatrix> operator-(const AnyMatrix &m1, const ZeroMatrix<nrows, ncols, Real> &m2) {
    return AnyMatrix(m1);
}
template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator-(const ZeroMatrix<nrows, ncols, Real> &m1, const ZeroMatrix<nrows, ncols, Real> &m2) {
    return ZeroMatrix<nrows, ncols, Real>();
}

template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator-(const ZeroMatrix<nrows, ncols, Real> &m) {
    return ZeroMatrix<nrows, ncols, Real>();
}

template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator*(const ZeroMatrix<nrows, ncols, Real> &m, const Real &r) {
    return m;
}
template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator*(const Real &r, const ZeroMatrix<nrows, ncols, Real> &m) {
    return m;
}

template <int nrows, int ncols, typename Real>
ZeroMatrix<nrows, ncols, Real> operator/(const ZeroMatrix<nrows, ncols, Real> &m, const Real &r) {
    return m;
}

}  // namespace linalg