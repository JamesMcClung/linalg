#pragma once

#include <memory>
#include <string>

#include "fullmatrix.hh"

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DESCRIPTION                             //
////////////////////////////////////////////////////////////////////////

// A BandedMatrix stores one or more bands of a matrix. The diagonal is
//   an example of a band and is always stored.
//
// The lower bandwidth (lbw) is the number of bands below the diagonal.
// The upper bandwidth (ubw) is the number of bands above the diagonal.
// An element at index (r,c) is stored iff -lbw <= c-r <= ubw.

////////////////////////////////////////////////////////////////////////
//                              UTILITY                               //
////////////////////////////////////////////////////////////////////////

namespace _bandedmatrix_util {

constexpr int max(int a, int b) {
    return a < b ? b : a;
}

constexpr int min(int a, int b) {
    return a < b ? a : b;
}

constexpr int modPositive(int a, int mod) {
    int res = a % mod;
    if (a < 0)
        a += mod;
    return a;
}

const int index_not_on_stored_band = -1;

}  // namespace _bandedmatrix_util

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

template <int dim, int lbw, int ubw, typename Real>
class BandedMatrix : public Matrix<dim, dim, Real> {
   private:
    static constexpr int nbands = _bandedmatrix_util::min(1 + lbw + ubw, dim);
    using vals_t = FullMatrix<dim, nbands, Real>;

    std::unique_ptr<vals_t> vals_ptr;

    static int getValsCol(int r, int c);

   public:
    static constexpr int lowerBW = lbw;
    static constexpr int upperBW = ubw;

    BandedMatrix();
    BandedMatrix(const Real &r);
    BandedMatrix(const BandedMatrix &bm);
    BandedMatrix(BandedMatrix &&bm);

    Real &operator()(int r, int c);
    Real operator()(int r, int c) const override;

    BandedMatrix &operator*=(const Real &r);
    BandedMatrix &operator/=(const Real &r);

    BandedMatrix &operator+=(const BandedMatrix &bm);
    BandedMatrix &operator-=(const BandedMatrix &bm);

    template <int nc>
    FullMatrix<dim, nc, Real> operator*(const FullMatrix<dim, nc, Real> &fm) const;
};

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

template <int dim, int lbw, int ubw, typename Real>
int BandedMatrix<dim, lbw, ubw, Real>::getValsCol(int r, int c) {
    int band = c - r;
    if (0 <= band && band <= ubw && r < dim - band) {
        return band;
    } else if (-lbw <= band && band < 0 && r >= -band) {
        return band + nbands;
    }
    return _bandedmatrix_util::index_not_on_stored_band;
}

template <int dim, int lbw, int ubw, typename Real>
BandedMatrix<dim, lbw, ubw, Real>::BandedMatrix() {
    vals_ptr.reset(new vals_t);
}

template <int dim, int lbw, int ubw, typename Real>
BandedMatrix<dim, lbw, ubw, Real>::BandedMatrix(const Real &r) {
    vals_ptr.reset(new vals_t(r));
}

template <int dim, int lbw, int ubw, typename Real>
BandedMatrix<dim, lbw, ubw, Real>::BandedMatrix(const BandedMatrix &bm) {
    vals_ptr.reset(new vals_t(*bm.vals_ptr));
}

template <int dim, int lbw, int ubw, typename Real>
BandedMatrix<dim, lbw, ubw, Real>::BandedMatrix(BandedMatrix &&bm) : BandedMatrix() {
    vals_ptr.swap(bm.vals_ptr);
}

template <int dim, int lbw, int ubw, typename Real>
Real BandedMatrix<dim, lbw, ubw, Real>::operator()(int r, int c) const {
    if (r < 0 || r >= dim || c < 0 || c >= dim) {
        throw std::out_of_range("Index out of bounds: (" + std::to_string(r) + ", " + std::to_string(c) + ") for dimensions " + std::to_string(dim) + "x" + std::to_string(dim));
    }

    int vr = r;
    int vc = getValsCol(r, c);

    if (vc == _bandedmatrix_util::index_not_on_stored_band) {
        return Real(0);
    }

    return vals_ptr->get(vr, vc);
}

template <int dim, int lbw, int ubw, typename Real>
Real &BandedMatrix<dim, lbw, ubw, Real>::operator()(int r, int c) {
    if (r < 0 || r >= dim || c < 0 || c >= dim) {
        throw std::out_of_range("Index out of bounds: (" + std::to_string(r) + ", " + std::to_string(c) + ") for dimensions " + std::to_string(dim) + "x" + std::to_string(dim));
    }

    int vr = r;
    int vc = getValsCol(r, c);

    if (vc == _bandedmatrix_util::index_not_on_stored_band) {
        throw std::out_of_range("Index not on band: (" + std::to_string(r) + ", " + std::to_string(c) + ")");
    }

    return (*vals_ptr)(vr, vc);
}

template <int dim, int lbw, int ubw, typename Real>
BandedMatrix<dim, lbw, ubw, Real> &BandedMatrix<dim, lbw, ubw, Real>::operator*=(const Real &r) {
    *vals_ptr *= r;
    return *this;
}

template <int dim, int lbw, int ubw, typename Real>
BandedMatrix<dim, lbw, ubw, Real> &BandedMatrix<dim, lbw, ubw, Real>::operator/=(const Real &r) {
    *vals_ptr /= r;
    return *this;
}

template <int dim, int lbw, int ubw, typename Real>
BandedMatrix<dim, lbw, ubw, Real> &BandedMatrix<dim, lbw, ubw, Real>::operator+=(const BandedMatrix<dim, lbw, ubw, Real> &bm) {
    *vals_ptr += *bm.vals_ptr;
    return *this;
}

template <int dim, int lbw, int ubw, typename Real>
BandedMatrix<dim, lbw, ubw, Real> &BandedMatrix<dim, lbw, ubw, Real>::operator-=(const BandedMatrix<dim, lbw, ubw, Real> &bm) {
    *vals_ptr -= *bm.vals_ptr;
    return *this;
}

template <int dim, int lbw, int ubw, typename Real>
template <int nc>
FullMatrix<dim, nc, Real> BandedMatrix<dim, lbw, ubw, Real>::operator*(const FullMatrix<dim, nc, Real> &fm) const {
    FullMatrix<dim, nc, Real> res(Real(0));
    for (int r = 0; r < res.nrows; r++) {
        for (int c = 0; c < res.ncols; c++) {
            int istart = _bandedmatrix_util::max(r - lbw, 0);
            int iend = _bandedmatrix_util::min(r + lbw + 1, dim);
            for (int i = istart; i < iend; i++) {
                res(r, c) += this->get(r, i) * fm.get(i, c);
            }
        }
    }
    return res;
}

}  // namespace linalg