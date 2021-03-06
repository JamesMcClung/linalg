#pragma once

#include <climits>
#include <cmath>
#include <cstdint>
#include <type_traits>

#include "complex.hh"
#include "fullmatrix.hh"

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DESCRIPTION                             //
////////////////////////////////////////////////////////////////////////

// An implementation of the Cooley-Tukey FFT algorithm and its inverse.

////////////////////////////////////////////////////////////////////////
//                              UTILITY                               //
////////////////////////////////////////////////////////////////////////

namespace _fft_util {

template <typename Uint>
Uint reverse_bits_impl(Uint i, Uint nbits, std::integral_constant<int, 16>) {
    i = (i & 0xaaaa) >> 1 | (i & 0x5555) << 1;
    i = (i & 0xcccc) >> 2 | (i & 0x3333) << 2;
    i = (i & 0xf0f0) >> 4 | (i & 0x0f0f) << 4;
    i = (i & 0xff00) >> 8 | (i & 0x00ff) << 8;
    return i >> (16 - nbits);
}

// Returns the last `nbits` bits of `i` reversed
template <typename Uint>
Uint reverse_bits(Uint i, Uint nbits) {
    return reverse_bits_impl(i, nbits, std::integral_constant<int, sizeof(Uint) * CHAR_BIT>());
}

// Returns the base-2 log, rounded down.
template <typename Uint>
Uint log(Uint n) {
    Uint log = 0;
    for (Uint i = 1; i < n; i *= 2)
        ++log;
    return log;
}

template <typename T, typename R>
using test_R = T;

template <typename T, typename R = T>
struct extract_real_t {
    using real_t = R;
};

template <typename T>
struct extract_real_t<T, test_R<T, typename T::real_t>> {
    using real_t = typename extract_real_t<typename T::real_t>::real_t;
};

template <class ColVec, typename Uint, int phase_sign>
auto fft_impl(const ColVec &u) {
    if constexpr (ColVec::nrows == 1) {
        return u;
    }

    using Real = typename _fft_util::extract_real_t<ColVec>::real_t;
    constexpr Uint dim = ColVec::nrows;
    const Uint log_dim = _fft_util::log(dim);

    FullMatrix<dim, 1, Complex<Real>> uhat;
    for (Uint i = 0; i < dim; i++) {
        uhat(i, 0) = u(_fft_util::reverse_bits(i, log_dim), 0);
    }

    for (Uint m = 1; m < dim; m *= 2) {
        Complex<Real> omega_m = Complex<Real>::phase(phase_sign * M_PI / m);
        for (Uint i = 0; i < dim; i += 2 * m) {
            Complex<Real> omega = 1;
            for (Uint j = 0; j < m; j++) {
                auto a = uhat(i + j, 0);
                auto b = omega * uhat(i + j + m, 0);
                uhat(i + j, 0) = a + b;
                uhat(i + j + m, 0) = a - b;
                omega *= omega_m;
            }
        }
    }

    return uhat;
}

}  // namespace _fft_util

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

template <class ColVec, typename Uint = uint16_t>
auto fft(const ColVec &u);

template <class ColVec, typename Uint = uint16_t>
auto ifft(const ColVec &uhat);

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

/**
 * Performs a fast Fourier transform.
 * @param u a matrix with ncols=1 containing values of a periodic
 * function on a discretized space
 * @return uhat, where uhat[i, 0] = coefficient of exp(i*x) for i <
 * nrows/2, or coefficient of exp((i - nrows) * x) for i >= nrows/2
 */
template <class ColVec, typename Uint>
auto fft(const ColVec &u) {
    return _fft_util::fft_impl<ColVec, Uint, -1>(u) / ColVec::nrows;
}

/**
 * Performs a fast inverse Fourier transform.
 * @param uhat a matrix of Fourier coefficients with ncols=1
 * @return u, the original periodic function
 */
template <class ColVec, typename Uint>
auto ifft(const ColVec &uhat) {
    return _fft_util::fft_impl<ColVec, Uint, 1>(uhat);
}

}  // namespace linalg