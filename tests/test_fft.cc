#include <gtest/gtest.h>

#include "linalg/fft.hh"

using namespace linalg;

template <typename T>
struct is_not_complex : std::true_type {};

template <typename T>
struct is_not_complex<Complex<T>> : std::false_type {};

template <typename Real>
Real mag2(Complex<Real> c) {
    return c.mag2();
}

template <typename Real, typename = std::enable_if_t<is_not_complex<Real>::value>>
Real mag2(Real r) {
    return r * r;
}

template <class M1, class M2>
auto max_diff2(const M1 &m1, const M2 &m2) {
    static_assert(M1::nrows == M2::nrows, "Cannot compare matrices of different dimensions.");
    static_assert(M1::ncols == M2::ncols, "Cannot compare matrices of different dimensions.");

    using Real = decltype(mag2(m1(0, 0) - m2(0, 0)));
    Real max_diff2 = 0;

    for (int i = 0; i < M1::nrows; i++)
        for (int j = 0; j < M1::ncols; j++)
            max_diff2 = std::max(max_diff2, mag2(m1(i, j) - m2(i, j)));

    return max_diff2;
}

TEST(FFT, Inverse) {
    double len_x = 2 * M_PI;
    FullMatrix<32, 1, double> u;
    for (int i = 0; i < u.nrows; i++)
        u(i, 0) = std::sin(i * len_x / u.nrows);

    EXPECT_LT(max_diff2(ifft(fft(u)), u), 1e-6);
}

TEST(FFT, FullSpectrum) {
    double len_x = 2 * M_PI;
    FullMatrix<32, 1, Complex<double>> u(0);
    FullMatrix<32, 1, Complex<double>> uhat_ref;
    for (int i = 0; i < u.nrows; i++) {
        for (int n = -u.nrows / 2; n < u.nrows / 2; n++)
            u(i, 0) += n * Complex<double>::phase(n * i * len_x / u.nrows);
        uhat_ref(i, 0) = i < u.nrows / 2 ? i : i - u.nrows;
    }

    EXPECT_LT(max_diff2(fft(u), uhat_ref), 1e-6);
}
