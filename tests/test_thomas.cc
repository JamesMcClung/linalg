#include <gtest/gtest.h>

#include "bandedmatrix.hh"
#include "testutil.hh"
#include "thomas.hh"

using namespace linalg;

TEST(Thomas, Solve_Identity) {
    BandedMatrix<8, 0, 0, double> a(1);
    FullMatrix<8, 1, double> b(2);

    auto prod = a * thomas(a, b);
    EXPECT_MATRIX_APEQ(prod, b);
}

TEST(Thomas, Solve_Subdiag) {
    BandedMatrix<8, 1, 0, double> a(1);
    for (int i = 1; i < a.nrows; i++) {
        a(i, i - 1) = -1;
    }
    FullMatrix<8, 1, double> b(2);

    auto prod = a * thomas(a, b);
    EXPECT_MATRIX_APEQ(prod, b);
}

TEST(Thomas, Solve_Tridiag) {
    BandedMatrix<8, 1, 1, double> a(1);
    for (int i = 1; i < a.nrows; i++) {
        a(i, i - 1) = -1;
        a(i - 1, i) = 2;
    }
    FullMatrix<8, 1, double> b(2);

    auto prod = a * thomas(a, b);
    EXPECT_MATRIX_APEQ(prod, b);
}

TEST(Thomas, Solve_Multicol) {
    BandedMatrix<8, 1, 1, double> a(1);
    for (int i = 1; i < a.nrows; i++) {
        a(i, i - 1) = -1;
        a(i - 1, i) = 2;
    }
    FullMatrix<8, 3, double> b;
    for (int i = 0; i < b.nrows; i++) {
        for (int j = 0; j < b.ncols; j++) {
            b(i, j) = j + 1;
        }
    }

    auto prod = a * thomas(a, b);
    EXPECT_MATRIX_APEQ(prod, b);
}