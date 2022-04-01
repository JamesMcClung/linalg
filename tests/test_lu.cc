#include <gtest/gtest.h>

#include <stdexcept>

#include "linalg/fullmatrix.hh"
#include "linalg/lu_decomp.hh"
#include "linalg/recursivematrix.hh"
#include "testutil.hh"

using namespace linalg;

TEST(LU, ConstructorNoError) {
    FullMatrix<2, 2, double> a({{2, 1}, {0, 1}});
    LU_Decomp lu1(a);

    RecursiveMatrix<2, 2, 1, 1, double> b(a);
    LU_Decomp lu2(b);
    LU_Decomp lu3(std::move(b));
}

TEST(LU, LUCorrect) {
    FullMatrix<3, 3, double> a({{1, 1, 1}, {2, 3, 5}, {4, 6, 8}});
    LU_Decomp lu(a);

    FullMatrix<3, 3, double> l_ref({{1, 0, 0}, {2, 1, 0}, {4, 2, 1}});
    FullMatrix<3, 3, double> u_ref({{1, 1, 1}, {0, 1, 3}, {0, 0, -2}});

    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            ASSERT_EQ(lu.u(r, c), u_ref(r, c));
            ASSERT_EQ(lu.l(r, c), l_ref(r, c));
        }
    }
}

TEST(LU, Solve3x3) {
    FullMatrix<3, 3, double> a({{1, 1, 1}, {2, 3, 5}, {4, 6, 8}});
    FullMatrix<3, 1, double> b({{0}, {-3}, {-2}});

    LU_Decomp lu(a);
    ASSERT_EQ(a * lu.solve(b), b);
}

TEST(LU, Solve2x2) {
    FullMatrix<2, 2, double> a({{2, 1}, {0, 1}});
    FullMatrix<2, 1, double> b({{4}, {-1}});

    LU_Decomp lu(a);
    ASSERT_EQ(a * lu.solve(b), b);
}

TEST(LU, Invert) {
    RecursiveMatrix<32, 32, 4, 4, double> mat, id, inv;
    for (int r = 0; r < mat.nrows; r++) {
        int c = mat.ncols - 1 - r;
        mat(r, c) = r;
        mat(r, r) = 1;
        id(r, r) = 1;
    }
    for (int r = 0; r < mat.nrows; r++) {
        int c = mat.ncols - 1 - r;
        double a = mat(c, r), b = mat(r, c);
        inv(r, r) = -1 / (a * b - 1);
        inv(r, c) = a / (a * b - 1);
    }

    auto res = mat * LU_Decomp(mat).solve(id);
    EXPECT_MATRIX_APEQ(res, id);
}

TEST(LU, DontConstruct) {
    auto throwError = []() {
        FullMatrix a({{0, 1}, {1, 1}});
        LU_Decomp lu(a);
    };
    EXPECT_THROW({ throwError(); }, std::domain_error);
}