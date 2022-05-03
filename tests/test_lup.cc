#include <gtest/gtest.h>

#include "linalg/fullmatrix.hh"
#include "linalg/lup_decomp.hh"
#include "linalg/recursivematrix.hh"
#include "testutil.hh"

using namespace linalg;

TEST(LUP, ConstructorNoError) {
    FullMatrix<2, 2, double> a({{2, 1}, {0, 1}});
    LUP_Decomp lu1(a);

    RecursiveMatrix<2, 2, 1, 1, double> b(a);
    LUP_Decomp lu2(b);
    LUP_Decomp lu3(std::move(b));
}

TEST(LUP, LUCorrect) {
    FullMatrix<3, 3, double> a({{1, 1, 1}, {2, 3, 5}, {4, 6, 8}});
    LUP_Decomp lup(a);

    FullMatrix<3, 3, double> l_ref({{1, 0, 0}, {.25, 1, 0}, {.5, 0, 1}});
    FullMatrix<3, 3, double> u_ref({{4, 6, 8}, {0, -.5, -1}, {0, 0, 1}});
    FullMatrix<3, 3, double> p_ref({{0, 0, 1}, {1, 0, 0}, {0, 1, 0}});

    ASSERT_EQ(lup.get_L(), l_ref);
    ASSERT_EQ(lup.get_U(), u_ref);
    ASSERT_EQ(lup.get_P(), p_ref);
}

TEST(LUP, Solve3x3) {
    FullMatrix<3, 3, double> a({{1, 1, 1}, {2, 3, 5}, {4, 6, 8}});
    FullMatrix<3, 1, double> b({{0}, {-3}, {-2}});

    LUP_Decomp lup(a);
    ASSERT_EQ(a * lup.solve(b), b);
}

TEST(LUP, Solve2x2) {
    FullMatrix<2, 2, double> a({{2, 1}, {0, 1}});
    FullMatrix<2, 1, double> b({{4}, {-1}});

    LUP_Decomp lup(a);
    ASSERT_EQ(a * lup.solve(b), b);
}

TEST(LUP, Invert) {
    constexpr int dim = 32;
    RecursiveMatrix<dim, dim, 4, 4, double> mat, id;
    for (int r = 0; r < mat.nrows; r++) {
        int c = mat.ncols - 1 - r;
        mat(r, c) = r;
        mat(r, r) = 1;
        id(r, r) = 1;
    }
    for (int r = 0; r < mat.nrows; r++) {
        int c = mat.ncols - 1 - r;
        double a = mat(c, r), b = mat(r, c);
    }

    auto res = LUP_Decomp(mat).solve(id) * mat;
    EXPECT_MATRIX_APEQ(res, id);
}