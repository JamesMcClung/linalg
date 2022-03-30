#include <gtest/gtest.h>

#include "fullmatrix.hh"

using namespace linalg;

TEST(FullMatrix, ConstructorArray) {
    double vals[][2] = {{2, 3.14}, {1, 4}, {0, 5}};
    FullMatrix f1(vals), f2(vals);

    EXPECT_EQ(f1, f2);
}

TEST(FullMatrix, ConstructorDefault) {
    FullMatrix<10, 12, double> f;

    bool els_are_random = false;
    for (int i = 0; i < f.nrows; i++) {
        for (int j = 0; j < f.ncols; j++) {
            els_are_random = els_are_random || f(0, 0) != f(i, j);
        }
    }
    EXPECT_TRUE(els_are_random);
}

TEST(FullMatrix, ConstructorReal) {
    FullMatrix<4, 2, double> f(3.14);

    for (int i = 0; i < f.nrows; i++) {
        for (int j = 0; j < f.ncols; j++) {
            EXPECT_EQ(f(i, j), 3.14);
        }
    }
}

TEST(FullMatrix, ConstructorCopy) {
    FullMatrix<3, 2, double> f(3.14);
    FullMatrix copy_lval(f);
    FullMatrix copy_rval(f * 1.);

    EXPECT_EQ(f, copy_lval);
    EXPECT_EQ(f, copy_rval);
    f(0, 0) = 0;
    EXPECT_NE(f, copy_lval);
    EXPECT_NE(f, copy_rval);
}

TEST(FullMatrix, Indexing) {
    FullMatrix<3, 2, double> f({{0, 1}, {2, 3}, {4, 5}});

    for (int i = 0; i < f.nrows; i++) {
        for (int j = 0; j < f.ncols; j++) {
            EXPECT_EQ(f(i, j), i * f.ncols + j);
        }
    }

    f(1, 1) = 3.14;
    EXPECT_EQ(f(1, 1), 3.14);
    f(0, 0) += 1;
    EXPECT_EQ(f(0, 0), 1);
}

TEST(FullMatrix, EqualsFMat) {
    double vals_a[][2] = {{0, 1}, {2, 3}, {4, 5}};
    double vals_b[][2] = {{0, 1}, {2, 3.14}, {4, 5}};
    FullMatrix f_a1(vals_a), f_a2(vals_a), f_b(vals_b);

    EXPECT_EQ(f_a1, f_a1);
    EXPECT_EQ(f_a1, f_a2);
    EXPECT_FALSE(f_a1 == f_b);
}

TEST(FullMatrix, NotEqualsFMat) {
    double vals_a[][2] = {{0, 1}, {2, 3}, {4, 5}};
    double vals_b[][2] = {{0, 1}, {2, 3.14}, {4, 5}};
    FullMatrix f_a1(vals_a), f_a2(vals_a), f_b(vals_b);

    EXPECT_FALSE(f_a1 != f_a1);
    EXPECT_FALSE(f_a1 != f_a2);
    EXPECT_NE(f_a1, f_b);
}

TEST(FullMatrix, MultiplyAssignReal) {
    FullMatrix f1({{0, 1}, {2, 3}, {4, 5}});
    FullMatrix f2({{0, 2}, {4, 6}, {8, 10}});

    EXPECT_NE(f1, f2);
    EXPECT_EQ(&(f1 *= 2), &f1);
    EXPECT_EQ(f1, f2);
}

TEST(FullMatrix, DivideAssignReal) {
    FullMatrix f1({{0, 1}, {2, 3}, {4, 5}});
    FullMatrix f2({{0, 2}, {4, 6}, {8, 10}});

    EXPECT_NE(f1, f2);
    EXPECT_EQ(&(f2 /= 2), &f2);
    EXPECT_EQ(f1, f2);
}

TEST(FullMatrix, MultiplyReal) {
    FullMatrix f1({{0, 1}, {2, 3}, {4, 5}});
    FullMatrix f2({{0, 2}, {4, 6}, {8, 10}});

    EXPECT_NE(f1, f2);
    EXPECT_EQ(f1 * 2, f2);
    EXPECT_EQ(2 * f1, f2);
    EXPECT_NE(f1, f2);
}

TEST(FullMatrix, DivideReal) {
    FullMatrix f1({{0, 1}, {2, 3}, {4, 5}});
    FullMatrix f2({{0, 2}, {4, 6}, {8, 10}});

    EXPECT_NE(f1, f2);
    EXPECT_EQ(f2 / 2, f1);
    EXPECT_NE(f1, f2);
}

TEST(FullMatrix, MultiplyFMat) {
    FullMatrix<3, 2, double> f1({{0, 1}, {2, 3}, {4, 5}});
    FullMatrix<2, 1, double> f2({{-1}, {10}});
    FullMatrix<3, 1, double> ref({{10}, {28}, {46}});

    EXPECT_EQ(f1 * f2, ref);
}

TEST(FullMatrix, AddAssignFMat) {
    FullMatrix f1({{0, 1}, {2, 3}, {4, 5}});
    FullMatrix f2({{0, -1}, {2, -3}, {4, 0}});
    FullMatrix ref({{0, 0}, {4, 0}, {8, 5}});

    EXPECT_NE(f1, ref);
    EXPECT_EQ(&(f1 += f2), &f1);
    EXPECT_EQ(f1, ref);
}

TEST(FullMatrix, SubtractAssignFMat) {
    FullMatrix f1({{0, 1}, {2, 3}, {4, 5}});
    FullMatrix f2({{0, -1}, {2, -3}, {4, 0}});
    FullMatrix ref({{0, 2}, {0, 6}, {0, 5}});

    EXPECT_NE(f1, ref);
    EXPECT_EQ(&(f1 -= f2), &f1);
    EXPECT_EQ(f1, ref);
}

TEST(FullMatrix, AddFMat) {
    FullMatrix f1({{0, 1}, {2, 3}, {4, 5}});
    FullMatrix f2({{0, -1}, {2, -3}, {4, 0}});
    FullMatrix ref({{0, 0}, {4, 0}, {8, 5}});

    EXPECT_NE(f1, ref);
    EXPECT_EQ(f1 + f2, ref);
}

TEST(FullMatrix, SubtractFMat) {
    FullMatrix f1({{0, 1}, {2, 3}, {4, 5}});
    FullMatrix f2({{0, -1}, {2, -3}, {4, 0}});
    FullMatrix ref({{0, 2}, {0, 6}, {0, 5}});

    EXPECT_NE(f1, ref);
    EXPECT_EQ(f1 - f2, ref);
}

TEST(FullMatrix, Negate) {
    FullMatrix f({{0, 1}, {2, 3}, {4, -5}});
    FullMatrix ref({{0, -1}, {-2, -3}, {-4, 5}});

    EXPECT_EQ(-f, ref);
}