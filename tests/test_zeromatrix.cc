#include <gtest/gtest.h>

#include "linalg/fullmatrix.hh"
#include "linalg/zeromatrix.hh"

using namespace linalg;

TEST(ZeroMatrix, ReallyBigZMatDoesntSegfault) {
    ZeroMatrix<1000000, 1000000, double> z;
}

TEST(ZeroMatrix, Indexing) {
    ZeroMatrix<3, 6, int> z;

    for (int i = 0; i < z.nrows; i++) {
        for (int j = 0; j < z.ncols; j++) {
            EXPECT_EQ(z(i, j), 0);
        }
    }
}

TEST(ZeroMatrix, Equals) {
    ZeroMatrix<3, 4, double> z1, z2;
    FullMatrix<3, 4, double> f1(0), f2(1);

    EXPECT_EQ(z1, z1);
    EXPECT_EQ(z1, z2);
    EXPECT_EQ(z1, f1);
    EXPECT_FALSE(z1 == f2);
}

TEST(ZeroMatrix, NotEquals) {
    ZeroMatrix<3, 4, double> z1, z2;
    FullMatrix<3, 4, double> f(1);

    EXPECT_NE(z1, f);
    EXPECT_FALSE(z1 != z1);
    EXPECT_FALSE(z1 != z2);
}

TEST(ZeroMatrix, MultiplyAssignReal) {
    ZeroMatrix<3, 4, double> z1, z2;

    EXPECT_EQ(&(z1 *= 2.0), &z1);
    (z1 *= 3.5) *= 0;
    EXPECT_EQ(z1, z2);
}

TEST(ZeroMatrix, DivideAssignReal) {
    ZeroMatrix<3, 4, double> z1, z2;

    EXPECT_EQ(&(z1 /= 2.0), &z1);
    (z1 /= 3.5) /= 1;
    EXPECT_EQ(z1, z2);
}

TEST(ZeroMatrix, MultiplyReal) {
    ZeroMatrix<3, 2, double> z1;
    auto z2 = z1 * 2.0;
    auto z3 = 3.0 * z1;

    EXPECT_NE(&z1, &z2);
    EXPECT_EQ(z2, z1);

    EXPECT_NE(&z1, &z3);
    EXPECT_EQ(z3, z1);
}

TEST(ZeroMatrix, DivideReal) {
    ZeroMatrix<3, 2, double> z1;
    auto z2 = z1 / 2.;

    EXPECT_NE(&z1, &z2);
    EXPECT_EQ(z2, z1);
}

TEST(ZeroMatrix, MultiplyFZ) {
    FullMatrix<3, 2, double> f({{0, 1}, {2, 3}, {4, 5}});
    ZeroMatrix<2, 1, double> z;
    ZeroMatrix<3, 1, double> ref;

    EXPECT_EQ(f * z, ref);
}

TEST(ZeroMatrix, MultiplyZF) {
    ZeroMatrix<3, 2, double> z;
    FullMatrix<2, 1, double> f({{3}, {1}});
    ZeroMatrix<3, 1, double> ref;

    EXPECT_EQ(z * f, ref);
}

TEST(ZeroMatrix, MultiplyZZ) {
    ZeroMatrix<3, 2, double> z1;
    ZeroMatrix<2, 1, double> z2;
    ZeroMatrix<3, 1, double> ref;

    EXPECT_EQ(z1 * z2, ref);
}

TEST(ZeroMatrix, Add) {
    ZeroMatrix<3, 2, double> z1, z2;
    FullMatrix<3, 2, double> f(3.14);

    EXPECT_EQ(f, f + z1);
    EXPECT_EQ(f, z1 + f);
    EXPECT_EQ(z1, z1 + z2);
}

TEST(ZeroMatrix, Subtract) {
    ZeroMatrix<3, 2, double> z1, z2;
    FullMatrix<3, 2, double> f(3.14);

    EXPECT_EQ(f, f - z1);
    EXPECT_EQ(-f, z1 - f);
    EXPECT_EQ(z1, z1 - z2);
}

TEST(ZeroMatrix, Negate) {
    ZeroMatrix<3, 2, double> z;

    EXPECT_EQ(z, -z);
}