#include <gtest/gtest.h>

#include "linalg/bandedmatrix.hh"

using namespace linalg;

TEST(BandedMatrix, ConstructorDefault_Compiles) {
    BandedMatrix<20, 3, 4, double> b;
}

template <int dim, int lbw, int ubw, typename Real>
void testBandedMatrixElsEq(const BandedMatrix<dim, lbw, ubw, Real> &bm, Real val) {
    for (int r = 0; r < dim; r++) {
        for (int c = 0; c < dim; c++) {
            if (-lbw <= c - r && c - r <= ubw)
                EXPECT_EQ(bm.get(r, c), val);
            else
                EXPECT_EQ(bm.get(r, c), Real(0));
        }
    }
}

TEST(BandedMatrix, ConstructorReal) {
    BandedMatrix<4, 1, 0, int> b1(1);
    testBandedMatrixElsEq(b1, 1);
    BandedMatrix<5, 3, 3, int> b2(2);
    testBandedMatrixElsEq(b2, 2);
    BandedMatrix<4, 0, 2, int> b3(3);
    testBandedMatrixElsEq(b3, 3);
    BandedMatrix<6, 0, 0, int> b4(4);
    testBandedMatrixElsEq(b4, 4);
}

TEST(BandedMatrix, ConstructorCopy) {
    BandedMatrix<5, 2, 1, double> b1(3.14);
    BandedMatrix<5, 2, 1, double> b2(b1);
    EXPECT_EQ(b1, b2);
    BandedMatrix<5, 0, 0, double> b3(3.15);
    BandedMatrix<5, 0, 0, double> b4(b3);
    EXPECT_EQ(b3, b4);
    BandedMatrix<5, 3, 3, double> b5(3.16);
    BandedMatrix<5, 3, 3, double> b6(b5);
    EXPECT_EQ(b5, b6);
}

TEST(BandedMatrix, ConstructorMove) {
    BandedMatrix<5, 2, 1, double> b1(3.14), b1ref(3.14);
    BandedMatrix<5, 2, 1, double> b2(std::move(b1));
    EXPECT_NE(b1, b2);
    EXPECT_EQ(b1ref, b2);
    BandedMatrix<5, 0, 0, double> b3(3.15), b3ref(3.15);
    BandedMatrix<5, 0, 0, double> b4(std::move(b3));
    EXPECT_NE(b3, b4);
    EXPECT_EQ(b3ref, b4);
    BandedMatrix<5, 3, 3, double> b5(3.16), b5ref(3.16);
    BandedMatrix<5, 3, 3, double> b6(std::move(b5));
    EXPECT_NE(b5, b6);
    EXPECT_EQ(b5ref, b6);
}

TEST(BandedMatrix, MultiplyAssign) {
    BandedMatrix<5, 2, 1, int> b1(1), b1_ref(2);
    EXPECT_NE(b1, b1_ref);
    EXPECT_EQ(&(b1 *= 2), &b1);
    EXPECT_EQ(b1, b1_ref);
    BandedMatrix<5, 3, 3, int> b2(4), b2_ref(-8);
    EXPECT_NE(b2, b2_ref);
    EXPECT_EQ(&(b2 *= -2), &b2);
    EXPECT_EQ(b2, b2_ref);
}

TEST(BandedMatrix, DivideAssign) {
    BandedMatrix<5, 2, 1, int> b1(1), b1_ref(0);
    EXPECT_NE(b1, b1_ref);
    EXPECT_EQ(&(b1 /= 2), &b1);
    EXPECT_EQ(b1, b1_ref);
    BandedMatrix<5, 3, 3, int> b2(4), b2_ref(-2);
    EXPECT_NE(b2, b2_ref);
    EXPECT_EQ(&(b2 /= -2), &b2);
    EXPECT_EQ(b2, b2_ref);
}

TEST(BandedMatrix, AddAssign) {
    BandedMatrix<5, 2, 1, int> b1(1), b1_add(2), b1_ref(3);
    EXPECT_NE(b1, b1_ref);
    EXPECT_EQ(&(b1 += b1_add), &b1);
    EXPECT_EQ(b1, b1_ref);
    BandedMatrix<5, 3, 3, int> b2(4), b2_add(-3), b2_ref(1);
    EXPECT_NE(b2, b2_ref);
    EXPECT_EQ(&(b2 += b2_add), &b2);
    EXPECT_EQ(b2, b2_ref);
}

TEST(BandedMatrix, SubtractAssign) {
    BandedMatrix<5, 2, 1, int> b1(1), b1_sub(2), b1_ref(-1);
    EXPECT_NE(b1, b1_ref);
    EXPECT_EQ(&(b1 -= b1_sub), &b1);
    EXPECT_EQ(b1, b1_ref);
    BandedMatrix<5, 3, 3, int> b2(4), b2_sub(-3), b2_ref(7);
    EXPECT_NE(b2, b2_ref);
    EXPECT_EQ(&(b2 -= b2_sub), &b2);
    EXPECT_EQ(b2, b2_ref);
}

TEST(BandedMatrix, MultiplyBMatFMat) {
    BandedMatrix<4, 2, 1, double> b(1), b_ref(b);
    FullMatrix<4, 2, double> f({{0, 1}, {2, 3}, {4, 5}, {6, 7}}), f_ref(f);
    FullMatrix<4, 2, double> bf_ref({{2, 4}, {6, 9}, {12, 16}, {12, 15}});

    EXPECT_EQ(b * f, bf_ref);
    EXPECT_EQ(b, b_ref);
    EXPECT_EQ(f, f_ref);
}