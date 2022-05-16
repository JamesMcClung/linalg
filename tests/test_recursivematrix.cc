#include <gtest/gtest.h>

#include "linalg/recursivematrix.hh"

using namespace linalg;

template <typename M, typename Real>
void test_elsEq(M& a, Real val) {
    for (int i = 0; i < a.nrows; i++) {
        for (int j = 0; j < a.ncols; j++) {
            EXPECT_EQ(a.get(i, j), val);
        }
    }
}

TEST(RecursiveMatrix, ConstructorDefault_Small) {
    RecursiveMatrix<4, 2, 2, 2, double> a;
    test_elsEq(a, 0.0);
}

TEST(RecursiveMatrix, ConstructorDefault_Big) {
    RecursiveMatrix<MAX_FMAT_DIM, 1, 4, 1, double> a;
    test_elsEq(a, 0.0);
}

TEST(RecursiveMatrix, ConstructorReal_Small) {
    RecursiveMatrix<4, 2, 2, 2, double> a(3.14);
    test_elsEq(a, 3.14);
}

TEST(RecursiveMatrix, ConstructorReal_Big) {
    RecursiveMatrix<MAX_FMAT_DIM, 1, 4, 1, double> a(3.14);
    test_elsEq(a, 3.14);
}

TEST(RecursiveMatrix, ConstructorMatrix) {
    FullMatrix<1, 4, double> f({{0, 1, 2, 3}});
    RecursiveMatrix<1, 4, 1, 2, double> r(f);
    EXPECT_EQ(f, r);
}

TEST(RecursiveMatrix, ConstructorZMat) {
    ZeroMatrix<MAX_FMAT_DIM, 1, double> z;
    RecursiveMatrix<MAX_FMAT_DIM, 1, 4, 1, double> r(z);
    EXPECT_EQ(z, r);
}

TEST(RecursiveMatrix, ConstructorCopy) {
    RecursiveMatrix<MAX_FMAT_DIM, 2, 8, 2, double> r1, ref;
    r1(0, 1) = 3.14;
    ref(0, 1) = 3.14;
    RecursiveMatrix<MAX_FMAT_DIM, 2, 8, 2, double> r2(r1);

    EXPECT_EQ(r1, ref);
    EXPECT_EQ(r2, ref);
}

TEST(RecursiveMatrix, ConstructorMove) {
    RecursiveMatrix<MAX_FMAT_DIM, 2, 8, 2, double> r1, ref;
    r1(0, 1) = 3.14;
    ref(0, 1) = 3.14;
    RecursiveMatrix<MAX_FMAT_DIM, 2, 8, 2, double> r2(std::move(r1));

    EXPECT_NE(r1, ref);
    EXPECT_EQ(r2, ref);
}

TEST(RecursiveMatrix, AssignCopy) {
    RecursiveMatrix<MAX_FMAT_DIM, 2, 8, 2, double> r1, ref;
    r1(0, 1) = 3.14;
    ref(0, 1) = 3.14;
    RecursiveMatrix<MAX_FMAT_DIM, 2, 8, 2, double> r2;

    ASSERT_EQ(r1, ref);
    EXPECT_NE(r2, ref);
    r2 = r1;
    EXPECT_EQ(r1, ref);
    EXPECT_EQ(r2, ref);
}

TEST(RecursiveMatrix, AssignMove) {
    RecursiveMatrix<MAX_FMAT_DIM, 2, 8, 2, double> r1, ref;
    r1(0, 1) = 3.14;
    ref(0, 1) = 3.14;
    RecursiveMatrix<MAX_FMAT_DIM, 2, 8, 2, double> r2;

    ASSERT_EQ(r1, ref);
    EXPECT_NE(r2, ref);
    r2 = std::move(r1);
    EXPECT_NE(r1, ref);
    EXPECT_EQ(r2, ref);
}

template <typename M>
void test_indexing(M& a) {
    int r = 0;
    int c = 0;
    EXPECT_EQ(a.get(r, c), 0.0);
    a(r, c) = 3.14;
    EXPECT_EQ(a.get(r, c), 3.14);

    r = a.nrows - 1;
    c = a.ncols / 4;
    EXPECT_EQ(a(r, c), 0.0);
    a(r, c) += 1;
    EXPECT_EQ(a(r, c), 1.0);

    r--;
    EXPECT_EQ(a.get(r, c), 0.0);
    a(r, c) -= 1;
    a(r, c) *= 2;
    EXPECT_EQ(a.get(r, c), -2.0);
}

TEST(RecursiveMatrix, Indexing_Small) {
    RecursiveMatrix<4, 2, 2, 2, double> a;
    test_indexing(a);
}

TEST(RecursiveMatrix, Indexing_Big) {
    RecursiveMatrix<MAX_FMAT_DIM, 128, 4, 4, double> a;
    test_indexing(a);
}

template <typename M1, typename M2>
void test_eq(M1& m1, const M2& m2) {
    EXPECT_EQ(m1, m2);
    EXPECT_EQ(m2, m1);
    EXPECT_FALSE(m1 != m2);
    EXPECT_FALSE(m2 != m1);

    m1(0, 0) = -1;
    EXPECT_FALSE(m1 == m2);
    EXPECT_FALSE(m2 == m1);
    EXPECT_NE(m1, m2);
    EXPECT_NE(m2, m1);

    m1(0, 0) = m2(0, 0);
    EXPECT_EQ(m1, m2);
    EXPECT_EQ(m2, m1);
    EXPECT_FALSE(m1 != m2);
    EXPECT_FALSE(m2 != m1);
}

TEST(RecursiveMatrix, EqualsRMat_Small) {
    RecursiveMatrix<128, 128, 4, 4, double> a;
    const RecursiveMatrix<128, 128, 2, 4, double> b;

    EXPECT_EQ(a, a);
    test_eq(a, b);
}

TEST(RecursiveMatrix, EqualsZMat_Small) {
    RecursiveMatrix<128, 128, 4, 4, double> r;
    const ZeroMatrix<128, 128, double> z;

    test_eq(r, z);
}

TEST(RecursiveMatrix, EqualsFMat) {
    RecursiveMatrix<16, 16, 4, 4, double> r(1);
    const FullMatrix<16, 16, double> f(1);

    test_eq(r, f);
}

TEST(RecursiveMatrix, EqualsRMat_BigSlow) {
    RecursiveMatrix<8, MAX_FMAT_DIM, 1, 4, double> a;
    const RecursiveMatrix<8, MAX_FMAT_DIM, 1, 8, double> b;

    EXPECT_EQ(a, a);
    test_eq(a, b);
}

TEST(RecursiveMatrix, EqualsRMat_BigFast) {
    RecursiveMatrix<8, MAX_FMAT_DIM, 1, 8, double> a;
    const RecursiveMatrix<8, MAX_FMAT_DIM, 1, 8, double> b;

    EXPECT_EQ(a, a);
    test_eq(a, b);
}

TEST(RecursiveMatrix, EqualsZMat_Big) {
    RecursiveMatrix<8, MAX_FMAT_DIM, 2, 8, double> r;
    const ZeroMatrix<8, MAX_FMAT_DIM, double> z;

    test_eq(r, z);
}

TEST(RecursiveMatrix, AddRMat) {
    RecursiveMatrix<16, MAX_FMAT_DIM, 16, 8, double> r1, r2, ref;
    // r1(0, 0) = 0;
    // r1(1, 0) = 0;
    // r1(2, 0) = 0;
    r1(3, 0) = 0;
    r1(4, 0) = 0;
    r1(5, 0) = 0;
    r1(6, 0) = 3.5;
    r1(7, 0) = 3.5;
    r1(8, 0) = 3.5;

    // r2(0, 0) = 0;
    r2(1, 0) = 0;
    r2(2, 0) = 3.5;
    // r2(3, 0) = 0;
    r2(4, 0) = 0;
    r2(5, 0) = 3.5;
    // r2(6, 0) = 0;
    r2(7, 0) = 0;
    r2(8, 0) = 3.5;

    ref(2, 0) = 3.5;
    ref(5, 0) = 3.5;
    ref(6, 0) = 3.5;
    ref(7, 0) = 3.5;
    ref(8, 0) = 7;

    EXPECT_NE(r1, ref);
    EXPECT_EQ(r1 + r2, ref);
}

TEST(RecursiveMatrix, SubtractRMat) {
    RecursiveMatrix<16, MAX_FMAT_DIM, 16, 8, double> r1, r2, ref;
    // r1(0, 0) = 0;
    // r1(1, 0) = 0;
    // r1(2, 0) = 0;
    r1(3, 0) = 0;
    r1(4, 0) = 0;
    r1(5, 0) = 0;
    r1(6, 0) = 3.5;
    r1(7, 0) = 3.5;
    r1(8, 0) = 3.5;

    // r2(0, 0) = 0;
    r2(1, 0) = 0;
    r2(2, 0) = 3.5;
    // r2(3, 0) = 0;
    r2(4, 0) = 0;
    r2(5, 0) = 3.5;
    // r2(6, 0) = 0;
    r2(7, 0) = 0;
    r2(8, 0) = 3.5;

    ref(2, 0) = -3.5;
    ref(5, 0) = -3.5;
    ref(6, 0) = 3.5;
    ref(7, 0) = 3.5;

    EXPECT_NE(r1, ref);
    EXPECT_EQ(r1 - r2, ref);
}

TEST(RecursiveMatrix, Negate) {
    RecursiveMatrix<4, MAX_FMAT_DIM, 4, 8, double> r, ref;
    // r(0, 0) = 0;
    r(1, 0) = 0;
    r(2, 0) = 3.14;

    ref(2, 0) = -3.14;

    EXPECT_NE(r, ref);
    EXPECT_EQ(-r, ref);
}

TEST(RecursiveMatrix, AddZMat) {
    RecursiveMatrix<16, MAX_FMAT_DIM, 16, 8, double> r, ref;
    ZeroMatrix<16, MAX_FMAT_DIM, double> z;
    // r(0, 0) = 0;
    r(1, 0) = 0;
    r(2, 0) = 3.14;

    ref(2, 0) = 3.14;

    EXPECT_EQ(r, ref);
    EXPECT_EQ(r + z, ref);
    EXPECT_EQ(z + r, ref);
}

TEST(RecursiveMatrix, SubtractZMat) {
    RecursiveMatrix<16, MAX_FMAT_DIM, 16, 8, double> r, posref, negref;
    ZeroMatrix<16, MAX_FMAT_DIM, double> z;
    // r(0, 0) = 0;
    r(1, 0) = 0;
    r(2, 0) = 3.14;

    posref(2, 0) = 3.14;
    negref(2, 0) = -3.14;

    EXPECT_EQ(r, posref);
    EXPECT_EQ(r - z, posref);
    EXPECT_NE(r, negref);
    EXPECT_EQ(z - r, negref);
}

TEST(RecursiveMatrix, AddAssignRMat) {
    RecursiveMatrix<16, MAX_FMAT_DIM, 16, 8, double> r1, r2, ref;
    // r1(0, 0) = 0;
    // r1(1, 0) = 0;
    // r1(2, 0) = 0;
    r1(3, 0) = 0;
    r1(4, 0) = 0;
    r1(5, 0) = 0;
    r1(6, 0) = 3.5;
    r1(7, 0) = 3.5;
    r1(8, 0) = 3.5;

    // r2(0, 0) = 0;
    r2(1, 0) = 0;
    r2(2, 0) = 3.5;
    // r2(3, 0) = 0;
    r2(4, 0) = 0;
    r2(5, 0) = 3.5;
    // r2(6, 0) = 0;
    r2(7, 0) = 0;
    r2(8, 0) = 3.5;

    ref(2, 0) = 3.5;
    ref(5, 0) = 3.5;
    ref(6, 0) = 3.5;
    ref(7, 0) = 3.5;
    ref(8, 0) = 7;

    EXPECT_NE(r1, ref);
    EXPECT_EQ(&(r1 += r2), &r1);
    EXPECT_EQ(r1, ref);
}

TEST(RecursiveMatrix, SubtractAssignRMat) {
    RecursiveMatrix<16, MAX_FMAT_DIM, 16, 8, double> r1, r2, ref;
    // r1(0, 0) = 0;
    // r1(1, 0) = 0;
    // r1(2, 0) = 0;
    r1(3, 0) = 0;
    r1(4, 0) = 0;
    r1(5, 0) = 0;
    r1(6, 0) = 3.5;
    r1(7, 0) = 3.5;
    r1(8, 0) = 3.5;

    // r2(0, 0) = 0;
    r2(1, 0) = 0;
    r2(2, 0) = 3.5;
    // r2(3, 0) = 0;
    r2(4, 0) = 0;
    r2(5, 0) = 3.5;
    // r2(6, 0) = 0;
    r2(7, 0) = 0;
    r2(8, 0) = 3.5;

    ref(2, 0) = -3.5;
    ref(5, 0) = -3.5;
    ref(6, 0) = 3.5;
    ref(7, 0) = 3.5;

    EXPECT_NE(r1, ref);
    EXPECT_EQ(&(r1 -= r2), &r1);
    EXPECT_EQ(r1, ref);
}

TEST(RecursiveMatrix, MultiplyReal) {
    RecursiveMatrix<4, MAX_FMAT_DIM, 16, 8, double> r, ref;
    // r1(0, 0) = 0;
    r(1, 0) = 0;
    r(2, 0) = 3.5;

    ref(2, 0) = 7;

    EXPECT_NE(r, ref);
    EXPECT_EQ(r * 2, ref);
    EXPECT_EQ(2.0 * r, ref);
}

TEST(RecursiveMatrix, DivideReal) {
    RecursiveMatrix<4, MAX_FMAT_DIM, 16, 8, double> r, ref;
    // r1(0, 0) = 0;
    r(1, 0) = 0;
    r(2, 0) = 7;

    ref(2, 0) = 3.5;

    EXPECT_NE(r, ref);
    EXPECT_EQ(r / 2, ref);
}

TEST(RecursiveMatrix, MultiplyAssignReal) {
    RecursiveMatrix<4, MAX_FMAT_DIM, 16, 8, double> r, ref;
    // r1(0, 0) = 0;
    r(1, 0) = 0;
    r(2, 0) = 3.5;

    ref(2, 0) = 7;

    EXPECT_NE(r, ref);
    EXPECT_EQ(&(r *= 2), &r);
    EXPECT_EQ(r, ref);
}

TEST(RecursiveMatrix, DivideAssignReal) {
    RecursiveMatrix<4, MAX_FMAT_DIM, 16, 8, double> r, ref;
    // r1(0, 0) = 0;
    r(1, 0) = 0;
    r(2, 0) = 7;

    ref(2, 0) = 3.5;

    EXPECT_NE(r, ref);
    EXPECT_EQ(&(r /= 2), &r);
    EXPECT_EQ(r, ref);
}

TEST(RecursiveMatrix, MultiplyRMat) {
    RecursiveMatrix<MAX_FMAT_DIM, MAX_FMAT_DIM, 16, 16, double> rmat;
    RecursiveMatrix<MAX_FMAT_DIM, 1, 16, 1, double> rvec, ref;

    for (int i = 16; i < rvec.nrows; i++)
        rvec(i, 0) = i;
    for (int i = 32; i < rmat.nrows; i++) {
        rmat(i, i) = i;
        if (i > 0)
            rmat(i, i - 1) = 1;
        ref(i, 0) = (i - 1) + (i * i);
    }

    EXPECT_NE(rvec, ref);
    EXPECT_EQ(rmat * rvec, ref);
}