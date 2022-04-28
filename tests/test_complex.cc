#include <gtest/gtest.h>

#include <type_traits>

#include "linalg/complex.hh"

using namespace linalg;
using C = Complex<double>;

TEST(Complex, ImplicitConversions) {
    C c(1.0);

    static_assert(std::is_same_v<decltype(c * 2), C>, "Expect implicit conversion to work.");
    static_assert(std::is_same_v<decltype(c / 2), C>, "Expect implicit conversion to work.");
    static_assert(std::is_same_v<decltype(c + 1), C>, "Expect implicit conversion to work.");
    static_assert(std::is_same_v<decltype(c - 1), C>, "Expect implicit conversion to work.");

    static_assert(std::is_same_v<decltype(2 * c), C>, "Expect implicit conversion to work.");
    static_assert(std::is_same_v<decltype(2 / c), C>, "Expect implicit conversion to work.");
    static_assert(std::is_same_v<decltype(1 + c), C>, "Expect implicit conversion to work.");
    static_assert(std::is_same_v<decltype(1 - c), C>, "Expect implicit conversion to work.");
}

TEST(Complex, Equals) {
    C c1(1, 2), c2(-1, 3);

    EXPECT_EQ(c1, c1);
    EXPECT_FALSE(c1 == c2);
}

TEST(Complex, NotEquals) {
    C c1(1, 2), c2(-1, 3);

    EXPECT_NE(c1, c2);
    EXPECT_FALSE(c1 != c1);
}

TEST(Complex, Constructors1) {
    C c1, c2(0), c3(0, 0), c4(c1);

    EXPECT_EQ(c1, c2);
    EXPECT_EQ(c1, c3);
    EXPECT_EQ(c1, c4);
}

TEST(Complex, Constructors2) {
    C c1(1), c2(1, 0), c3(c1);

    EXPECT_EQ(c1, c2);
    EXPECT_EQ(c1, c3);

    EXPECT_NE(c1, C(0, 1));
}

TEST(Complex, Addition) {
    C c1(1, 2), c2(-1, 3);

    EXPECT_EQ(c1 + c2, C(0, 5));
}

TEST(Complex, Subtraction) {
    C c1(1, 2), c2(-1, 3);

    EXPECT_EQ(c1 - c2, C(2, -1));
}

TEST(Complex, Negation) {
    C c(1, 2);

    EXPECT_EQ(-c, C(-1, -2));
}

TEST(Complex, Multiplication) {
    C c1(1, 2), c2(-1, 3);

    EXPECT_EQ(c1 * c2, C(-7, 1));
    EXPECT_EQ(c1 * 2, C(2, 4));
}

TEST(Complex, Division) {
    C c1(-7, 1), c2(-1, 3);

    EXPECT_EQ(c1 / c2, C(1, 2));
    EXPECT_EQ(c1 / 1, c1);
}

TEST(Complex, Conjugation) {
    C c(-7, 2);

    EXPECT_EQ(c.conj(), C(-7, -2));
}

TEST(Complex, MagnitudeSquare) {
    C c(-7, 2);

    EXPECT_EQ(c.mag2(), 53);
}

TEST(Complex, Phase) {
    C c1 = C::phase(0);
    C c2 = C::phase(M_PI / 2);

    EXPECT_EQ(c1, 1);
    EXPECT_NEAR((c2 - C(0, 1)).mag2(), 0, 1e-16);
}

TEST(Complex, AddAssign) {
    C c1(1, 2), c2(-1, 3), ref(0, 5);

    EXPECT_NE(c1, ref);
    EXPECT_EQ(&(c1 += c2), &c1);
    EXPECT_EQ(c1, ref);
}

TEST(Complex, SubtractAssign) {
    C c1(1, 2), c2(-1, 3), ref(2, -1);

    EXPECT_NE(c1, ref);
    EXPECT_EQ(&(c1 -= c2), &c1);
    EXPECT_EQ(c1, ref);
}

TEST(Complex, MultiplyAssign) {
    C c1(1, 2), c2(-1, 3), ref(-7, 1);

    EXPECT_NE(c1, ref);
    EXPECT_EQ(&(c1 *= c2), &c1);
    EXPECT_EQ(c1, ref);
}

TEST(Complex, DivideAssign) {
    C c1(-7, 1), c2(1, 2), ref(-1, 3);

    EXPECT_NE(c1, ref);
    EXPECT_EQ(&(c1 /= c2), &c1);
    EXPECT_EQ(c1, ref);
}