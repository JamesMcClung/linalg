#include <gtest/gtest.h>

#include "linalg/fullmatrix.hh"
#include "linalg/permutation.hh"

using namespace linalg;

TEST(PermutationMatrix, Constructor) {
    auto id = FullMatrix({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    PermutationMatrix<3, int> p_id;

    EXPECT_EQ(id, p_id);
}

TEST(PermutationMatrix, SwapRows) {
    auto f = FullMatrix({{0, 1, 0}, {0, 0, 1}, {1, 0, 0}});
    PermutationMatrix<3, int> p;
    p.swapRows(0, 2);
    p.swapRows(1, 0);

    EXPECT_EQ(f, p);
}

TEST(PermutationMatrix, SwapCols) {
    auto f = FullMatrix({{0, 1, 0}, {0, 0, 1}, {1, 0, 0}});
    PermutationMatrix<3, int> p;
    p.swapCols(0, 2);
    p.swapCols(1, 2);

    EXPECT_EQ(f, p);
}

TEST(PermutationMatrix, SwapRowsAndCols) {
    auto f = FullMatrix({{0, 1, 0}, {0, 0, 1}, {1, 0, 0}});
    PermutationMatrix<3, int> p;
    p.swapCols(0, 2);
    p.swapRows(0, 1);

    EXPECT_EQ(f, p);
}

TEST(PermutedMatrix, PermuteRows) {
    auto f = FullMatrix({{0, 1, 2}, {3, 4, 5}, {6, 7, 8}});
    PermutationMatrix<3, int> p;
    p.swapRows(0, 2);

    auto pf = FullMatrix({{6, 7, 8}, {3, 4, 5}, {0, 1, 2}});
    EXPECT_EQ(p * f, pf);
}

TEST(PermutedMatrix, PermuteCols) {
    auto f = FullMatrix({{0, 1, 2}, {3, 4, 5}, {6, 7, 8}});
    PermutationMatrix<3, int> p;
    p.swapCols(0, 2);

    auto fp = FullMatrix({{2, 1, 0}, {5, 4, 3}, {8, 7, 6}});
    EXPECT_EQ(f * p, fp);
}

TEST(PermutedMatrix, MutateOriginalMatrix) {
    auto f = FullMatrix({{0, 1, 2}, {3, 4, 5}, {6, 7, 8}});
    auto f_ref = FullMatrix({{0, 1, 2}, {3, 4, 5}, {-1, 7, 8}});
    PermutationMatrix<3, int> p;
    p.swapRows(0, 2);

    // TODO I consider this unexpected behavior; use different syntax
    auto pf = p * f;
    pf(0, 0) = -1;  // like, this shouldn't change f, right?
    EXPECT_EQ(f, f_ref);
}