#include <gtest/gtest.h>

#include "linalg/fullmatrix.hh"
#include "linalg/permutation.hh"
#include "linalg/recursivematrix.hh"

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

TEST(PermutationMatrix, MultiplyLeft) {
    auto f = FullMatrix({{0, 1, 2}, {3, 4, 5}, {6, 7, 8}});
    PermutationMatrix<3, int> p;
    p.swapRows(0, 2);

    auto pf = FullMatrix({{6, 7, 8}, {3, 4, 5}, {0, 1, 2}});
    EXPECT_EQ(p * f, pf);
}

TEST(PermutationMatrix, MultiplyRight) {
    auto f = FullMatrix({{0, 1, 2}, {3, 4, 5}, {6, 7, 8}});
    PermutationMatrix<3, int> p;
    p.swapCols(0, 2);

    auto fp = FullMatrix({{2, 1, 0}, {5, 4, 3}, {8, 7, 6}});
    EXPECT_EQ(f * p, fp);
}

TEST(PermutationMatrix, Multiply_Id) {
    auto f = FullMatrix({{0, 1, 2}, {3, 4, 5}, {6, 7, 8}});
    PermutationMatrix<3, int> p;

    EXPECT_EQ(f * p, f);
    EXPECT_EQ(p * f, f);
}

TEST(PermutationMatrix, MultiplyLeft_Big) {
    constexpr int dim = 256;
    RecursiveMatrix<dim, dim / 2, 4, 2, double> m, m_row_reversed;
    RecursiveMatrix<dim, dim, 4, 4, double> p_ref;
    PermutationMatrix<dim, double> p;

    for (int r = 0; r < m.nrows; r++) {
        for (int c = 0; c < m.ncols; c++) {
            m(r, c) = r * m.nrows + c;
            m_row_reversed(dim - 1 - r, c) = r * m.nrows + c;
        }
    }

    for (int r = 0; r < dim; r++)
        p_ref(r, dim - 1 - r) = 1;

    EXPECT_EQ(p * m, m);

    for (int r = 0; r < dim / 2; r++)
        p.swapRows(r, dim - 1 - r);

    EXPECT_EQ(p, p_ref);
    ASSERT_EQ(p_ref * m, m_row_reversed);
    EXPECT_EQ(p * m, m_row_reversed);
}

TEST(PermutationMatrix, Multiply_WrongWay) {
    auto f = FullMatrix({{0, 1, 2}, {3, 4, 5}, {6, 7, 8}});
    PermutationMatrix<3, int> p1, p2;

    p1.swapRows(0, 1);
    p2.swapCols(0, 1);

    ASSERT_EQ(p1, p2);

    EXPECT_EQ(p1 * f, p2 * f);
    EXPECT_EQ(f * p1, f * p2);
}