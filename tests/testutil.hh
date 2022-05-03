#pragma once
#include <gtest/gtest.h>

#define EXPECT_MATRIX_APEQ(mat1, mat2)       \
    ASSERT_EQ(mat1.nrows, mat2.nrows);       \
    ASSERT_EQ(mat1.ncols, mat2.ncols);       \
    for (int r = 0; r < mat1.nrows; r++)     \
        for (int c = 0; c < mat1.ncols; c++) \
            ASSERT_NEAR(mat1(r, c), mat2(r, c), 1e-6);