#include <gtest/gtest.h>

#include <chrono>
#include <iostream>

#include "linalg/recursivematrix.hh"

using namespace linalg;

#define TIME std::chrono::high_resolution_clock::now()

template <class Matrix_T, class Vector_T>
auto getMultTime(int batch, int n0s_vec) {
    Matrix_T mat(0);
    Vector_T vec(0), ref(0);

    for (int i = n0s_vec; i < vec.nrows; i++)
        vec(i, 0) = i;
    for (int i = 2 * n0s_vec; i < mat.nrows; i++) {
        mat(i, i) = i;
        if (i > 0)
            mat(i, i - 1) = batch;
        ref(i, 0) = batch * (i - 1) + (i * i);
    }

    auto start = TIME;
    auto res = mat * vec;
    auto finish = TIME;

    EXPECT_EQ(res, ref);
    return std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count();
}

TEST(Efficiency, SparseMultiplication_RMatVsFMat) {
    auto f_time = getMultTime<FullMatrix<MAX_FMAT_DIM, MAX_FMAT_DIM, int>, FullMatrix<MAX_FMAT_DIM, 1, int>>(0, 16);
    auto r_time = getMultTime<RecursiveMatrix<MAX_FMAT_DIM, MAX_FMAT_DIM, 16, 16, int>, RecursiveMatrix<MAX_FMAT_DIM, 0, 16, 1, int>>(1, 16);

    const int nbatches = 8;
    for (int batch = 1; batch < nbatches; batch++) {
        f_time += getMultTime<FullMatrix<MAX_FMAT_DIM, MAX_FMAT_DIM, int>, FullMatrix<MAX_FMAT_DIM, 1, int>>(batch, 16);
        r_time += getMultTime<RecursiveMatrix<MAX_FMAT_DIM, MAX_FMAT_DIM, 16, 16, int>, RecursiveMatrix<MAX_FMAT_DIM, 1, 16, 1, int>>(batch, 16);
    }

    ASSERT_GE(f_time, r_time);

    // most recent benchmark on 2.6 GHz 6-Core Intel Core i7: ~5.8x speedup
    // std::cout << "fmat mean time: " << f_time / nbatches << " ns\nrmat mean time: " << r_time / nbatches << " ns\nspeedup: " << (float)f_time / r_time << "x\n";
}