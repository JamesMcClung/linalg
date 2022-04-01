#include <gtest/gtest.h>

#include "linalg/zeromatrix.hh"

using namespace linalg;

TEST(MatrixPolymorphism, Get) {
    FullMatrix<4, 4, double> fmat(3.14);
    Matrix<4, 4, double> *mat_ptr;

    mat_ptr = new ZeroMatrix<4, 4, double>;
    ASSERT_EQ(mat_ptr->get(0, 0), 0.0);
    mat_ptr = &fmat;
    ASSERT_EQ(mat_ptr->get(0, 0), 3.14);
}

TEST(MatrixPolymorphism, IndexingConst) {
    FullMatrix<4, 4, double> fmat(3.14);
    Matrix<4, 4, double> *mat_ptr;

    mat_ptr = new ZeroMatrix<4, 4, double>;
    ASSERT_EQ(mat_ptr->operator()(0, 0), 0.0);
    mat_ptr = &fmat;
    ASSERT_EQ(mat_ptr->operator()(0, 0), 3.14);
}