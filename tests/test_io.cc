#include <gtest/gtest.h>

#include "config.hh"
#include "matrix_io.hh"

using namespace linalg;

template <typename Real>
Real getPrepMatrixEl(int r, int c, int nrows, int ncols) {
    if (r == 0 && c == 0)
        return Real(3.14);
    return Real(r * ncols + c);
}

template <int nrows, int ncols, typename Real>
void prepMatrix(FullMatrix<nrows, ncols, Real> &f) {
    for (int r = 0; r < nrows; r++)
        for (int c = 0; c < ncols; c++)
            f(r, c) = getPrepMatrixEl<Real>(r, c, nrows, ncols);
}

TEST(MatrixIO, WriteRead) {
    auto file = TEST_OUT_DIRECTORY "MatrixIO.WriteRead.csv";
    const int nrows = 6, ncols = 9;

    FullMatrix<nrows, ncols, double> f_saved;
    prepMatrix(f_saved);
    saveMatrix(file, f_saved);

    auto f_loaded = loadMatrix<nrows, ncols, double>(file);
    ASSERT_EQ(f_saved, f_loaded);
}

TEST(MatrixIO, WriteWriteRead) {
    auto file = TEST_OUT_DIRECTORY "MatrixIO.WriteWriteRead.csv";
    const int nrows = 6, ncols = 9;

    FullMatrix<nrows, ncols, double> f_saved1(3.14);
    saveMatrix(file, f_saved1);

    FullMatrix<nrows, ncols, double> f_saved2;
    prepMatrix(f_saved2);
    saveMatrix(file, f_saved2);

    auto f_loaded = loadMatrix<nrows, ncols, double>(file);
    ASSERT_EQ(f_saved2, f_loaded);
}

TEST(MatrixIO, WriteReadRows) {
    auto file = TEST_OUT_DIRECTORY "MatrixIO.WriteReadRows.csv";
    const int nrows = 6, ncols = 9;

    FullMatrix<nrows, ncols, double> f_saved;
    prepMatrix(f_saved);
    saveMatrix(file, f_saved);

    FullMatrix<1, ncols, double> frows_loaded[nrows];
    for (int r = 0; r < nrows; r++)
        frows_loaded[r] = loadMatrix<1, ncols, double>(file, r);
    for (int r = 0; r < nrows; r++)
        for (int c = 0; c < ncols; c++)
            ASSERT_EQ(frows_loaded[r](0, c), getPrepMatrixEl<double>(r, c, nrows, ncols));
}

TEST(MatrixIO, WriteRowsRead) {
    auto file = TEST_OUT_DIRECTORY "MatrixIO.WriteRowsRead.csv";
    const int nrows = 6, ncols = 9;

    FullMatrix<1, ncols, double> frows_saved[nrows];
    for (int r = 0; r < nrows; r++)
        for (int c = 0; c < ncols; c++)
            frows_saved[r](0, c) = getPrepMatrixEl<double>(r, c, nrows, ncols);
    for (int r = 0; r < nrows; r++)
        saveMatrix(file, frows_saved[r], r > 0);

    auto f_loaded = loadMatrix<nrows, ncols, double>(file);
    for (int r = 0; r < nrows; r++)
        for (int c = 0; c < ncols; c++)
            ASSERT_EQ(f_loaded(r, c), getPrepMatrixEl<double>(r, c, nrows, ncols));
}