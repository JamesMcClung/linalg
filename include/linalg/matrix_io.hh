#pragma once

#include <fstream>
#include <string>

#include "fullmatrix.hh"

#define MATRIX_IO_DELIM " "

namespace linalg {

////////////////////////////////////////////////////////////////////////
//                            DECLARATIONS                            //
////////////////////////////////////////////////////////////////////////

template <int nrows, int ncols, typename Real>
void saveMatrix(const std::string &file_path, const Matrix<nrows, ncols, Real> &m, bool append = false);

template <int nrows, int ncols, typename Real>
FullMatrix<nrows, ncols, Real> loadMatrix(const std::string &file_path, unsigned int startLine = 0);

template <int nrows, int ncols, typename Real>
void saveMatrix_transpose(const std::string &file_path, const Matrix<nrows, ncols, Real> &m, bool append = false);

template <int nrows, int ncols, typename Real>
FullMatrix<nrows, ncols, Real> loadMatrix_transpose(const std::string &file_path, unsigned int startLine = 0);

////////////////////////////////////////////////////////////////////////
//                            DEFINITIONS                             //
////////////////////////////////////////////////////////////////////////

template <int nrows, int ncols, typename Real>
void saveMatrix(const std::string &file_path, const Matrix<nrows, ncols, Real> &m, bool append) {
    std::ofstream file(file_path, append ? std::fstream::app : 0);

    for (int r = 0; r < nrows; r++) {
        file << m(r, 0);
        for (int c = 1; c < ncols; c++) {
            file << MATRIX_IO_DELIM << m(r, c);
        }
        file << '\n';
    }

    file.close();
}

template <int nrows, int ncols, typename Real>
FullMatrix<nrows, ncols, Real> loadMatrix(const std::string &file_path, unsigned int startLine) {
    FullMatrix<nrows, ncols, Real> m;
    std::ifstream file(file_path);

    for (int line = 0; line < startLine; line++) {
        file.ignore(10000, '\n');
    }

    for (int r = 0; r < nrows; r++) {
        for (int c = 0; c < ncols; c++) {
            file >> m(r, c);
        }
    }

    file.close();
    return m;
}

template <int nrows, int ncols, typename Real>
void saveMatrix_transpose(const std::string &file_path, const Matrix<nrows, ncols, Real> &m, bool append) {
    std::ofstream file(file_path, append ? std::fstream::app : 0);

    for (int c = 0; c < ncols; c++) {
        file << m(0, c);
        for (int r = 1; r < nrows; r++) {
            file << MATRIX_IO_DELIM << m(r, c);
        }
        file << '\n';
    }

    file.close();
}

template <int nrows, int ncols, typename Real>
FullMatrix<nrows, ncols, Real> loadMatrix_transpose(const std::string &file_path, unsigned int startLine) {
    FullMatrix<nrows, ncols, Real> m;
    std::ifstream file(file_path);

    for (int line = 0; line < startLine; line++) {
        file.ignore(10000, '\n');
    }

    for (int c = 0; c < ncols; c++) {
        for (int r = 0; r < nrows; r++) {
            file >> m(r, c);
        }
    }

    file.close();
    return m;
}

}  // namespace linalg