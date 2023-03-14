#ifndef MATRIXMATH_H
#define MATRIXMATH_H

#include "CoolPropTools.h"
#include "Exceptions.h"

#include <vector>
#include <string>
#include <numeric>  // inner_product
#include <sstream>
#include "float.h"

#include "Eigen/Core"

/// A wrapper around std::vector
/** This wrapper makes the standard vector multi-dimensional.
 *  A useful thing even though we might not need it that
 *  much. However, it makes the code look better and the
 *  polynomial class really is a mess...
 *  Source: http://stackoverflow.com/questions/13105514/n-dimensional-vector
 */
template <size_t dimcount, typename T>
struct VectorNd
{
    typedef std::vector<typename VectorNd<dimcount - 1, T>::type> type;
};
template <typename T>
struct VectorNd<0, T>
{
    typedef T type;
};

namespace CoolProp {

/// Some shortcuts and regularly needed operations
template <class T>
std::size_t num_rows(std::vector<T> const& in) {
    return in.size();
}
template <class T>
std::size_t num_rows(std::vector<std::vector<T>> const& in) {
    return in.size();
}

template <class T>
std::size_t max_cols(std::vector<std::vector<T>> const& in) {
    std::size_t cols = 0;
    std::size_t col = 0;
    for (std::size_t i = 0; i < in.size(); i++) {
        col = in[i].size();
        if (cols < col) {
            cols = col;
        }
    }
    return cols;
};
template <class T>
bool is_squared(std::vector<std::vector<T>> const& in) {
    std::size_t cols = max_cols(in);
    if (cols != num_rows(in)) {
        return false;
    } else {
        for (std::size_t i = 0; i < in.size(); i++) {
            if (cols != in[i].size()) {
                return false;
            }
        }
    }
    return true;
};

template <class T>
std::size_t num_cols(std::vector<T> const& in) {
    return 1;
}
template <class T>
std::size_t num_cols(std::vector<std::vector<T>> const& in) {
    if (num_rows(in) > 0) {
        if (is_squared(in)) {
            return in[0].size();
        } else {
            return max_cols(in);
        }
    } else {
        return 0;
    }
};

/// Convert vectors and matrices
/** Conversion functions for the different kinds of object-like
 *  parameters. This might be obsolete at a later stage, but now
 *  it is still needed.
 *  @param coefficients matrix containing the ordered coefficients
 *  @param axis axis along which to extract
 */
template <typename T>
std::vector<T> eigen_to_vec1D(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& coefficients, int axis = 0) {
    std::vector<T> result;
    size_t r = coefficients.rows(), c = coefficients.cols();
    if (axis == 0) {
        if (c != 1) throw ValueError(format("Your matrix has the wrong dimensions: %d,%d", r, c));
        result.resize(r);
        for (size_t i = 0; i < r; ++i) {
            result[i] = coefficients(i, 0);
        }
    } else if (axis == 1) {
        if (r != 1) throw ValueError(format("Your matrix has the wrong dimensions: %d,%d", r, c));
        result.resize(c);
        for (size_t i = 0; i < c; ++i) {
            result[i] = coefficients(0, i);
        }
    } else {
        throw ValueError(format("You have to provide axis information: %d is not valid. ", axis));
    }
    return result;
}
/// @param coefficients matrix containing the ordered coefficients
template <class T>
std::vector<std::vector<T>> eigen_to_vec(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& coefficients) {
    // Eigen uses columns as major axis, this might be faster than the row iteration.
    // However, the 2D vector stores things differently, no idea what is faster...
    std::vector<std::vector<T>> result;
    size_t r = coefficients.rows(), c = coefficients.cols();
    result.resize(r, std::vector<T>(c, 0));  // extends vector if necessary
    for (size_t i = 0; i < r; ++i) {
        result[i].resize(c, 0);
        for (size_t j = 0; j < c; ++j) {
            result[i][j] = coefficients(i, j);
        }
    }
    return result;
}

/// @param coefficients matrix containing the ordered coefficients
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vec_to_eigen(const std::vector<std::vector<T>>& coefficients) {
    size_t nRows = num_rows(coefficients), nCols = num_cols(coefficients);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(nRows, nCols);
    for (size_t i = 0; i < nCols; ++i) {
        for (size_t j = 0; j < nRows; ++j) {
            result(j, i) = coefficients[j][i];
        }
    }
    return result;
}
/**
 * @param coefficients matrix containing the ordered coefficients
 * @param axis axis along which to extract data
 */
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vec_to_eigen(const std::vector<T>& coefficients, int axis = 0) {
    size_t nRows = num_rows(coefficients);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result;
    if (axis == 0)
        result.resize(nRows, 1);
    else if (axis == 1)
        result.resize(1, nRows);
    else
        throw ValueError(format("You have to provide axis information: %d is not valid. ", axis));
    for (size_t i = 0; i < nRows; ++i) {
        if (axis == 0) result(i, 0) = coefficients[i];
        if (axis == 1) result(0, i) = coefficients[i];
    }
    return result;
}
/// @param coefficient
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vec_to_eigen(const T& coefficient) {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(1, 1);
    result(0, 0) = coefficient;
    return result;
}

/// Convert 1D matrix to vector
/** Returns either a row- or a column-based
 *  vector. By default, Eigen prefers column
 *  major ordering, just like Fortran.
 */

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1> makeColVector(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix) {
    std::size_t r = matrix.rows();
    std::size_t c = matrix.cols();
    Eigen::Matrix<T, Eigen::Dynamic, 1> vector;
    if (r == 1 && c >= 1) {  // Check passed, matrix can be transformed
        vector = matrix.transpose().block(0, 0, c, r);
    } else if (r >= 1 && c == 1) {  // Check passed, matrix can be transformed
        vector = matrix.block(0, 0, r, c);
    } else {  // Check failed, throw error
        throw ValueError(format("Your matrix (%d,%d) cannot be converted into a vector (x,1).", r, c));
    }
    return vector;
}
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1> makeVector(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix) {
    return makeColVector(matrix);
}
template <class T>
Eigen::Matrix<T, 1, Eigen::Dynamic> makeRowVector(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix) {
    std::size_t r = matrix.rows();
    std::size_t c = matrix.cols();
    Eigen::Matrix<T, 1, Eigen::Dynamic> vector;
    if (r == 1 && c >= 1) {  // Check passed, matrix can be transformed
        vector = matrix.block(0, 0, r, c);
    } else if (r >= 1 && c == 1) {  // Check passed, matrix can be transformed
        vector = matrix.transpose().block(0, 0, c, r);
    } else {  // Check failed, throw error
        throw ValueError(format("Your matrix (%d,%d) cannot be converted into a vector (1,x).", r, c));
    }
    return vector;
}

/// Remove rows and columns from matrices
/** A set of convenience functions inspired by http://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
 *  but altered to respect templates.
 */
template <class T>
void removeRow(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix, std::size_t rowToRemove) {
    //template<class T> void removeRow(Eigen::MatrixXd& matrix, std::size_t rowToRemove){
    //void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove){
    //template <typename Derived> void removeRow(Eigen::MatrixBase<Derived> &matrix, std::size_t rowToRemove){
    std::size_t numRows = matrix.rows() - 1;
    std::size_t numCols = matrix.cols();
    if (rowToRemove < numRows) {
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
    } else {
        if (rowToRemove > numRows) {
            throw ValueError(format("Your matrix does not have enough rows, %d is not greater or equal to %d.", numRows, rowToRemove));
        }
        // Do nothing, resize removes the last row
    }
    matrix.conservativeResize(numRows, numCols);
}

template <class T>
void removeColumn(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix, std::size_t colToRemove) {
    //template<class T> void removeColumn(Eigen::MatrixXd& matrix, std::size_t colToRemove){
    //void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove){
    //template <typename Derived> void removeColumn(Eigen::MatrixBase<Derived> &matrix, std::size_t colToRemove){
    std::size_t numRows = matrix.rows();
    std::size_t numCols = matrix.cols() - 1;
    if (colToRemove < numCols) {
        matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
    } else {
        if (colToRemove > numCols) {
            throw ValueError(format("Your matrix does not have enough columns, %d is not greater or equal to %d.", numCols, colToRemove));
        }
        // Do nothing, resize removes the last column
    }
    matrix.conservativeResize(numRows, numCols);
}

///// @param coefficients matrix containing the ordered coefficients
//template <class T> Eigen::Matrix<T, Eigen::Dynamic,Eigen::Dynamic> convert(const std::vector<std::vector<T> > &coefficients){
//    size_t nRows = num_rows(coefficients), nCols = num_cols(coefficients);
//    Eigen::MatrixBase<T> result(nRows,nCols);
//    for (size_t i = 0; i < nCols; ++i) {
//        for (size_t j = 0; j < nRows; ++j) {
//            result(j,i) = coefficients[j][i];
//        }
//    }
//    return result;
//}
//
///// @param coefficients matrix containing the ordered coefficients
//template <class T, int R, int C> void convert(const std::vector<std::vector<T> > &coefficients, Eigen::Matrix<T,R,C> &result){
//    size_t nRows = num_rows(coefficients), nCols = num_cols(coefficients);
//    //Eigen::MatrixBase<T> result(nRows,nCols);
//    for (size_t i = 0; i < nCols; ++i) {
//        for (size_t j = 0; j < nRows; ++j) {
//            result(j,i) = coefficients[j][i];
//        }
//    }
//    //return result;
//}

//
//template <class T> void convert(const std::vector<std::vector<T> > &coefficients, Eigen::MatrixBase<T> &result){
//    size_t nRows = num_rows(coefficients), nCols = num_cols(coefficients);
//    //Eigen::MatrixBase<T> result;
//    //if ((R!=nRows) || (C!=nCols))
//    result.resize(nRows,nCols);
//    for (size_t i = 0; i < nCols; ++i) {
//        for (size_t j = 0; j < nRows; ++j) {
//            result(j,i) = coefficients[j][i];
//        }
//    }
//    //return result;
//}

//template <class Derived>
//inline void func1(MatrixBase<Derived> &out_mat ){
//  // Do something then return a matrix
//  out_mat = ...
//}

//template <class Derived>
//Eigen::Matrix<class Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>
//Multiply(const Eigen::MatrixBase<DerivedA>& p1,
//    const Eigen::MatrixBase<DerivedB>& p2)
//{
//    return (p1 * p2);
//}
//
//
//template <typename DerivedA, typename DerivedB>
//Eigen::Matrix<typename DerivedA::Scalar, DerivedA::RowsAtCompileTime, DerivedB::ColsAtCompileTime>
//Multiply(const Eigen::MatrixBase<DerivedA>& p1,
//    const Eigen::MatrixBase<DerivedB>& p2)
//{
//    return (p1 * p2);
//}

/// Templates for printing numbers, vectors and matrices
static const char* stdFmt = "%8.3f";

///Templates for turning vectors (1D-matrices) into strings
template <class T>
std::string vec_to_string(const std::vector<T>& a, const char* fmt) {
    if (a.size() < 1) return std::string("");
    std::stringstream out;
    out << "[ " << format(fmt, a[0]);
    for (size_t j = 1; j < a.size(); j++) {
        out << ", " << format(fmt, a[j]);
    }
    out << " ]";
    return out.str();
};
template <class T>
std::string vec_to_string(const std::vector<T>& a) {
    return vec_to_string(std::vector<double>(a.begin(), a.end()), stdFmt);
};
///Templates for turning vectors (1D-matrices) into strings
inline std::string stringvec_to_string(const std::vector<std::string>& a) {
    if (a.size() < 1) return std::string("");
    std::stringstream out;
    out << "[ " << format("%s", a[0].c_str());
    for (size_t j = 1; j < a.size(); j++) {
        out << ", " << format("%s", a[j].c_str());
    }
    out << " ]";
    return out.str();
};

/// Templates for turning numbers (0D-matrices) into strings
template <class T>
std::string vec_to_string(const T& a, const char* fmt) {
    std::vector<T> vec;
    vec.push_back(a);
    return vec_to_string(vec, fmt);
};
template <class T>
std::string vec_to_string(const T& a) {
    return vec_to_string((double)a, stdFmt);
};

///Templates for turning 2D-matrices into strings
template <class T>
std::string vec_to_string(const std::vector<std::vector<T>>& A, const char* fmt) {
    if (A.size() < 1) return std::string("");
    std::stringstream out;
    out << "[ " << vec_to_string(A[0], fmt);
    for (size_t j = 1; j < A.size(); j++) {
        out << ", " << std::endl << "  " << vec_to_string(A[j], fmt);
    }
    out << " ]";
    return out.str();
};
template <class T>
std::string vec_to_string(const std::vector<std::vector<T>>& A) {
    return vec_to_string(A, stdFmt);
};

///Templates for turning Eigen matrices into strings
template <class T>
std::string mat_to_string(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, const char* fmt) {
    //std::string mat_to_string(const Eigen::MatrixXd &A, const char *fmt) {
    std::size_t r = A.rows();
    std::size_t c = A.cols();
    if ((r < 1) || (c < 1)) return std::string("");
    std::stringstream out;
    out << "[ ";
    if (r == 1) {
        out << format(fmt, A(0, 0));
        for (size_t j = 1; j < c; j++) {
            out << ", " << format(fmt, A(0, j));
        }
    } else {
        out << mat_to_string(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(A.row(0)), fmt);
        for (size_t i = 1; i < r; i++) {
            out << ", " << std::endl << "  " << mat_to_string(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(A.row(i)), fmt);
        }
    }
    out << " ]";
    return out.str();
};
template <class T>
std::string mat_to_string(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A) {
    //std::string vec_to_string(const Eigen::MatrixXd &A) {
    return mat_to_string(A, stdFmt);
};

///// Templates for printing numbers, vectors and matrices
//static const char* stdFmt = "%8.3f";
//
/////Templates for turning vectors (1D-matrices) into strings
//template<class T> std::string vec_to_string(const             std::vector<T>   &a, const char *fmt) {
//    if (a.size()<1) return std::string("");
//    std::stringstream out;
//    out << "[ " << format(fmt,a[0]);
//    for (size_t j = 1; j < a.size(); j++) {
//        out << ", " << format(fmt, a[j]);
//    }
//    out << " ]";
//    return out.str();
//};
//template<class T> std::string vec_to_string(const             std::vector<T>   &a) {
//    return vec_to_string(a, stdFmt);
//};
//
///// Templates for turning numbers (0D-matrices) into strings
//template<class T> std::string vec_to_string(const                         T    &a, const char *fmt) {
//    std::vector<T> vec;
//    vec.push_back(a);
//    return vec_to_string(vec, fmt);
//};
//template<class T> std::string vec_to_string(const                         T    &a) {
//    return vec_to_string(a, stdFmt);
//};
//
/////Templates for turning 2D-matrices into strings
//template<class T> std::string vec_to_string(const std::vector<std::vector<T> > &A, const char *fmt) {
//    if (A.size()<1) return std::string("");
//    std::stringstream out;
//    out << "[ " << vec_to_string(A[0], fmt);
//    for (size_t j = 1; j < A.size(); j++) {
//        out << ", " << std::endl << "  " << vec_to_string(A[j], fmt);
//    }
//    out << " ]";
//    return out.str();
//};
//template<class T> std::string vec_to_string(const std::vector<std::vector<T> > &A) {
//    return vec_to_string(A, stdFmt);
//};
//
/////Templates for turning Eigen matrices into strings
//template <class T>  std::string mat_to_string(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &A, const char *fmt) {
////std::string mat_to_string(const Eigen::MatrixXd &A, const char *fmt) {
//    std::size_t r = A.rows();
//    std::size_t c = A.cols();
//    if ((r<1)||(c<1)) return std::string("");
//    std::stringstream out;
//    out << "[ ";
//    if (r==1) {
//        out << format(fmt, A(0,0));
//        for (size_t j = 1; j < c; j++) {
//            out << ", " << format(fmt, A(0,j));
//        }
//    } else {
//        out << mat_to_string(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>(A.row(0)), fmt);
//        for (size_t i = 1; i < r; i++) {
//            out << ", " << std::endl << "  " << mat_to_string(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>(A.row(i)), fmt);
//        }
//    }
//    out << " ]";
//    return out.str();
//};
//template <class T> std::string mat_to_string(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &A) {
////std::string vec_to_string(const Eigen::MatrixXd &A) {
//    return mat_to_string(A, stdFmt);
//};

/// Template class for turning numbers (0D-matrices) into strings
//template<class T> std::string vec_to_string(const             T  &a){
//    return vec_to_string(a, stdFmt);
//    std::stringstream out;
//    out << format("[ %7.3f ]",a);
//    return out.str();
//};
//template<class T> std::string vec_to_string(const VectorNd<0, T> &a){
//    return vec_to_string(a, stdFmt);
//};
//template<class T> std::string vec_to_string(const VectorNd<0, T> &a, const char *fmt) {
//    VectorNd<1, T> vec;
//    vec.push_back(a);
//    return vec_to_string(vec, fmt);
//};
//
/////Template classes for turning vectors (1D-matrices) into strings
//template<class T> std::string vec_to_string(const VectorNd<1, T> &a) {
//    return vec_to_string(a, stdFmt);
//};
//template<class T> std::string vec_to_string(const VectorNd<1, T> &a, const char *fmt) {
//    if (a.size()<1) {
//        return std::string("");
//    } else {
//        std::stringstream out;
//        out << "[ ";
//        out << format(fmt,a[0]);
//        for (size_t j = 1; j < a.size(); j++) {
//            out << ", ";
//            out << format(fmt,a[j]);
//        }
//        out << " ]";
//        return out.str();
//    }
//};
//
/////Template classes for turning 2D-matrices into strings
//template<class T> std::string vec_to_string(const VectorNd<2, T> &A) {
//    return vec_to_string(A, stdFmt);
//}
//template<class T> std::string vec_to_string(const VectorNd<2, T> &A, const char *fmt) {
//    if (A.size()<1) return std::string("");
//    std::stringstream out;
//    out << "[ " << format(fmt,A[0]);
//    for (size_t j = 1; j < A.size(); j++) {
//        out << ", " << std::endl << "  " << vec_to_string(A[j], fmt);
//    }
//    out << " ]";
//    return out.str();
//}

///// Publish the linear algebra solver
//template<class T> std::vector<T>   linsolve(std::vector<std::vector<T> > const& A,             std::vector<T>   const& b);
//template<class T> std::vector<std::vector<T> > linsolve(std::vector<std::vector<T> > const& A, std::vector<std::vector<T> > const& B);
//
///// Some shortcuts and regularly needed operations
//template<class T> std::size_t         num_rows  (std::vector<std::vector<T> > const& in);
//template<class T> std::size_t         num_cols  (std::vector<std::vector<T> > const& in);
//template<class T> std::size_t         max_cols  (std::vector<std::vector<T> > const& in);
//template<class T> std::vector<T> get_row   (std::vector<std::vector<T> > const& in, size_t row);
//template<class T> std::vector<T> get_col   (std::vector<std::vector<T> > const& in, size_t col);
//template<class T> bool                is_squared(std::vector<std::vector<T> > const& in);
//template<class T> std::vector<std::vector<T> > make_squared(std::vector<std::vector<T> > const& in);
//
///// Define some basic math operations for vectors
//template<class T> T    multiply(            std::vector<T>   const& A,             std::vector<T>   const& B);
//template<class T>             std::vector<T>   multiply(std::vector<std::vector<T> > const& A,             std::vector<T>   const& B);
//template<class T> std::vector<std::vector<T> > multiply(std::vector<std::vector<T> > const& A, std::vector<std::vector<T> > const& B);
//
//template<class T> T              dot_product(std::vector<T> const& a, std::vector<T> const& b);
//template<class T> std::vector<T> cross_product(std::vector<T> const& a, std::vector<T> const& b);
//
//template<class T> std::vector<std::vector<T> > transpose(std::vector<std::vector<T> > const& in);
//template<class T> std::vector<std::vector<T> >    invert(std::vector<std::vector<T> > const& in);
//
//template<class T> std::string vec_to_string(                        T    const& a);
//template<class T> std::string vec_to_string(            std::vector<T>   const& a);
//template<class T> std::string vec_to_string(std::vector<std::vector<T> > const& A);
//
//template<class T> std::string vec_to_string(            std::vector<T>   const& a, const char *fmt);
//template<class T> std::string vec_to_string(std::vector<std::vector<T> > const& A, const char *fmt);

/*
Owe a debt of gratitude to http://sole.ooz.ie/en - very clear treatment of GJ
*/
template <typename T>
void swap_rows(std::vector<std::vector<T>>* A, size_t row1, size_t row2) {
    for (size_t col = 0; col < (*A)[0].size(); col++) {
        std::swap((*A)[row1][col], (*A)[row2][col]);
    }
};
template <typename T>
void subtract_row_multiple(std::vector<std::vector<T>>* A, size_t row, T multiple, size_t pivot_row) {
    for (size_t col = 0; col < (*A)[0].size(); col++) {
        (*A)[row][col] -= multiple * (*A)[pivot_row][col];
    }
};
template <typename T>
void divide_row_by(std::vector<std::vector<T>>* A, size_t row, T value) {
    for (size_t col = 0; col < (*A)[0].size(); col++) {
        (*A)[row][col] /= value;
    }
};

template <typename T>
size_t get_pivot_row(std::vector<std::vector<T>>* A, size_t col) {
    std::size_t index = col;
    T max = 0, val;

    for (size_t row = col; row < (*A).size(); row++) {
        val = (*A)[row][col];
        if (std::abs(val) > max) {
            max = std::abs(val);
            index = row;
        }
    }
    return index;
};

template <typename T>
std::vector<std::vector<T>> linsolve_Gauss_Jordan(std::vector<std::vector<T>> const& A, std::vector<std::vector<T>> const& B) {
    std::vector<std::vector<T>> AB;
    std::vector<std::vector<T>> X;
    size_t pivot_row;
    T pivot_element;

    size_t NrowA = num_rows(A);
    size_t NrowB = num_rows(B);
    size_t NcolA = num_cols(A);
    size_t NcolB = num_cols(B);

    if (NrowA != NrowB) throw ValueError(format("You have to provide matrices with the same number of rows: %d is not %d. ", NrowA, NrowB));

    AB.resize(NrowA, std::vector<T>(NcolA + NcolB, 0));
    X.resize(NrowA, std::vector<T>(NcolB, 0));

    // Build the augmented matrix
    for (size_t row = 0; row < NrowA; row++) {
        for (size_t col = 0; col < NcolA; col++) {
            AB[row][col] = A[row][col];
        }
        for (size_t col = NcolA; col < NcolA + NcolB; col++) {
            AB[row][col] = B[row][col - NcolA];
        }
    }

    for (size_t col = 0; col < NcolA; col++) {
        // Find the pivot value
        pivot_row = get_pivot_row(&AB, col);

        if (std::abs(AB[pivot_row][col]) < 10 * DBL_EPSILON) {
            throw ValueError(format("Zero occurred in row %d, the matrix is singular. ", pivot_row));
        }

        if (pivot_row >= col) {
            // Swap pivot row and current row
            swap_rows(&AB, col, pivot_row);
        }
        // Get the pivot element
        pivot_element = AB[col][col];
        // Divide the pivot row by the pivot element
        divide_row_by(&AB, col, pivot_element);

        if (col < NrowA - 1) {
            // All the rest of the rows, subtract the value of the [r][c] combination
            for (size_t row = col + 1; row < NrowA; row++) {
                subtract_row_multiple(&AB, row, AB[row][col], col);
            }
        }
    }
    for (std::size_t col = NcolA - 1; col > 0; col--) {
        for (int row = static_cast<int>(col) - 1; row >= 0; row--) {
            subtract_row_multiple(&AB, row, AB[row][col], col);
        }
    }
    // Set the output value
    for (size_t row = 0; row < NrowA; row++) {
        for (size_t col = 0; col < NcolB; col++) {
            X[row][col] = AB[row][NcolA + col];
        }
    }
    return X;
};

//std::vector<std::vector<double> > linsolve_Gauss_Jordan_reimpl(std::vector<std::vector<double> > const& A, std::vector<std::vector<double> > const& B) {
//    std::vector<std::vector<double> > AB;
//    std::vector<std::vector<double> > X;
//    size_t pivot_row;
//    double pivot_element;
//    double tmp_element;
//
//    size_t NrowA = num_rows(A);
//    size_t NrowB = num_rows(B);
//    size_t NcolA = num_cols(A);
//    size_t NcolB = num_cols(B);
//
//    if (NrowA!=NrowB) throw ValueError(format("You have to provide matrices with the same number of rows: %d is not %d. ",NrowA,NrowB));
//
//    AB.resize(NrowA, std::vector<double>(NcolA+NcolB, 0));
//     X.resize(NrowA, std::vector<double>(NcolB, 0));
//
//    // Build the augmented matrix
//    for (size_t row = 0; row < NrowA; row++){
//        for (size_t col  = 0; col < NcolA; col++){
//            AB[row][col] = A[row][col];
//        }
//        for (size_t col  = NcolA; col < NcolA+NcolB; col++){
//            AB[row][col] = B[row][col-NcolA];
//        }
//    }
//
//    for (size_t col = 0; col < NcolA; col++){
//        // Find the pivot row
//        pivot_row     = 0;
//        pivot_element = 0.0;
//        for (size_t row = col; row < NrowA; row++){
//            tmp_element = std::abs(AB[row][col]);
//            if (tmp_element>pivot_element) {
//                pivot_element = tmp_element;
//                pivot_row = row;
//            }
//        }
//        // Check for errors
//        if (AB[pivot_row][col]<1./_HUGE) throw ValueError(format("Zero occurred in row %d, the matrix is singular. ",pivot_row));
//        // Swap the rows
//        if (pivot_row>col) {
//            for (size_t colInt = 0; colInt < NcolA; colInt++){
//                std::swap(AB[pivot_row][colInt],AB[pivot_row][colInt]);
//            }
//        }
//        // Process the entries below current element
//        for (size_t row = col; row < NrowA; row++){
//            // Entries to the right of current element (until end of A)
//            for (size_t colInt = col+1; colInt < NcolA; colInt++){
//                // All entries in augmented matrix
//                for (size_t colFull = col; colFull < NcolA+NcolB; colFull++){
//                    AB[colInt][colFull] -= AB[col][colFull] * AB[colInt][col] / AB[col][col];
//                }
//                AB[colInt][col] = 0.0;
//            }
//        }
//    }
//    return AB;
//}

template <class T>
std::vector<std::vector<T>> linsolve(std::vector<std::vector<T>> const& A, std::vector<std::vector<T>> const& B) {
    return linsolve_Gauss_Jordan(A, B);
};

template <class T>
std::vector<T> linsolve(std::vector<std::vector<T>> const& A, std::vector<T> const& b) {
    std::vector<std::vector<T>> B;
    for (size_t i = 0; i < b.size(); i++) {
        B.push_back(std::vector<T>(1, b[i]));
    }
    B = linsolve(A, B);
    B[0].resize(B.size(), 0.0);
    for (size_t i = 1; i < B.size(); i++) {
        B[0][i] = B[i][0];
    }
    return B[0];
};

template <class T>
std::vector<T> get_row(std::vector<std::vector<T>> const& in, size_t row) {
    return in[row];
};
template <class T>
std::vector<T> get_col(std::vector<std::vector<T>> const& in, size_t col) {
    std::size_t sizeX = in.size();
    if (sizeX < 1) throw ValueError(format("You have to provide values, a vector length of %d is not valid. ", sizeX));
    size_t sizeY = in[0].size();
    if (sizeY < 1) throw ValueError(format("You have to provide values, a vector length of %d is not valid. ", sizeY));
    std::vector<T> out;
    for (std::size_t i = 0; i < sizeX; i++) {
        sizeY = in[i].size();
        if (sizeY - 1 < col)
            throw ValueError(format("Your matrix does not have enough entries in row %d, last index %d is less than %d. ", i, sizeY - 1, col));
        out.push_back(in[i][col]);
    }
    return out;
};

template <class T>
std::vector<std::vector<T>> make_squared(std::vector<std::vector<T>> const& in) {
    std::size_t cols = max_cols(in);
    std::size_t rows = num_rows(in);
    std::size_t maxVal = 0;
    std::vector<std::vector<T>> out;
    std::vector<T> tmp;

    if (cols > rows) {
        maxVal = cols;
    } else {
        maxVal = rows;
    }
    out.clear();
    for (std::size_t i = 0; i < in.size(); i++) {
        tmp.clear();
        for (std::size_t j = 0; j < in[i].size(); j++) {
            tmp.push_back(in[i][j]);
        }
        while (maxVal > tmp.size()) {
            tmp.push_back(0.0);
        }
        out.push_back(tmp);
    }
    // Check rows
    tmp.clear();
    tmp.resize(maxVal, 0.0);
    while (maxVal > out.size()) {
        out.push_back(tmp);
    }
    return out;
};

template <class T>
T multiply(std::vector<T> const& a, std::vector<T> const& b) {
    return dot_product(a, b);
};
template <class T>
std::vector<T> multiply(std::vector<std::vector<T>> const& A, std::vector<T> const& b) {
    std::vector<std::vector<T>> B;
    for (size_t i = 0; i < b.size(); i++) {
        B.push_back(std::vector<T>(1, b[i]));
    }
    B = multiply(A, B);
    B[0].resize(B.size(), 0.0);
    for (size_t i = 1; i < B.size(); i++) {
        B[0][i] = B[i][0];
    }
    return B[0];
}

template <class T>
std::vector<std::vector<T>> multiply(std::vector<std::vector<T>> const& A, std::vector<std::vector<T>> const& B) {
    if (num_cols(A) != num_rows(B)) {
        throw ValueError(format("You have to provide matrices with the same columns and rows: %d is not equal to %d. ", num_cols(A), num_rows(B)));
    }
    size_t rows = num_rows(A);
    size_t cols = num_cols(B);
    T tmp;
    std::vector<std::vector<T>> outVec;
    std::vector<T> tmpVec;
    outVec.clear();
    for (size_t i = 0; i < rows; i++) {
        tmpVec.clear();
        for (size_t j = 0; j < cols; j++) {
            tmp = 0.0;
            for (size_t k = 0; k < num_cols(A); k++) {
                tmp += A[i][k] * B[k][j];
            }
            tmpVec.push_back(tmp);
        }
        outVec.push_back(tmpVec);
    }
    return outVec;
};

template <class T>
T dot_product(std::vector<T> const& a, std::vector<T> const& b) {
    if (a.size() == b.size()) {
        return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
    }
    throw ValueError(format("You have to provide vectors with the same length: %d is not equal to %d. ", a.size(), b.size()));
};

template <class T>
std::vector<T> cross_product(std::vector<T> const& a, std::vector<T> const& b) {
    throw NotImplementedError("The cross product function has not been implemented, yet");
};

template <class T>
std::vector<std::vector<T>> transpose(std::vector<std::vector<T>> const& in) {
    size_t sizeX = in.size();
    if (sizeX < 1) throw ValueError(format("You have to provide values, a vector length of %d is not a valid. ", sizeX));
    size_t sizeY = in[0].size();
    size_t sizeYOld = sizeY;
    if (sizeY < 1) throw ValueError(format("You have to provide values, a vector length of %d is not a valid. ", sizeY));
    std::vector<std::vector<T>> out(sizeY, std::vector<T>(sizeX));
    for (size_t i = 0; i < sizeX; ++i) {
        sizeY = in[i].size();
        if (sizeY != sizeYOld) throw ValueError(format("You have to provide a rectangular matrix: %d is not equal to %d. ", sizeY, sizeYOld));
        for (size_t j = 0; j < sizeY; ++j) {
            out[j][i] = in[i][j];
        }
    }
    return out;
};

template <class T>
std::vector<std::vector<T>> invert(std::vector<std::vector<T>> const& in) {
    if (!is_squared(in)) throw ValueError(format("Only square matrices can be inverted: %d is not equal to %d. ", num_rows(in), num_cols(in)));
    std::vector<std::vector<T>> identity;
    // Build the identity matrix
    size_t dim = num_rows(in);
    identity.resize(dim, std::vector<T>(dim, 0));
    for (size_t row = 0; row < dim; row++) {
        identity[row][row] = 1.0;
    }
    return linsolve(in, identity);
};

inline void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove) {
    unsigned int numRows = static_cast<unsigned int>(matrix.rows()) - 1;
    unsigned int numCols = static_cast<unsigned int>(matrix.cols());

    if (rowToRemove <= numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
    else {
        throw ValueError(format("Trying to remove row index [%d] greater than max index [%d] ", rowToRemove, numRows));
    }
    matrix.conservativeResize(numRows, numCols);
};

inline void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove) {
    unsigned int numRows = static_cast<unsigned int>(matrix.rows());
    unsigned int numCols = static_cast<unsigned int>(matrix.cols()) - 1;

    if (colToRemove <= numCols)
        matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
    else {
        throw ValueError(format("Trying to remove column index [%d] greater than max index [%d] ", colToRemove, numCols));
    }
    matrix.conservativeResize(numRows, numCols);
};
template <typename Derived>
inline Eigen::MatrixXd minor_matrix(const Eigen::MatrixBase<Derived>& A, std::size_t i, std::size_t j) {
    Eigen::MatrixXd Am = A;
    removeRow(Am, static_cast<unsigned int>(i));
    removeColumn(Am, static_cast<unsigned int>(j));
    return Am;
};

template <typename Derived>
static Eigen::MatrixXd adjugate(const Eigen::MatrixBase<Derived>& A) {
    std::size_t N = A.rows();
    if (N == 1) {
        Eigen::MatrixXd Aadj(1, 1);
        Aadj << 1;
        return Aadj;
    }
    Eigen::MatrixXd Aadj(N, N);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            int negative_1_to_the_i_plus_j = ((i + j) % 2 == 0) ? 1 : -1;
            Aadj(i, j) = negative_1_to_the_i_plus_j * minor_matrix(A, j, i).determinant();
        }
    }
    return Aadj;
}

template <typename Derived>
static Eigen::MatrixXd adjugate_derivative(const Eigen::MatrixBase<Derived>& A, const Eigen::MatrixBase<Derived>& dAdt) {
    std::size_t N = A.rows();
    Eigen::MatrixXd Aadj(N, N);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            int negative_1_to_the_i_plus_j = ((i + j) % 2 == 0) ? 1 : -1;
            Eigen::MatrixXd mm = minor_matrix(A, j, i);
            Aadj(i, j) = negative_1_to_the_i_plus_j * (adjugate(minor_matrix(A, j, i)) * minor_matrix(dAdt, j, i)).trace();
        }
    }
    return Aadj;
}

}; /* namespace CoolProp */
#endif
