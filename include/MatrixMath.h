#ifndef MATRIXMATH_H
#define MATRIXMATH_H

#include "CoolPropTools.h"
#include "Exceptions.h"

#include <vector>
#include <string>
#include <numeric> // inner_product
#include <sstream>
#include "float.h"

namespace CoolProp{

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

// Forward definitions
template<class T> std::size_t         num_rows  (std::vector<std::vector<T> > const& in);
template<class T> std::size_t         max_cols  (std::vector<std::vector<T> > const& in);

template<class T> bool                is_squared(std::vector<std::vector<T> > const& in){
	std::size_t cols = max_cols(in);
	if (cols!=num_rows(in)) { return false;}
	else {
		for (std::size_t i = 0; i < in.size(); i++) {
			if (cols!=in[i].size()) {return false; }
		}
	}
	return true;
};
template<class T> std::size_t         max_cols  (std::vector<std::vector<T> > const& in){
	std::size_t cols = 0;
	std::size_t col  = 0;
	for (std::size_t i = 0; i < in.size(); i++) {
		col = in[i].size();
		if (cols<col) {cols = col;}
    }
	return cols;
};
/// Some shortcuts and regularly needed operations
template<class T> std::size_t         num_rows  (std::vector<std::vector<T> > const& in){ return in.size(); }
template<class T> std::size_t         num_cols  (std::vector<std::vector<T> > const& in){
	if (num_rows(in)>0) {
		if (is_squared(in)) {
			return in[0].size();
		} else {
			return max_cols(in);
		}
	} else {
		return 0;
	}
};


/*
Owe a debt of gratitude to http://sole.ooz.ie/en - very clear treatment of GJ
*/
template<typename T> void swap_rows(std::vector<std::vector<T> > *A, size_t row1, size_t row2)
{
	for (size_t col = 0; col < (*A)[0].size(); col++){
		std::swap((*A)[row1][col],(*A)[row2][col]);
	}
};
template<typename T> void subtract_row_multiple(std::vector<std::vector<T> > *A, size_t row, T multiple, size_t pivot_row)
{
	for (size_t col = 0; col < (*A)[0].size(); col++){
		(*A)[row][col] -= multiple*(*A)[pivot_row][col];
	}
};
template<typename T> void divide_row_by(std::vector<std::vector<T> > *A, size_t row, T value)
{
	for (size_t col = 0; col < (*A)[0].size(); col++){
		(*A)[row][col] /= value;
	}
};

template<typename T> size_t get_pivot_row(std::vector<std::vector<T> > *A, size_t col)
{
	int index = col;
	T max = 0, val;

	for (size_t row = col; row < (*A).size(); row++)
	{
		val = (*A)[row][col];
		if (fabs(val) > max)
		{
			max = fabs(val);
			index = row;
		}
	}
	return index;
};


template<typename T> std::vector<std::vector<T> > linsolve_Gauss_Jordan(std::vector<std::vector<T> > const& A, std::vector<std::vector<T> > const& B) {
	std::vector<std::vector<T> > AB;
	std::vector<std::vector<T> > X;
	size_t pivot_row;
	T pivot_element;

	size_t NrowA = num_rows(A);
	size_t NrowB = num_rows(B);
	size_t NcolA = num_cols(A);
	size_t NcolB = num_cols(B);

	if (NrowA!=NrowB) throw ValueError(format("You have to provide matrices with the same number of rows: %d is not %d. ",NrowA,NrowB));

	AB.resize(NrowA, std::vector<T>(NcolA+NcolB, 0));
	 X.resize(NrowA, std::vector<T>(NcolB, 0));

	// Build the augmented matrix
	for (size_t row = 0; row < NrowA; row++){
		for (size_t col  = 0; col < NcolA; col++){
			AB[row][col] = A[row][col];
		}
		for (size_t col  = NcolA; col < NcolA+NcolB; col++){
			AB[row][col] = B[row][col-NcolA];
		}
	}

	for (size_t col = 0; col < NcolA; col++){
		// Find the pivot value
		pivot_row = get_pivot_row(&AB, col);

		if (fabs(AB[pivot_row][col]) < 10*DBL_EPSILON){ throw ValueError(format("Zero occurred in row %d, the matrix is singular. ",pivot_row));}

		if (pivot_row>=col){
			// Swap pivot row and current row
			swap_rows(&AB, col, pivot_row);
		}
		// Get the pivot element
		pivot_element = AB[col][col];
		// Divide the pivot row by the pivot element
		divide_row_by(&AB,col,pivot_element);

		if (col < NrowA-1)
		{
			// All the rest of the rows, subtract the value of the [r][c] combination
			for (size_t row = col + 1; row < NrowA; row++)
			{
				subtract_row_multiple(&AB,row,AB[row][col],col);
			}
		}
	}
	for (int col = NcolA - 1; col > 0; col--)
	{
		for (int row = col - 1; row >=0; row--)
		{
			subtract_row_multiple(&AB,row,AB[row][col],col);
		}
	}
	// Set the output value
	for (size_t row = 0; row < NrowA; row++){
		for (size_t col  = 0; col < NcolB; col++){
			X[row][col] = AB[row][NcolA+col];
		}
	}
	return X;
};


//std::vector<std::vector<double> > linsolve_Gauss_Jordan_reimpl(std::vector<std::vector<double> > const& A, std::vector<std::vector<double> > const& B) {
//	std::vector<std::vector<double> > AB;
//	std::vector<std::vector<double> > X;
//	size_t pivot_row;
//	double pivot_element;
//	double tmp_element;
//
//	size_t NrowA = num_rows(A);
//	size_t NrowB = num_rows(B);
//	size_t NcolA = num_cols(A);
//	size_t NcolB = num_cols(B);
//
//	if (NrowA!=NrowB) throw ValueError(format("You have to provide matrices with the same number of rows: %d is not %d. ",NrowA,NrowB));
//
//	AB.resize(NrowA, std::vector<double>(NcolA+NcolB, 0));
//	 X.resize(NrowA, std::vector<double>(NcolB, 0));
//
//	// Build the augmented matrix
//	for (size_t row = 0; row < NrowA; row++){
//		for (size_t col  = 0; col < NcolA; col++){
//			AB[row][col] = A[row][col];
//		}
//		for (size_t col  = NcolA; col < NcolA+NcolB; col++){
//			AB[row][col] = B[row][col-NcolA];
//		}
//	}
//
//	for (size_t col = 0; col < NcolA; col++){
//		// Find the pivot row
//		pivot_row     = 0;
//		pivot_element = 0.0;
//		for (size_t row = col; row < NrowA; row++){
//			tmp_element = fabs(AB[row][col]);
//			if (tmp_element>pivot_element) {
//				pivot_element = tmp_element;
//				pivot_row = row;
//			}
//		}
//		// Check for errors
//		if (AB[pivot_row][col]<1./_HUGE) throw ValueError(format("Zero occurred in row %d, the matrix is singular. ",pivot_row));
//		// Swap the rows
//		if (pivot_row>col) {
//			for (size_t colInt = 0; colInt < NcolA; colInt++){
//				std::swap(AB[pivot_row][colInt],AB[pivot_row][colInt]);
//			}
//		}
//		// Process the entries below current element
//		for (size_t row = col; row < NrowA; row++){
//			// Entries to the right of current element (until end of A)
//			for (size_t colInt = col+1; colInt < NcolA; colInt++){
//				// All entries in augmented matrix
//				for (size_t colFull = col; colFull < NcolA+NcolB; colFull++){
//					AB[colInt][colFull] -= AB[col][colFull] * AB[colInt][col] / AB[col][col];
//				}
//				AB[colInt][col] = 0.0;
//			}
//		}
//	}
//	return AB;
//}






template<class T> std::vector<std::vector<T> > linsolve(std::vector<std::vector<T> > const& A, std::vector<std::vector<T> > const& B){
	return linsolve_Gauss_Jordan(A, B);
};

template<class T> std::vector<T>               linsolve(std::vector<std::vector<T> > const& A,             std::vector<T>   const& b){
	std::vector<std::vector<T> > B;
	for (size_t i = 0; i < b.size(); i++){
		B.push_back(std::vector<T>(1,b[i]));
	}
	B = linsolve(A, B);
	B[0].resize(B.size(),0.0);
	for (size_t i = 1; i < B.size(); i++){
		B[0][i] = B[i][0];
	}
	return B[0];
};



template<class T> std::vector<T> get_row(std::vector< std::vector<T> > const& in, size_t row) { return in[row]; };
template<class T> std::vector<T> get_col(std::vector< std::vector<T> > const& in, size_t col) {
	std::size_t sizeX  = in.size();
	if (sizeX<1) throw ValueError(format("You have to provide values, a vector length of %d is not valid. ",sizeX));
	size_t sizeY = in[0].size();
	if (sizeY<1) throw ValueError(format("You have to provide values, a vector length of %d is not valid. ",sizeY));
	std::vector<T> out;
	for (std::size_t i = 0; i < sizeX; i++) {
		sizeY = in[i].size();
		if (sizeY-1<col) throw ValueError(format("Your matrix does not have enough entries in row %d, last index %d is less than %d. ",i,sizeY-1,col));
		out.push_back(in[i][col]);
	}
	return out;
};

template<class T> std::vector<std::vector<T> > make_squared(std::vector<std::vector<T> > const& in){
	std::size_t cols   = max_cols(in);
	std::size_t rows   = num_rows(in);
	std::size_t maxVal = 0;
	std::vector<std::vector<T> > out;
	            std::vector<T>   tmp;

	if (cols>rows) {maxVal = cols; }
	else           {maxVal = rows; }
	out.clear();
	for (std::size_t i = 0; i < in.size(); i++) {
		tmp.clear();
		for (std::size_t j = 0; j < in[i].size(); j++) {
			tmp.push_back(in[i][j]);
		}
		while (maxVal>tmp.size()) {
			tmp.push_back(0.0);
		}
		out.push_back(tmp);
    }
	// Check rows
	tmp.clear();
	tmp.resize(maxVal,0.0);
	while (maxVal>out.size()) {
		out.push_back(tmp);
	}
	return out;
};

template<class T>                         T    multiply(            std::vector<T>   const& a,             std::vector<T>   const& b){
    return dot_product(a,b);

};
template<class T>             std::vector<T>   multiply(std::vector<std::vector<T> > const& A,             std::vector<T>   const& b){
	std::vector<std::vector<T> > B;
	for (size_t i = 0; i < b.size(); i++){
		B.push_back(std::vector<T>(1,b[i]));
	}
	B = multiply(A, B);
	B[0].resize(B.size(),0.0);
	for (size_t i = 1; i < B.size(); i++){
		B[0][i] = B[i][0];
	}
	return B[0];
}

template<class T> std::vector<std::vector<T> > multiply(std::vector<std::vector<T> > const& A, std::vector<std::vector<T> > const& B){
	if (num_cols(A) != num_rows(B)){
		throw ValueError(format("You have to provide matrices with the same columns and rows: %d is not equal to %d. ",num_cols(A),num_rows(B)));
	}
	size_t rows = num_rows(A);
	size_t cols = num_cols(B);
	T tmp;
	std::vector<std::vector<T> > outVec;
		        std::vector<T>   tmpVec;
    outVec.clear();
	for (size_t i = 0; i < rows; i++){
		tmpVec.clear();
		for (size_t j = 0; j < cols; j++){
			tmp = 0.0;
			for (size_t k = 0; k < num_cols(A); k++){
				tmp += A[i][k] * B[k][j];
			}
			tmpVec.push_back(tmp);
		}
		outVec.push_back(tmpVec);
	}
	return outVec;
};

template<class T> T              dot_product(std::vector<T> const& a, std::vector<T> const& b){
	if (a.size()==b.size()){
		return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
	}
	throw ValueError(format("You have to provide vectors with the same length: %d is not equal to %d. ",a.size(),b.size()));
};

template<class T> std::vector<T> cross_product(std::vector<T> const& a, std::vector<T> const& b){
	throw NotImplementedError("The cross product function has not been implemented, yet");
};

template<class T> std::vector< std::vector<T> > transpose(std::vector<std::vector<T> > const& in){
	size_t sizeX = in.size();
	if (sizeX<1) throw ValueError(format("You have to provide values, a vector length of %d is not a valid. ",sizeX));
	size_t sizeY    = in[0].size();
	size_t sizeYOld = sizeY;
	if (sizeY<1) throw ValueError(format("You have to provide values, a vector length of %d is not a valid. ",sizeY));
	std::vector< std::vector<T> > out(sizeY,std::vector<T>(sizeX));
	for (size_t i = 0; i < sizeX; ++i){
		sizeY = in[i].size();
		if (sizeY!=sizeYOld) throw ValueError(format("You have to provide a rectangular matrix: %d is not equal to %d. ",sizeY,sizeYOld));
		for (size_t j = 0; j < sizeY; ++j){
			out[j][i] = in[i][j];
		}
	}
	return out;
};

template<class T> std::vector< std::vector<T> >    invert(std::vector<std::vector<T> > const& in){
	if (!is_squared(in)) throw ValueError(format("Only square matrices can be inverted: %d is not equal to %d. ",num_rows(in),num_cols(in)));
	std::vector<std::vector<T> > identity;
	// Build the identity matrix
	size_t dim = num_rows(in);
	identity.resize(dim, std::vector<T>(dim, 0));
	for (size_t row = 0; row < dim; row++){
		identity[row][row] = 1.0;
	}
	return linsolve(in,identity);
};

template<class T> std::string vec_to_string(                        T    const& a){
	std::stringstream out;
	out << format("[ %7.3f ]",a);
	return out.str();
};

template<class T> std::string vec_to_string(            std::vector<T>   const& a) {
	return vec_to_string(a,"%7.3g");
};
template<class T> std::string vec_to_string(            std::vector<T>   const& a, const char *fmt) {
	if (a.size()<1) {
		return std::string("");
	} else {
		std::stringstream out;
		out << format("[ ");
		out << format(fmt,a[0]);
		for (size_t j = 1; j < a.size(); j++) {
			out << ", ";
			out << format(fmt,a[j]);
		}
		out << " ]";
		return out.str();
	}
};

template<class T> std::string vec_to_string(std::vector<std::vector<T> > const& A) {
	return vec_to_string(A, "%7.3g");
}

template<class T> std::string vec_to_string(std::vector<std::vector<T> > const& A, const char *fmt) {
	std::stringstream out;
	for (size_t j = 0; j < A.size(); j++) {
		out << vec_to_string(A[j], fmt);
    }
    return out.str();
}

}; /* namespace CoolProp */
#endif
