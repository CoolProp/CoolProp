
#include "MatrixMath.h"

#include "CoolPropTools.h"
#include "Exceptions.h"

#include <string>
#include <sstream>
#include <vector>
#include <numeric>
#include <math.h>

namespace CoolProp {}; /* namespace CoolProp */

#ifdef ENABLE_CATCH
#    include <math.h>
#    include <iostream>
#    include <catch2/catch_all.hpp>

TEST_CASE("Internal consistency checks and example use cases for MatrixMath.h", "[MatrixMath]") {
    bool PRINT = false;

    /// Test case for "SylthermXLT" by "Dow Chemicals"
    std::vector<double> cHeat;
    cHeat.clear();
    cHeat.push_back(+1.1562261074E+03);
    cHeat.push_back(+2.0994549103E+00);
    cHeat.push_back(+7.7175381057E-07);
    cHeat.push_back(-3.7008444051E-20);

    std::vector<std::vector<double>> cHeat2D;
    cHeat2D.push_back(cHeat);
    cHeat2D.push_back(cHeat);

    SECTION("Pretty printing tests") {

        Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(4, 1);
        std::string tmpStr;
        if (PRINT) std::cout << std::endl;

        CHECK_NOTHROW(tmpStr = CoolProp::vec_to_string(cHeat[0]));
        if (PRINT) std::cout << tmpStr << std::endl;
        CHECK_NOTHROW(tmpStr = CoolProp::vec_to_string(cHeat));
        if (PRINT) std::cout << tmpStr << std::endl;
        CHECK_NOTHROW(tmpStr = CoolProp::vec_to_string(cHeat2D));
        if (PRINT) std::cout << tmpStr << std::endl;

        CHECK_NOTHROW(tmpStr = CoolProp::mat_to_string(CoolProp::vec_to_eigen(cHeat[0])));
        if (PRINT) std::cout << tmpStr << std::endl;
        CHECK_NOTHROW(tmpStr = CoolProp::mat_to_string(CoolProp::vec_to_eigen(cHeat, 1)));
        if (PRINT) std::cout << tmpStr << std::endl;
        CHECK_THROWS(tmpStr = CoolProp::mat_to_string(CoolProp::vec_to_eigen(cHeat, 2)));
        if (PRINT) std::cout << tmpStr << std::endl;
        CHECK_NOTHROW(tmpStr = CoolProp::mat_to_string(CoolProp::vec_to_eigen(cHeat2D)));
        if (PRINT) std::cout << tmpStr << std::endl;
    }

    SECTION("Matrix modifications") {
        Eigen::MatrixXd matrix = CoolProp::vec_to_eigen(cHeat2D);

        if (PRINT) std::cout << CoolProp::mat_to_string(matrix) << std::endl;

        CHECK_NOTHROW(CoolProp::removeColumn(matrix, 1));
        if (PRINT) std::cout << CoolProp::mat_to_string(matrix) << std::endl;

        CHECK_NOTHROW(CoolProp::removeRow(matrix, 1));
        if (PRINT) std::cout << CoolProp::mat_to_string(matrix) << std::endl;

        CHECK_THROWS(CoolProp::removeColumn(matrix, 10));
        CHECK_THROWS(CoolProp::removeRow(matrix, 10));
    }

    SECTION("std::vector to Eigen::Matrix and back") {
        std::vector<std::vector<double>> vec2D(cHeat2D);
        Eigen::MatrixXd matrix = CoolProp::vec_to_eigen(vec2D);
        for (size_t i = 0; i < static_cast<size_t>(matrix.cols()); ++i) {
            for (size_t j = 0; j < static_cast<size_t>(matrix.rows()); ++j) {
                CHECK(std::abs(matrix(j, i) - vec2D[j][i]) <= 1e-10);
            }
        }
        vec2D = CoolProp::eigen_to_vec(matrix);
        for (size_t i = 0; i < static_cast<size_t>(matrix.cols()); ++i) {
            for (size_t j = 0; j < static_cast<size_t>(matrix.rows()); ++j) {
                CHECK(std::abs(matrix(j, i) - vec2D[j][i]) <= 1e-10);
            }
        }
        std::vector<double> vec1D(cHeat);
        matrix = CoolProp::vec_to_eigen(vec1D);
        for (size_t i = 0; i < static_cast<size_t>(matrix.cols()); ++i) {
            for (size_t j = 0; j < static_cast<size_t>(matrix.rows()); ++j) {
                CHECK(std::abs(matrix(j, i) - vec1D[j]) <= 1e-10);
            }
        }
        vec1D = CoolProp::eigen_to_vec1D(matrix);
        for (size_t i = 0; i < static_cast<size_t>(matrix.cols()); ++i) {
            for (size_t j = 0; j < static_cast<size_t>(matrix.rows()); ++j) {
                CHECK(std::abs(matrix(j, i) - vec1D[j]) <= 1e-10);
            }
        }
    }
}

#endif /* ENABLE_CATCH */

//#include <unsupported/Eigen/Polynomials>
//#include <iostream>
////using namespace Eigen;
////using namespace std;
//
//#include <vector>
//#include <string>
//#include <MatrixMath.h>
//
//int main()
//{
//Eigen::Vector4d roots = Eigen::Vector4d::Random();
//std::cout << "Roots: " << roots.transpose() << std::endl;
//Eigen::Matrix<double,5,1> polynomial;
//Eigen::roots_to_monicPolynomial( roots, polynomial );
//std::cout << "Polynomial: ";
//for( int i=0; i<4; ++i ){ std::cout << polynomial[i] << ".x^" << i << "+ "; }
//std::cout << polynomial[4] << ".x^4" << std::endl;
//Eigen::Vector4d evaluation;
//for( int i=0; i<4; ++i ){
//evaluation[i] = Eigen::poly_eval( polynomial, roots[i] ); }
//std::cout << "Evaluation of the polynomial at the roots: " << evaluation.transpose() << std::endl;
//std::cout << std::endl;
////
////Eigen::MatrixXd coeffs = Eigen::MatrixXd::Random(5,1);
////Eigen::MatrixXd input  = Eigen::MatrixXd::Random(2,1)*1e0;
//Eigen::Vector4d coeffs = Eigen::Vector4d::Random()*1e2;
//double input  = 1.9e0;
//std::cout << "Coeffs: " << std::endl << coeffs.transpose() << std::endl;
//double eval = Eigen::poly_eval( coeffs, input);
//std::cout << "Evaluation of the polynomial at " << input << std::endl;
//std::cout << eval << std::endl;
//
//double vec0 = 0.1;
//std::vector<double> vec1(2,0.44);
//std::vector< std::vector<double> > vec2;
//vec2.push_back(std::vector<double>(2,0.2));
//vec2.push_back(std::vector<double>(2,0.3));
//
//std::cout << CoolProp::vec_to_string(vec0) << std::endl;
//std::cout << CoolProp::vec_to_string(vec1) << std::endl;
//std::cout << CoolProp::vec_to_string(vec2) << std::endl;
//
//Eigen::Matrix<double,2,2> mat;
//mat.setConstant(2,2,0.25);
//std::vector< std::vector<double> > vec;
//
//CoolProp::convert(mat, vec);
//std::cout << CoolProp::vec_to_string(vec) << std::endl;
//
////Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> mat;
////mat.resize(6,2);
//
//Eigen::Matrix<double,2,2> mat2;
//CoolProp::convert(vec2, mat2);
//CoolProp::convert(mat2, vec);
//std::cout << CoolProp::vec_to_string(vec) << std::endl;
//
//Eigen::Matrix<double,2,1> mat1;
//CoolProp::convert(vec1, mat1);
//std::vector<double> vec3;
//CoolProp::convert(mat1, vec);
//std::cout << CoolProp::vec_to_string(vec) << std::endl;
//
////std::vector< std::vector<double> > vec(vec2);
////CoolProp::convert(mat,vec);
//
////std::cout << CoolProp::vec_to_string() << std::endl;
//
////Eigen::Matrix2d mat2 = CoolProp::convert(vec2);
//
////Eigen::MatrixXd mat2(10,10);
////CoolProp::convert(vec2, mat2);
//
////std::cout << CoolProp::vec_to_string(CoolProp::convert(mat2)) << std::endl;
//
//}
