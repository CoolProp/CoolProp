#include <unsupported/Eigen/Polynomials>
#include <iostream>
//using namespace Eigen;
//using namespace std;

#include <vector>
#include <string>
#include <MatrixMath.h>

int main()
{
Eigen::Vector4d roots = Eigen::Vector4d::Random();
std::cout << "Roots: " << roots.transpose() << std::endl;
Eigen::Matrix<double,5,1> polynomial;
Eigen::roots_to_monicPolynomial( roots, polynomial );
std::cout << "Polynomial: ";
for( int i=0; i<4; ++i ){ std::cout << polynomial[i] << ".x^" << i << "+ "; }
std::cout << polynomial[4] << ".x^4" << std::endl;
Eigen::Vector4d evaluation;
for( int i=0; i<4; ++i ){
evaluation[i] = Eigen::poly_eval( polynomial, roots[i] ); }
std::cout << "Evaluation of the polynomial at the roots: " << evaluation.transpose() << std::endl;
std::cout << std::endl;
//
//Eigen::MatrixXd coeffs = Eigen::MatrixXd::Random(5,1);
//Eigen::MatrixXd input  = Eigen::MatrixXd::Random(2,1)*1e0;
Eigen::Vector4d coeffs = Eigen::Vector4d::Random()*1e2;
double input  = 1.9e0;
std::cout << "Coeffs: " << std::endl << coeffs.transpose() << std::endl;
double eval = Eigen::poly_eval( coeffs, input);
std::cout << "Evaluation of the polynomial at " << input << std::endl;
std::cout << eval << std::endl;

double vec0 = 0.1;
std::vector<double> vec1(2,0.2);
std::vector< std::vector<double> > vec2;
vec2.push_back(std::vector<double>(2,0.2));
vec2.push_back(std::vector<double>(2,0.3));

std::cout << CoolProp::vec_to_string(vec0) << std::endl;
std::cout << CoolProp::vec_to_string(vec1) << std::endl;
std::cout << CoolProp::vec_to_string(vec2) << std::endl;

}
