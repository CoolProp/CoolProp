
#include "MatrixMath.h"

#include "CoolPropTools.h"
#include "Exceptions.h"

#include <string>
#include <sstream>
#include <vector>
#include <numeric>
#include <math.h>

namespace CoolProp{


}; /* namespace CoolProp */




#ifdef ENABLE_CATCH
#include <math.h>
#include "catch.hpp"

TEST_CASE("Internal consistency checks and example use cases for MatrixMath.h","[MatrixMath]")
{
	/// Test case for "SylthermXLT" by "Dow Chemicals"
	std::vector<double> cHeat;
	cHeat.clear();
	cHeat.push_back(+1.1562261074E+03);
	cHeat.push_back(+2.0994549103E+00);
	cHeat.push_back(+7.7175381057E-07);
	cHeat.push_back(-3.7008444051E-20);

	SECTION("Eigen::Vector from std::vector") {
		Eigen::Matrix<double,2,1> matrix;
		CoolProp::convert(cHeat, matrix);
	}
}

#endif /* CATCH_ENABLED */
