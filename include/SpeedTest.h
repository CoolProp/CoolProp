#ifndef SPEEDTEST_H
#define SPEEDTEST_H

#include <string>

namespace CoolProp{

void compare_REFPROP_and_CoolProp(std::string fluid, int inputs, double val1, double val2, std::size_t N);

} /* namespace CoolProp */

#endif