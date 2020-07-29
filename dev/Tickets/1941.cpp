
#include <vector>
#include <string>
#include <chrono>

#include "CoolProp.h"

int main(int argc, const char* argv[]) {

    //double v_ref = 1.0/CoolProp::PropsSI("Dmass", "P", 72e5, "T", 31+273.15, "HEOS::CarbonDioxide");

    std::vector<std::string> Outputs;
    std::string Name1;
    std::vector<double> Prop1;
    std::string Name2;
    std::vector<double> Prop2;
    std::string backend;
    std::vector<std::string> fluids;
    std::vector<double> fractions;

    std::vector<std::vector<double> > out = CoolProp::PropsSImulti(Outputs, Name1, Prop1, Name2, Prop2, backend, fluids, fractions);


    return 0;
}
