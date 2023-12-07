#include "CoolProp.h"
#include <iostream>
#include <stdlib.h>
using namespace CoolProp;
int main() {

    std::vector<std::string> fluids = {"CarbonDioxide", "Nitrogen"};
    std::vector<double> z = {0.97,0.03};
    double temp = 264.65; // 263.15 is fine
    double pres = 6.51e6;
    std::vector<std::string> outputs = {"DMASS", "H", "Z", "PHASE", "VISCOSITY"};

    std::vector<double> T(1, temp), p(1, pres);
    std::vector<std::vector<double>> coolsol =PropsSImulti(outputs, "T", T, "P", p, "HEOS", fluids, z);
    for(int i=0;i<size(outputs);i++)
    {
        std::cout << outputs[i] <<" "<< coolsol[0][i] << std::endl;
    }
    // All done return
    return EXIT_SUCCESS;
}
