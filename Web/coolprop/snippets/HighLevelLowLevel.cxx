#include "CoolPropLib.h"
#include "CoolPropTools.h"
#include <vector>
#include <time.h>

int main() {
    double t1, t2;
    const long buffersize = 500;
    long errcode = 0;
    char buffer[buffersize];
    long handle = AbstractState_factory("BICUBIC&HEOS", "Water", &errcode, buffer, buffersize);
    long _HmassP = get_input_pair_index("HmassP_INPUTS");
    long _Dmass = get_param_index("Dmass");
    long len = 20000;
    std::vector<double> h = linspace(700000.0, 1500000.0, len);
    std::vector<double> p = linspace(2.8e6, 3.0e6, len);
    double summer = 0;
    t1 = clock();
    for (long i = 0; i < len; ++i) {
        AbstractState_update(handle, _HmassP, h[i], p[i], &errcode, buffer, buffersize);
        summer += AbstractState_keyed_output(handle, _Dmass, &errcode, buffer, buffersize);
    }
    t2 = clock();
    std::cout << format("value(all): %0.13g, %g us/call\n", summer, ((double)(t2 - t1)) / CLOCKS_PER_SEC / double(len) * 1e6);
    return EXIT_SUCCESS;
}