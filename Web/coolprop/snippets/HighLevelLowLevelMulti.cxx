#include "CoolPropLib.h"
#include "CoolPropTools.h"
#include <vector>
#include <time.h>

int main() {
    const long buffer_size = 1000, length = 100000;
    long ierr;
    char herr[buffer_size];
    long handle = AbstractState_factory("BICUBIC&HEOS", "Water", &ierr, herr, buffer_size);
    std::vector<double> T(length), p(length), rhomolar(length), hmolar(length), smolar(length);
    std::vector<double> input1 = linspace(700000.0, 1500000.0, length);
    std::vector<double> input2 = linspace(2.8e6, 3.0e6, length);
    long input_pair = get_input_pair_index("HmassP_INPUTS");
    double t1 = clock();
    AbstractState_update_and_common_out(handle, input_pair, &(input1[0]), &(input2[0]), length, &(T[0]), &(p[0]), &(rhomolar[0]), &(hmolar[0]),
                                        &(smolar[0]), &ierr, herr, buffer_size);
    double t2 = clock();
    std::cout << format("value(commons): %g us/call\n", ((double)(t2 - t1)) / CLOCKS_PER_SEC / double(length) * 1e6);

    std::vector<long> outputs(5);
    outputs[0] = get_param_index("T");
    outputs[1] = get_param_index("P");
    outputs[2] = get_param_index("Dmolar");
    outputs[3] = get_param_index("Hmolar");
    outputs[4] = get_param_index("Smolar");
    std::vector<double> out1(length), out2(length), out3(length), out4(length), out5(length);
    t1 = clock();
    AbstractState_update_and_5_out(handle, input_pair, &(input1[0]), &(input2[0]), length, &(outputs[0]), &(out1[0]), &(out2[0]), &(out3[0]),
                                   &(out4[0]), &(out5[0]), &ierr, herr, buffer_size);
    t2 = clock();
    std::cout << format("value(user-specified): %g us/call\n", ((double)(t2 - t1)) / CLOCKS_PER_SEC / double(length) * 1e6);
}