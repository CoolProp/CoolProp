
#include <vector>
#include <string>

#include "CoolProp.h"
#include "HumidAirProp.h"

//std::map<std::string, double> generate_values(double T, double R, double P = 101325) {
//    double psi_w = HumidAir::HAPropsSI("psi_w", "T", T, "R", R, "P", P);
//    std::vector<std::string> other_output_keys = { "T_wb","T_dp","Hda","Sda","Vda","Omega" };
//    std::map<std::string, double> outputs;
//    outputs.insert_or_assign("psi_w", psi_w);
//    outputs.insert_or_assign("T", T);
//    outputs.insert_or_assign("P", P);
//    outputs.insert_or_assign("R", R);
//    for (auto k : other_output_keys) {
//        outputs.insert_or_assign(k, HumidAir::HAPropsSI(k, "T", T, "R", R, "P", P));
//    }
//    return outputs;
//}

std::vector<std::pair<std::string, double>> generate_values(double T, double R, double P = 101325) {
    double psi_w = HumidAir::HAPropsSI("psi_w", "T", T, "R", R, "P", P);
    std::vector<std::string> other_output_keys = { "T_wb","T_dp","Hda","Sda","Vda","Omega" };
    std::vector<std::pair<std::string, double>> outputs;
    outputs.push_back(std::pair<std::string, double>("psi_w", psi_w));
    outputs.push_back(std::pair<std::string, double>("T", T));
    outputs.push_back(std::pair<std::string, double>("P", P));
    outputs.push_back(std::pair<std::string, double>("R", R));
    for (auto k : other_output_keys) {
        outputs.push_back(std::pair<std::string, double>(k, HumidAir::HAPropsSI(k, "T", T, "R", R, "P", P)));
    }
    return outputs;
}

std::vector<std::pair<std::string, std::string>> get_supported_input_pairs() {
    std::vector<std::pair<std::string, std::string>> good_ones;
    auto inputs = generate_values(300, 0.5);
    std::string k1, k2;
    double v1 = -_HUGE, v2 = -_HUGE, p = -_HUGE;
    for (std::size_t i = 0; i < inputs.size(); i++) {
        k1 = inputs[i].first;
        if (k1.compare("P") == 0) {
            p = inputs[i].second;
        }
    }
    for (std::size_t i = 0; i < inputs.size(); i++) {
        for (std::size_t j = 0; j < inputs.size(); j++) {
            k1 = inputs[i].first;
            k2 = inputs[j].first;
            if (k1.compare("P") == 0 || k2.compare("P") == 0 || k1.compare(k2) == 0) {
                continue;
            }
            v1 = inputs[i].second;
            v2 = inputs[j].second;
            try {
                double psi_w_new = HumidAir::HAPropsSI(
                    "psi_w",
                    k1, v1,
                    k2, v2,
                    "P", p);
                good_ones.push_back(std::pair<std::string, std::string>(k1, k2));
            }
            catch (...) {
                double _nu = -_HUGE;
            }
        }
    }
    return good_ones;
}

void calculate(std::vector<std::pair<std::string, double>> inputs) {
    //auto errors = []

    auto supported_pairs = get_supported_input_pairs();
    std::string k1, k2;
    double psi_w_input = -_HUGE, P_input = -_HUGE, v1 = -_HUGE, v2 = -_HUGE;

    for (std::size_t i = 0; i < inputs.size(); i++) {
        k1 = inputs[i].first;
        if (k1.compare("psi_w") == 0) {
            psi_w_input = inputs[i].second;
        }
        if (k1.compare("P") == 0) {
            P_input = inputs[i].second;
        }
    }

    for (std::size_t i = 0; i < supported_pairs.size(); i++) {
        k1 = supported_pairs[i].first;
        k2 = supported_pairs[i].second;

        for (std::size_t j = 0; j < inputs.size(); j++) {
            if (inputs[j].first.compare(k1) == 0) {
                v1 = inputs[j].second;
            }
            else if (inputs[j].first.compare(k2) == 0) {
                v2 = inputs[j].second;
            }
        }

        try {
            double psi_w_new = HumidAir::HAPropsSI(
                "psi_w",
                k1, v1,
                k2, v2,
                "P", P_input);
        }
        catch (...) {
            double _nu = -_HUGE;
        }
    }
}

int main(int argc, const char* argv[]) {

    //CoolProp::set_debug_level(1);

    //double h = HumidAir::HAPropsSI("H", "T", 298.15, "P", 101325, "R", 0.5);
    //double T = HumidAir::HAPropsSI("T", "P", 101325, "H", h, "R", 1.0);
    //T = HumidAir::HAPropsSI("T", "H", h, "R", 1.0, "P", 101325);

    std::size_t num = 31;
    std::vector<double> T(num), R(num);
    for (std::size_t i = 0; i < num; i++) {
        T[i] = ((360.0 - 240.0) * i / double(num-1) + 240.0);
        R[i] = ((1.0 - 0.0) * i / double(num-1) + 0.0);
    }
    for (std::size_t i = 0; i < num; i++) {
        for (std::size_t j = 0; j < num; j++) {
            double Tv = T[i];
            double Rv = R[j];
            auto input_values = generate_values(Tv, Rv);
            calculate(input_values);
        }
    }
}


// # Humid air example from Sphinx
// from CoolProp.HumidAirProp import HAPropsSI
// h = HAPropsSI("H","T",298.15,"P",101325,"R",0.5); print(h)
// T = HAPropsSI("T","P",101325,"H",h,"R",1.0); print(T)
// T = HAPropsSI("T","H",h,"R",1.0,"P",101325); print(T)

// import sys
// sys.exit()

// # Verification script 
// import CoolProp.CoolProp as CP
// import numpy as np
// import itertools
// from multiprocessing import Pool

// def generate_values(TR,P=101325):
    // """ Starting with T,R as inputs, generate all other values """
    // T,R = TR
    // psi_w = CP.HAPropsSI("psi_w","T",T,"R",R,"P",P)
    // other_output_keys = ["T_wb","T_dp","Hda","Sda","Vda","Omega"]
    // outputs = {"psi_w":psi_w,"T":T,"P":P,"R":R}
    // for k in other_output_keys:
        // outputs[k] = CP.HAPropsSI(k,"T",T,"R",R,"P",P)
    // return outputs

// def get_supported_input_pairs():
    // """ Determine which input pairs are supported """
    // good_ones = []
    // inputs = generate_values((300, 0.5))
    // for k1, k2 in itertools.product(inputs.keys(), inputs.keys()):
        // if "P" in [k1,k2] or k1==k2:
            // continue
        // args = ("psi_w", k1, inputs[k1], k2, inputs[k2], "P", inputs["P"])
        // try:
            // psi_w_new = CP.HAPropsSI(*args)
            // good_ones.append((k1,k2))
        // except BaseException as BE:
            // pass
            // if "currently at least one of" in str(BE) or "cannot provide two inputs" in str(BE):
                // pass
            // else:
                // print(BE)
                // good_ones.append((k1,k2))
    // return good_ones

// def calculate(inputs):
    // """ For a given input, try all possible input pairs """
    // errors = []
    // supported_pairs = get_supported_input_pairs()
    // for k1, k2 in supported_pairs:
        // psi_w_input = inputs["psi_w"]
        // args = "psi_w",k1,inputs[k1],k2,inputs[k2],"P",inputs["P"]
        // try:
            // psi_w_new = CP.HAPropsSI(*args)
        // except BaseException as BE:
            // errors.append((str(BE),args, inputs))
    // return errors


// if __name__ == "__main__":
    // TR = itertools.product(np.linspace(240, 360, 31), np.linspace(0, 1, 31))
    // with Pool(processes=2) as pool:
        // input_values = pool.map(generate_values, TR)
        // errors = pool.map(calculate, input_values)
        // for err in itertools.chain.from_iterable(errors):
            // print(err)
