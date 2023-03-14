
#include <vector>
#include <string>
#include <chrono>

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
    std::vector<std::string> other_output_keys = {"T_wb", "T_dp", "Hda", "Sda", "Vda", "Omega"};
    std::vector<std::pair<std::string, double>> outputs;
    outputs.push_back(std::pair<std::string, double>("psi_w", psi_w));
    outputs.push_back(std::pair<std::string, double>("T", T));
    outputs.push_back(std::pair<std::string, double>("P", P));
    outputs.push_back(std::pair<std::string, double>("R", R));
    for (auto k : other_output_keys) {
        outputs.push_back(std::pair<std::string, double>(k, HumidAir::HAPropsSI(k, "T", T, "psi_w", psi_w, "P", P)));
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
                double psi_w_new = HumidAir::HAPropsSI("psi_w", k1, v1, k2, v2, "P", p);
                if (ValidNumber(psi_w_new)) {
                    good_ones.push_back(std::pair<std::string, std::string>(k1, k2));
                }
            } catch (std::exception& e) {
                std::cout << e.what();
            }
        }
    }
    return good_ones;
}

void calculate(std::vector<std::pair<std::string, double>> inputs, std::size_t& clc_count, std::size_t& err_count, std::size_t& acc_count,
               const std::vector<std::pair<std::string, std::string>>& supported_pairs) {
    //auto errors = []

    std::string k1, k2;
    double psi_w_input = -_HUGE, P_input = -_HUGE, v1 = -_HUGE, v2 = -_HUGE;

    for (const auto& kv : inputs) {
        if (kv.first == "psi_w") {
            psi_w_input = kv.second;
            break;
        }
    }
    for (const auto& kv : inputs) {
        if (kv.first == "P") {
            P_input = kv.second;
            break;
        }
    }

    for (std::size_t i = 0; i < supported_pairs.size(); i++) {
        k1 = supported_pairs[i].first;
        k2 = supported_pairs[i].second;

        for (std::size_t j = 0; j < inputs.size(); j++) {
            if (inputs[j].first.compare(k1) == 0) {
                v1 = inputs[j].second;
            } else if (inputs[j].first.compare(k2) == 0) {
                v2 = inputs[j].second;
            }
        }

        clc_count += 1;
        try {
            double psi_w_new = HumidAir::HAPropsSI("psi_w", k1, v1, k2, v2, "P", P_input);
            double delta = std::abs(psi_w_input - psi_w_new);
            if (delta > 1e-6) {
                acc_count += 1;
                HumidAir::HAPropsSI("psi_w", k1, v1, k2, v2, "P", P_input);
                std::cout << "deviation: " << delta << " @ HAPropsSI(\"psi_w\",\"" << k1 << "\"," << v1 << ",\"" << k2 << "\"," << v2
                          << ",\"P\",101325); error: " + CoolProp::get_global_param_string("errstring") << std::endl;

                //                std::cout << "\n-------------- Error --------------\n";
                //                std::cout << "delta = " << delta << "\n";
                //                std::cout << k1 << " = " << v1 << "\n";
                //                std::cout << k2 << " = " << v2 << "\n";
                //                std::cout << "P" << " = " << P_input << "\n";
            }
        } catch (std::exception& e) {
            err_count += 1;
            std::cout << e.what();
        }
    }
}

int main(int argc, const char* argv[]) {

    //CoolProp::set_debug_level(1);
    {
        //        for (auto R = 0.0; R < 1.0; R += 0.01){
        //            std::cout << R << " " << HumidAir::HAPropsSI("Hda", "T", 240, "R", R, "P", 101325) << "\n";
        //        }
        auto hh = HumidAir::HAPropsSI("psi_w", "R", 0.0333333, "Vda", 0.958997, "P", 101325);
        ;
        double h = HumidAir::HAPropsSI("S", "T", 240, "P", 101325, "R", 0);
        //        double T = HumidAir::HAPropsSI("W", "P", 101325, "S", h, "T", 240);
        //        T = HumidAir::HAPropsSI("T", "H", h, "R", 1.0, "P", 101325);
    }

    auto supported_pairs = get_supported_input_pairs();
    double time = 0, _time = 0;
    std::size_t err_count = 0, clc_count = 0, acc_count = 0;
    std::size_t _err_count = 0, _clc_count = 0, _acc_count = 0;
    std::size_t num = 31;
    std::vector<double> T(num), R(num);
    // Full range : -143.15 C to 350.0 C
    double T_lo = (-143.15 + 273.15) * 1.001;
    double T_hi = (350.00 + 273.15) * 0.999;
    // Full range : 0.0 to 1.0
    double R_lo = 0.0 * 1.001;
    double R_hi = 1.0 * 0.999;
    for (std::size_t i = 0; i < num; i++) {
        T[i] = ((T_hi - T_lo) * i / double(num - 1) + T_lo);
        R[i] = ((R_hi - R_lo) * i / double(num - 1) + R_lo);
    }
    for (std::size_t i = 0; i < num; i++) {
        _err_count = 0;
        _clc_count = 0;
        _acc_count = 0;
        _time = 0;
        auto tic = std::chrono::high_resolution_clock::now();
        double Tdb = T[i];
        for (std::size_t j = 0; j < num; j++) {

            double Rv = R[j];
            auto input_values = generate_values(Tdb, Rv);
            calculate(input_values, _clc_count, _err_count, _acc_count, supported_pairs);
        }
        auto toc = std::chrono::high_resolution_clock::now();
        _time = std::chrono::duration<double>(toc - tic).count();
        std::cout << "\n----- Errors for run " << i << " @ T_drybulb = " << T[i] << " K ----- \n";
        std::cout << "Exceptions: " << _err_count << " / " << _clc_count << " = " << _err_count * 100.0 / _clc_count << "% \n";
        std::cout << "Bad accuracy: " << _acc_count << " / " << _clc_count << " = " << _acc_count * 100.0 / _clc_count << "% \n";
        if (_clc_count != (R.size() * supported_pairs.size())) return 1;
        std::cout << "Time: " << _time << " s / " << _clc_count << " = " << _time / _clc_count * 1e3 << " ms per call \n";
        err_count += _err_count;
        clc_count += _clc_count;
        acc_count += _acc_count;
        time += _time;
    }
    std::cout << "\n----- Final Errors ----- \n";
    std::cout << "Exceptions: " << err_count << " / " << clc_count << " = " << err_count * 100.0 / clc_count << "% \n";
    std::cout << "Bad accuracy: " << acc_count << " / " << clc_count << " = " << acc_count * 100.0 / clc_count << "% \n";
    if (clc_count != (T.size() * R.size() * supported_pairs.size())) return 1;
    std::cout << "Time: " << time << " s / " << clc_count << " = " << time / clc_count * 1e3 << " ms per call \n";
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
