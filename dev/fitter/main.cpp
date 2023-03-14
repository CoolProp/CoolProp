#include <iostream>
#include "Eigen/Dense"
#include "time.h"
#include "Helmholtz.h"
#include "CoolProp.h"

class EOSFitter;
#include "Fitter.h"
#include "DataTypes.h"

int main() {
    double n[] = {0.0,           0.5586817e-3, 0.4982230e0,   0.2458698e-0,  0.8570145e-3,  0.4788584e-3,  -0.1800808e-1, 0.2671641e0,
                  -0.4781652e1,  0.1423987e1,  0.3324062e0,   -0.7485907e-2, 0.1017263e-3,  -0.5184567e+0, -0.8692288e-1, 0.2057144e+0,
                  -0.5000457e-2, 0.4603262e-3, -0.3497836e-2, 0.6995038e-2,  -0.1452184e-1, -0.1285458e-3};
    double d[] = {0, 2, 1, 3, 6, 6, 1, 1, 2, 5, 2, 2, 4, 1, 4, 1, 2, 4, 1, 5, 3, 10};
    double t[] = {0.0, -1.0 / 2.0, 0.0, 0.0, 0.0, 3.0 / 2.0, 3.0 / 2.0, 2.0,  2.0,  1.0,  3.0,
                  5.0, 1.0,        5.0, 5.0, 6.0, 10.0,      10.0,      10.0, 18.0, 22.0, 50.0};
    double c[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0};

    std::vector<double> nv(n, n + sizeof(n) / sizeof(double));

    double mm = Props1SI("R134a", "molemass");
    double rhoL, rhoV;
    bool supercritical_T;

    double Tr = Props1SI("R134a", "Treduce");
    EOSFitter* pEOS = new EOSFitterFixedForm(Props1SI("R134a", "Treduce"), Props1SI("R134a", "rhoreduce") / mm * 1000, 8.314471);
    EOSFitter& EOS = *pEOS;

    // ----------------------------
    // Generate "experimental" data
    // ----------------------------
    for (double T = 250; T < 500; T += 10) {
        if (T < Tr) {
            rhoL = PropsSI("D", "T", T, "Q", 0, "R134a");
            rhoV = PropsSI("D", "T", T, "Q", 1, "R134a");
            supercritical_T = false;
        } else {
            rhoL = -1;
            rhoV = -1;
            supercritical_T = true;
        }

        for (double rho = 1e-10; rho < 1200; rho *= 1.5) {
            if (!supercritical_T && (rho < rhoL && rho > rhoV)) {
                continue;
            }
            double p = PropsSI("P", "T", T, "D", rho, "R134a");
            double rhobar = rho / mm * 1000;
            double cp = PropsSI("C", "T", T, "D", rho, "R134a");  // [J/kg/K]; convert to J/mol/K by *mm/1000
            double variance = 1;                                  // TODO; change this
            EOS.linear_data_points.push_back(new PressureDataPoint(pEOS, T, rho / mm * 1000, p, variance));
            EOS.nonlinear_data_points.push_back(new SpecificHeatCPDataPoint(pEOS, T, rho / mm * 1000, cp * mm / 1000, variance * 100));
        }
    }

    // Setup the EOS
    EOS.alphar = phir_power(n, d, t, c, 1, 21, 22);

    static const double a0[] = {
      0.0,        //[0]
      -1.019535,  //[1]
      9.047135,   //[2]
      -1.629789,  //[3]
      -9.723916,  //[4]
      -3.927170   //[5]
    };
    static const double t0[] = {
      0.0,         //[0]
      0.0,         //[1]
      0.0,         //[2]
      0.0,         //[3]
      -1.0 / 2.0,  //[4]
      -3.0 / 4.0   //[5]
    };

    // phi0=log(delta)+a0[1]+a0[2]*tau+a0[3]*log(tau)+a0[4]*pow(tau,-1.0/2.0)+a0[5]*pow(tau,-3.0/4.0);
    EOS.alpha0.push_back(new phi0_lead(a0[1], a0[2]));
    EOS.alpha0.push_back(new phi0_logtau(a0[3]));
    EOS.alpha0.push_back(new phi0_power(a0, t0, 4, 5, 6));

    /*for (unsigned int i = 0; i < EOS.nonlinear_data_points.size();i++)
	{
		std::cout << EOS.nonlinear_data_points[i]->residual(nv) << std::endl;
	}*/

    // Set the coefficients in the preliminary EOS
    EOS.set_n(nv);
    std::cout << format("before fit x2 %g\n", EOS.sum_squares(nv, false));
    // Solve for n without nonlinear terms to get an approximate solution
    EOS.solve_for_n(nv, false);
    std::cout << format("solved for n x2 %g\n", EOS.sum_squares(nv, false));
    EOS.set_n(nv);
    std::cout << format("applied n x2 %g\n", EOS.sum_squares(nv, true));

    for (int iter = 0; iter < 5; iter++) {
        EOS.set_n(nv);

        // Turn on the nonlinear terms and try again
        EOS.solve_for_n(nv, true);

        std::cout << nv[1] << " " << nv[2] << std::endl;

        std::cout << format("iter: %d x2 %g\n", iter, EOS.sum_squares(nv, true));
    }

    double rr = 0;
}