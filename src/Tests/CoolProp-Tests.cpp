

#include "AbstractState.h"
#include "DataStructures.h"

// ############################################
//                      TESTS
// ############################################

#if defined(ENABLE_CATCH)

#include "catch.hpp"

namespace ViscosityValidation{

// A structure to hold the values for one validation call for viscosity
struct vel
{
public:
    std::string in1, in2, out, fluid;
    double v1, v2, tol, expected;
    vel(std::string fluid, std::string in1, double v1, std::string in2, double v2, std::string out, double expected, double tol)
    {
        this->in1 = in1; this->in2 = in2; this->fluid = fluid;
        this->v1 = v1; this->v2 = v2; this->expected = expected; this->tol = tol;
    };
};

vel viscosity_validation_data[] = {
// From Vogel, JPCRD, 1998
vel("Propane", "T", 90, "Dmolar", 16.52e3, "V", 7388e-6, 1e-3),
vel("Propane", "T", 150, "Dmolar", 15.14e3, "V", 656.9e-6, 5e-3),
vel("Propane", "T", 600, "Dmolar", 10.03e3, "V", 73.92e-6, 5e-3),
vel("Propane", "T", 280, "Dmolar", 11.78e3, "V", 117.4e-6,1e-3),

// From Michailidou, JPCRD, 2013
vel("Hexane", "T", 250, "Dmass", 1e-14, "V", 5.2584e-6, 1e-3),
vel("Hexane", "T", 400, "Dmass", 1e-14, "V", 8.4149e-6, 1e-3),
vel("Hexane", "T", 550, "Dmass", 1e-14, "V", 11.442e-6, 1e-3),
vel("Hexane", "T", 250, "Dmass", 700, "V", 528.2e-6, 1e-3),
vel("Hexane", "T", 400, "Dmass", 600, "V", 177.62e-6, 1e-3),
vel("Hexane", "T", 550, "Dmass", 500, "V", 95.002e-6, 1e-3),

// From Vesovic, JPCRD, 1990
vel("CO2", "T", 220, "Dmass", 2.440, "V", 11.06e-6, 1e-3),
vel("CO2", "T", 300, "Dmass", 1.773, "V", 15.02e-6, 1e-3),
vel("CO2", "T", 800, "Dmass", 0.662, "V", 35.09e-6, 1e-3),
vel("CO2", "T", 304, "Dmass", 254.3205, "V", 20.89e-6, 1e-3),
vel("CO2", "T", 220, "Dmass", 1194.86, "V", 274.22e-6, 1e-3),
vel("CO2", "T", 300, "Dmass", 1029.27, "V", 133.15e-6, 1e-3),
vel("CO2", "T", 800, "Dmass", 407.828, "V", 48.62e-6, 1e-3),

// From Lemmon and Jacobsen, JPCRD, 2004
vel("Nitrogen", "T", 100, "Dmolar", 1e-14, "V", 6.90349e-6, 1e-3),
vel("Nitrogen", "T", 300, "Dmolar", 1e-14, "V", 17.8771e-6, 1e-3),
vel("Nitrogen", "T", 100, "Dmolar", 25000, "V", 79.7418e-6, 1e-3),
vel("Nitrogen", "T", 200, "Dmolar", 10000, "V", 21.0810e-6, 1e-3),
vel("Nitrogen", "T", 300, "Dmolar", 5000, "V", 20.7430e-6, 1e-3),
vel("Nitrogen", "T", 126.195, "Dmolar", 11180, "V", 18.2978e-6, 1e-3),
vel("Argon", "T", 100, "Dmolar", 1e-14, "V", 8.18940e-6, 1e-3),
vel("Argon", "T", 300, "Dmolar", 1e-14, "V", 22.7241e-6, 1e-3),
vel("Argon", "T", 100, "Dmolar", 33000, "V", 184.232e-6, 1e-3),
vel("Argon", "T", 200, "Dmolar", 10000, "V", 25.5662e-6, 1e-3),
vel("Argon", "T", 300, "Dmolar", 5000, "V", 26.3706e-6, 1e-3),
vel("Argon", "T", 150.69, "Dmolar", 13400, "V", 27.6101e-6, 1e-3),
vel("Oxygen", "T", 100, "Dmolar", 1e-14, "V", 7.70243e-6, 1e-3),
vel("Oxygen", "T", 300, "Dmolar", 1e-14, "V", 20.6307e-6, 1e-3),
vel("Oxygen", "T", 100, "Dmolar", 35000, "V", 172.136e-6, 1e-3),
vel("Oxygen", "T", 200, "Dmolar", 10000, "V", 22.4445e-6, 1e-3),
vel("Oxygen", "T", 300, "Dmolar", 5000, "V", 23.7577e-6, 1e-3),
vel("Oxygen", "T", 150.69, "Dmolar", 13600, "V", 24.7898e-6, 1e-3),
vel("Air", "T", 100, "Dmolar", 1e-14, "V", 7.09559e-6, 1e-3),
vel("Air", "T", 300, "Dmolar", 1e-14, "V", 18.5230e-6, 1e-3),
vel("Air", "T", 100, "Dmolar", 28000, "V", 107.923e-6, 1e-3),
vel("Air", "T", 200, "Dmolar", 10000, "V", 21.1392e-6, 1e-3),
vel("Air", "T", 300, "Dmolar", 5000, "V", 21.3241e-6, 1e-3),
vel("Air", "T", 132.64, "Dmolar", 10400, "V", 17.7623e-6, 1e-3),

// From Fenhour, JPCRD, 1995
vel("Ammonia", "T", 200, "Dmolar", 3.9, "V", 6.95e-6, 1e-3),
vel("Ammonia", "T", 200, "Dmolar", 42754.4, "V", 507.28e-6, 1e-3),
vel("Ammonia", "T", 398, "Dmolar", 7044.7, "V", 17.67e-6, 1e-3),
vel("Ammonia", "T", 398, "Dmolar", 21066.7, "V", 43.95e-6, 1e-3),

// Huber, Energy & Fuels, 2004
vel("n-Dodecane", "T", 300, "Dmolar", 4411.5, "V", 1484.8e-6, 1e-3),
vel("n-Dodecane", "T", 500, "Dmolar", 3444.7, "V", 183.76e-6, 1e-3),

// Huber, I&ECR, 2006
vel("R125", "T", 300, "Dmolar", 10596.9998, "V", 177.37e-6, 1e-3),
vel("R125", "T", 400, "Dmolar", 30.631, "V", 17.070e-6, 1e-3),

// Huber, FPE, 2004
vel("n-Octane", "T", 300, "Dmolar", 6177.2, "V", 553.60e-6, 1e-3),
vel("n-Nonane", "T", 300, "Dmolar", 5619.1, "V", 709.53e-6, 1e-3),
vel("n-Decane", "T", 300, "Dmolar", 5150.4, "V", 926.44e-6, 1e-3),

};

class ViscosityValidationFixture
{
protected:
    long double actual, x1, x2;
    CoolProp::AbstractState *pState;
    int pair;
public:
    ViscosityValidationFixture(){ pState = NULL; }
    ~ViscosityValidationFixture(){ delete pState; }
    void set_backend(std::string backend, std::string fluid_name){
        pState = CoolProp::AbstractState::factory(backend, fluid_name);
    }
    void set_pair(std::string &in1, double v1, std::string &in2, double v2){ 
        double o1, o2;
        long iin1 = CoolProp::get_parameter_index(in1);
        long iin2 = CoolProp::get_parameter_index(in2);
        long pair = CoolProp::generate_update_pair(iin1, v1, iin2, v2, o1, o2);
        pState->update(pair, o1, o2);
    }
    void get_value()
    {
        actual = pState->viscosity();
    }
};

TEST_CASE_METHOD(ViscosityValidationFixture, "Compare viscosities against published data", "[viscosity]")
{
    int inputsN = sizeof(viscosity_validation_data)/sizeof(viscosity_validation_data[0]);
    for (int i = 0; i < inputsN; ++i)
    {
        vel el = viscosity_validation_data[i];
        CHECK_NOTHROW(set_backend("REFPROP", el.fluid));

        CAPTURE(el.fluid);
        CAPTURE(el.in1);
        CAPTURE(el.v1);
        CAPTURE(el.in2);
        CAPTURE(el.v2);
        CHECK_NOTHROW(set_pair(el.in1, el.v1, el.in2, el.v2));
        get_value();
        CAPTURE(el.expected);
        CAPTURE(actual);
        CHECK(fabs(actual/el.expected-1) < el.tol);
    }
}

}; /* namespace ViscosityValidation */

static int inputs[] = {
    CoolProp::DmolarT_INPUTS,
    CoolProp::SmolarT_INPUTS,
    CoolProp::HmolarT_INPUTS, 
    CoolProp::TUmolar_INPUTS,

    CoolProp::DmolarP_INPUTS, 
    CoolProp::DmolarHmolar_INPUTS, 
    CoolProp::DmolarSmolar_INPUTS, 
    CoolProp::DmolarUmolar_INPUTS,
        
    /*
    CoolProp::HmolarP_INPUTS,
    CoolProp::PSmolar_INPUTS,
    CoolProp::PUmolar_INPUTS, 
    */

    /*
    CoolProp::HmolarSmolar_INPUTS, 
    CoolProp::HmolarUmolar_INPUTS, 
    CoolProp::SmolarUmolar_INPUTS 
    */
};

class ConsistencyFixture
{
protected:
    long double hmolar, pmolar, smolar, umolar, rhomolar, T, p, x1, x2;
    CoolProp::AbstractState *pState;
    int pair;
public:
    ConsistencyFixture(){
        pState = NULL;
    }
    ~ConsistencyFixture(){
        delete pState;
    }
    void set_backend(std::string backend, std::string fluid_name){
        pState = CoolProp::AbstractState::factory(backend, fluid_name);
    }
    void set_pair(int pair){ 
        this->pair = pair;
    }
    void set_TP(long double T, long double p)
    {
        this->T = T; this->p = p;
        CoolProp::AbstractState &State = *pState;

        // Start with T,P as inputs, cycle through all the other pairs that are supported
        State.update(CoolProp::PT_INPUTS, p, T);
            
        // Set the other state variables
        rhomolar = State.rhomolar(); hmolar = State.hmolar(); smolar = State.smolar(); umolar = State.umolar();
    }
    void get_variables()
    {
        CoolProp::AbstractState &State = *pState;
            
        switch (pair)
        {
        /// In this group, T is one of the known inputs, iterate for the other one (easy)
        case CoolProp::HmolarT_INPUTS:
            x1 = hmolar; x2 = T;  break;
        case CoolProp::SmolarT_INPUTS:
            x1 = smolar; x2 = T; break;
        case CoolProp::TUmolar_INPUTS:
            x1 = T; x2 = umolar; break;
        case CoolProp::DmolarT_INPUTS:
            x1 = rhomolar; x2 = T; break;

        /// In this group, D is one of the known inputs, iterate for the other one (a little bit harder)
        case CoolProp::DmolarHmolar_INPUTS:
            x1 = rhomolar; x2 = hmolar; break;
        case CoolProp::DmolarSmolar_INPUTS:
            x1 = rhomolar; x2 = smolar; break;
        case CoolProp::DmolarUmolar_INPUTS:
            x1 = rhomolar; x2 = umolar; break;
        case CoolProp::DmolarP_INPUTS:
            x1 = rhomolar; x2 = p; break;

        /// In this group, p is one of the known inputs (a little less easy)
        case CoolProp::HmolarP_INPUTS:
            x1 = hmolar; x2 = p; break;
        case CoolProp::PSmolar_INPUTS:
            x1 = p; x2 = smolar; break;
        case CoolProp::PUmolar_INPUTS:
            x1 = p; x2 = umolar; break;

        case CoolProp::HmolarSmolar_INPUTS:
            x1 = hmolar; x2 = smolar; break;
        case CoolProp::SmolarUmolar_INPUTS:
            x1 = smolar; x2 = umolar; break;
        }
    }
    void single_phase_consistency_check()
    {
        CoolProp::AbstractState &State = *pState;
        State.update(pair, x1, x2);

        // Make sure we end up back at the same temperature and pressure we started out with
        if(fabs(T-State.T()) > 1e-2) throw CoolProp::ValueError(format("Error on T [%g K] is greater than 1e-2",fabs(State.T()-T)));
        if(fabs(p-State.p())/p*100 > 1e-2)  throw CoolProp::ValueError(format("Error on p [%g %%] is greater than 1e-2 %%",fabs(p-State.p())/p ));
    }
};

TEST_CASE_METHOD(ConsistencyFixture, "Test all input pairs for CO2 using all valid backends", "[]")
{
    CHECK_NOTHROW(set_backend("HEOS", "CO2"));
        
    int inputsN = sizeof(inputs)/sizeof(inputs[0]);
    for (double p = 600000; p < 800000000.0; p *= 5)
    {
        for (double T = 220; T < pState->Tmax(); T += 5)
        {
            CHECK_NOTHROW(set_TP(T, p));

            for (int i = 0; i < inputsN; ++i)
            {
                int pair = inputs[i];
                std::string pair_desc = CoolProp::get_input_pair_short_desc(pair);
                set_pair(pair);
                CAPTURE(pair_desc);
                CAPTURE(T);
                CAPTURE(p);
                get_variables();
                CAPTURE(x1);
                CAPTURE(x2);
                CHECK_NOTHROW(single_phase_consistency_check());
            }
        }
    }
}

#endif
