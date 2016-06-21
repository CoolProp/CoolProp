#include "CubicBackend.h"
#include "Solvers.h"
#include "Backends/Helmholtz/VLERoutines.h"

void CoolProp::AbstractCubicBackend::setup(bool generate_SatL_and_SatV){
    // Set the pure fluid flag
    is_pure_or_pseudopure = cubic->get_Tc().size() == 1;
    // Resize the vector
    resize(cubic->get_Tc().size());
	// Reset the residual Helmholtz energy class
	residual_helmholtz.reset(new CubicResidualHelmholtz(this));
	// If pure, set the mole fractions to be unity
	if (is_pure_or_pseudopure){
		mole_fractions = std::vector<CoolPropDbl>(1, 1.0);
        mole_fractions_double = std::vector<double>(1, 1.0);
	}
    else{
        mole_fractions.clear();
        mole_fractions_double.clear();
    }
	// Now set the reducing function for the mixture
    Reducing.reset(new ConstantReducingFunction(cubic->T_r, cubic->rho_r));

    // Top-level class can hold copies of the base saturation classes,
    // saturation classes cannot hold copies of the saturation classes
    if (generate_SatL_and_SatV)
    {
        bool SatLSatV = false;
        SatL.reset(this->get_copy(SatLSatV));
        SatL->specify_phase(iphase_liquid);
        linked_states.push_back(SatL);
        SatV.reset(this->get_copy(SatLSatV));
        SatV->specify_phase(iphase_gas);
        linked_states.push_back(SatV);
    }
}

void CoolProp::AbstractCubicBackend::get_linear_reducing_parameters(double &rhomolar_r, double &T_r){
    // In the case of models where the reducing temperature is not a function of composition (SRK, PR, etc.), 
    // we need to use an appropriate value for T_r and v_r, here we use a linear weighting
    T_r = 0; 
    double v_r = 0;
    const std::vector<double> Tc = cubic->get_Tc(), pc = cubic->get_pc();
    for (std::size_t i = 0; i < mole_fractions.size(); ++i){   
        T_r += mole_fractions[i]*Tc[i];
        // Curve fit from all the pure fluids in CoolProp (thanks to recommendation of A. Kazakov)
        double v_c_Lmol = 2.14107171795*(Tc[i]/pc[i]*1000)+0.00773144012514; // [L/mol]
        v_r += mole_fractions[i]*v_c_Lmol/1000.0;
    }
    rhomolar_r = 1/v_r;
}

void CoolProp::AbstractCubicBackend::get_critical_point_search_radii(double &R_delta, double &R_tau)
{
    // Get the starting values from the base class
    CoolProp::HelmholtzEOSMixtureBackend::get_critical_point_search_radii(R_delta, R_tau);

    // Now we scale them to get the appropriate search radii
    double Tr_GERGlike, rhor_GERGlike;
    get_linear_reducing_parameters(rhor_GERGlike, Tr_GERGlike);
    R_delta *= rhor_GERGlike/rhomolar_reducing()*5;
    R_tau *= T_reducing()/Tr_GERGlike*5;
}

bool CoolProp::AbstractCubicBackend::get_critical_is_terminated(double &delta, double &tau)
{
    // If the volume is less than the mixture covolume, stop.  The mixture covolume is the
    // smallest volume that is physically allowed for a cubic EOS
    double b = get_cubic()->bm_term(mole_fractions); // [m^3/mol]
    double v = 1/(delta*rhomolar_reducing()); //[m^3/mol]
    bool covolume_check = v < 1.1*b;
    
    return covolume_check;
}

void CoolProp::AbstractCubicBackend::get_critical_point_starting_values(double &delta0, double &tau0){
    
    // Get the starting values from the base class
    CoolProp::HelmholtzEOSMixtureBackend::get_critical_point_starting_values(delta0, tau0);

    // The starting tau and delta values can be thought of as being given with the 
    // GERG reducing values, or tau0 = Tr_GERG/T = 0.xxx and delta0 = rho/rhor_GERG = 0.xxx 
    // Then we need to multiply by 
    // 
    // Tr_cubic/Tr_GERGlike & rhor_GERGlike/rhor_cubic 
    // 
    // to get shifted values:
    //
    // delta0 = rho/rhor_GERG*(rhor_GERGlike/rhor_cubic)
    // tau0 = Tr_GERG/T*(Tr_cubic/Tr_GERGlike)
    // 
    double Tr_GERGlike, rhor_GERGlike;
    get_linear_reducing_parameters(rhor_GERGlike, Tr_GERGlike);
    delta0 *= rhor_GERGlike/rhomolar_reducing();
    tau0 *= T_reducing()/Tr_GERGlike;
}

CoolPropDbl CoolProp::AbstractCubicBackend::calc_alphar_deriv_nocache(const int nTau, const int nDelta, const std::vector<CoolPropDbl> & mole_fractions, const CoolPropDbl &tau, const CoolPropDbl &delta){
    bool cache_values = true;
    HelmholtzDerivatives derivs = residual_helmholtz->all(*this, mole_fractions, tau, delta, cache_values);
    switch (nTau){
        case 0:
        {
            switch (nDelta){
                case 0: return derivs.alphar;
                case 1: return derivs.dalphar_ddelta;
                case 2: return derivs.d2alphar_ddelta2;
                case 3: return derivs.d3alphar_ddelta3;
                case 4: return derivs.d4alphar_ddelta4;
                default: throw ValueError(format("nDelta (%d) is invalid",nDelta));
            }
            break;
        }
        case 1:
        {
            switch (nDelta){
                case 0: return derivs.dalphar_dtau;
                case 1: return derivs.d2alphar_ddelta_dtau;
                case 2: return derivs.d3alphar_ddelta2_dtau;
                case 3: return derivs.d4alphar_ddelta3_dtau;
                default: throw ValueError(format("nDelta (%d) is invalid",nDelta));
            }
            break;
        }
        case 2:
        {
            switch (nDelta){
                case 0: return derivs.d2alphar_dtau2;
                case 1: return derivs.d3alphar_ddelta_dtau2;
                case 2: return derivs.d4alphar_ddelta2_dtau2;
                default: throw ValueError(format("nDelta (%d) is invalid",nDelta));
            }
        }
        case 3:
        {
            switch (nDelta){
                case 0: return derivs.d3alphar_dtau3;
                case 1: return derivs.d4alphar_ddelta_dtau3;
                default: throw ValueError(format("nDelta (%d) is invalid",nDelta));
            }
        }
        case 4:
        {
            switch (nDelta){
                case 0: return derivs.d4alphar_dtau4;
                default: throw ValueError(format("nDelta (%d) is invalid",nDelta));
            }
        }
        default: throw ValueError(format("nTau (%d) is invalid",nTau));
    }
}

void CoolProp::AbstractCubicBackend::update(CoolProp::input_pairs input_pair, double value1, double value2){
    if (get_debug_level() > 10){std::cout << format("%s (%d): update called with (%d: (%s), %g, %g)",__FILE__,__LINE__, input_pair, get_input_pair_short_desc(input_pair).c_str(), value1, value2) << std::endl;}
    
    CoolPropDbl ld_value1 = value1, ld_value2 = value2;
    pre_update(input_pair, ld_value1, ld_value2);
    value1 = ld_value1; value2 = ld_value2;
    
    switch(input_pair)
    {
        case PT_INPUTS:
            _p = value1; _T = value2; _rhomolar = solver_rho_Tp(value2/*T*/, value1/*p*/); break;
        case QT_INPUTS:
            _Q = value1; _T = value2; saturation(input_pair); break;
        case PQ_INPUTS:
            _p = value1; _Q = value2; saturation(input_pair); break;
        case DmolarT_INPUTS:
        case SmolarT_INPUTS:
        case DmolarP_INPUTS:
        case DmolarHmolar_INPUTS:
        case DmolarSmolar_INPUTS:
        case DmolarUmolar_INPUTS:
        case HmolarP_INPUTS:
        case PSmolar_INPUTS:
        case PUmolar_INPUTS:
        case HmolarSmolar_INPUTS:

        case QSmolar_INPUTS:
        case HmolarQ_INPUTS:
        case DmolarQ_INPUTS:
            HelmholtzEOSMixtureBackend::update(input_pair, value1, value2); break;
        default:
            throw ValueError(format("This pair of inputs [%s] is not yet supported", get_input_pair_short_desc(input_pair).c_str()));
    }
    
    post_update();
}

void CoolProp::AbstractCubicBackend::rho_Tp_cubic(CoolPropDbl T, CoolPropDbl p, int &Nsolns, double &rho0, double &rho1, double &rho2){
    AbstractCubic *cubic = get_cubic().get();
    double R = cubic->get_R_u();
    double Delta_1 = cubic->get_Delta_1();
    double Delta_2 = cubic->get_Delta_2();
    double am = cubic->am_term(cubic->T_r/T, mole_fractions_double, 0);
    double bm = cubic->bm_term(mole_fractions);
    
    // ** SYMPY CODE ***
    // R,T,v,b,a,Delta_1,Delta_2,p,Z,rho = symbols('R,T,v,b,a,Delta_1,Delta_2,p,Z,rho')
    // eqn = (R*T/(v-b)-a/(v+Delta_1*b)/(v+Delta_2*b)-p).subs(v,1/rho)
    // eqn2 = eqn*(-b + 1/rho)*(Delta_1*b + 1/rho)*(Delta_2*b + 1/rho)
    // display(simplify(factor(expand(eqn2),rho)))
    //
    // yields:
    // (b*rho**3*(Delta_1*Delta_2*R*T*b + Delta_1*Delta_2*b**2*p + a) - p + rho**2*(-Delta_1*Delta_2*b**2*p + Delta_1*R*T*b + Delta_1*b**2*p + Delta_2*R*T*b + Delta_2*b**2*p - a) - rho*(Delta_1*b*p + Delta_2*b*p - R*T - b*p))/rho**3
    
    double crho0 = -p;
    double crho1 = -1*((Delta_1+Delta_2-1)*bm*p - R*T);
    double crho2 = -Delta_1*Delta_2*bm*bm*p + (Delta_1+Delta_2)*(R*T*bm + bm*bm*p) - am;
    double crho3 = bm*(Delta_1*Delta_2*(R*T*bm + bm*bm*p) + am);
    solve_cubic(crho3, crho2, crho1, crho0, Nsolns, rho0, rho1, rho2);
    sort3(rho0, rho1, rho2);
    return;
    
//    double A = cubic->am_term(cubic->T_r/T, mole_fractions_double, 0)*p/(POW2(R*T));
//    double B = cubic->bm_term(mole_fractions)*p/(R*T);
//    double Z0=0, Z1=0, Z2=0;
//    solve_cubic(1,
//                B*(Delta_1+Delta_2-1)-1,
//                A + B*B*(Delta_1*Delta_2-Delta_1-Delta_2) - B*(Delta_1+Delta_2),
//                -A*B-Delta_1*Delta_2*(POW2(B)+POW3(B)),
//                Nsolns, Z0, Z1, Z2);
//    if (Nsolns == 1){ rho0 = p/(Z0*R*T); }
//    else if (Nsolns == 3){
//        rho0 = p/(Z0*R*T);
//        rho1 = p/(Z1*R*T);
//        rho2 = p/(Z2*R*T);
//        sort3(rho0, rho1, rho2);
//    }
}

class SaturationResidual : public CoolProp::FuncWrapper1D{
public:
    CoolProp::AbstractCubicBackend *ACB;
    CoolProp::input_pairs inputs;
    double imposed_variable;
    double deltaL, deltaV;

    SaturationResidual(){};
    SaturationResidual(CoolProp::AbstractCubicBackend *ACB, CoolProp::input_pairs inputs, double imposed_variable) : ACB(ACB), inputs(inputs), imposed_variable(imposed_variable) {};
    
    double call(double value){
        int Nsolns = 0;
        double rho0 = -1, rho1 = -1, rho2 = -1, T, p;
        
        if (inputs == CoolProp::PQ_INPUTS){
            T = value; p = imposed_variable;
        }
        else if (inputs == CoolProp::QT_INPUTS){
            p = value; T = imposed_variable;
        }
        else{
            throw CoolProp::ValueError("Cannot have something other than PQ_INPUTS or QT_INPUTS here");
        }
        
        // Calculate the liquid and vapor densities
        ACB->rho_Tp_cubic(T, p, Nsolns, rho0, rho1, rho2);
        
        // -----------------------------------------------------
        // Calculate the difference in Gibbs between the phases
        // -----------------------------------------------------
        AbstractCubic *cubic = ACB->get_cubic().get();
        double rho_r = cubic->rho_r, T_r = cubic->T_r;
        double tau = T_r/T;
        // There are three density solutions, we know the highest is the liquid, the lowest is the vapor
        deltaL = rho2/rho_r; deltaV = rho0/rho_r;
        // From alpha0; all terms that are only a function of temperature drop out since TL=TV
        double DELTAgibbs = log(deltaV) - log(deltaL);
        // From alphar;
        DELTAgibbs += (cubic->alphar(tau, deltaV, ACB->get_mole_fractions_doubleref(), 0, 0)
                       -cubic->alphar(tau, deltaL, ACB->get_mole_fractions_doubleref(), 0, 0));
        // From delta*dalphar_dDelta
        DELTAgibbs += (deltaV*cubic->alphar(tau, deltaV, ACB->get_mole_fractions_doubleref(), 0, 1)
                       -deltaL*cubic->alphar(tau, deltaL, ACB->get_mole_fractions_doubleref(), 0, 1));
        return DELTAgibbs;
    };
};

std::vector<double> CoolProp::AbstractCubicBackend::spinodal_densities(){
    //// SPINODAL
    AbstractCubic *cubic = get_cubic().get();
    double tau = cubic->T_r / _T;
    std::vector<double> x(1, 1);
    double a = cubic->am_term(tau, x, 0);
    double b = cubic->bm_term(x);
    double R = cubic->get_R_u();
    double Delta_1 = cubic->get_Delta_1();
    double Delta_2 = cubic->get_Delta_2();

    double crho4 = -powInt(Delta_1*Delta_2, 2)*R*_T*powInt(b, 4) + a*powInt(b, 3)*(Delta_1 + Delta_2);
    double crho3 = -2 * ((Delta_1*Delta_1*Delta_2 + Delta_1*Delta_2*Delta_2)*R*_T*powInt(b, 3) + a*powInt(b, 2)*(Delta_1 + Delta_2 - 1));
    double crho2 = ((-Delta_1*Delta_1 - Delta_2*Delta_2 - 4 * Delta_1*Delta_2)*R*_T*powInt(b, 2) + a*b*(Delta_1 + Delta_2 - 4));
    double crho1 = -2 * (Delta_1 + Delta_2)*R*_T*b + 2 * a;
    double crho0 = -R*_T;
    double rho0, rho1, rho2, rho3;
    int Nsoln;
    solve_quartic(crho4, crho3, crho2, crho1, crho0, Nsoln, rho0, rho1, rho2, rho3);
    std::vector<double> roots;
    if (rho0 > 0 && 1 / rho0 > b) { roots.push_back(rho0); }
    if (rho1 > 0 && 1 / rho1 > b) { roots.push_back(rho1); }
    if (rho2 > 0 && 1 / rho2 > b) { roots.push_back(rho2); }
    if (rho3 > 0 && 1 / rho3 > b) { roots.push_back(rho3); }
    return roots;
}

void CoolProp::AbstractCubicBackend::saturation(CoolProp::input_pairs inputs){
    AbstractCubic *cubic = get_cubic().get();
    double Tc = cubic->get_Tc()[0], pc = cubic->get_pc()[0], acentric = cubic->get_acentric()[0];
    double rhoL=-1, rhoV=-1;
    if (inputs == PQ_INPUTS){
        if (is_pure_or_pseudopure){
            // Estimate temperature from the acentric factor relationship
            double theta = -log10(_p/pc)*(1/0.7-1)/(acentric+1);
            double Ts_est = Tc/(theta+1);
            SaturationResidual resid(this, inputs, _p);
            static std::string errstr;
            double Ts = CoolProp::Secant(resid, Ts_est, -0.1, 1e-10, 100);
            _T = Ts;
            rhoL = resid.deltaL*cubic->T_r;
            rhoV = resid.deltaV*cubic->T_r;
        }
        else{
            HelmholtzEOSMixtureBackend::update(PQ_INPUTS, _p, _Q);
            return;
        }
    }
    else if (inputs == QT_INPUTS){
        if (is_pure_or_pseudopure){
            SaturationResidual resid(this, inputs, _T);
            static std::string errstr;
            std::vector<double> roots = spinodal_densities();

            // Estimate pressure from the acentric factor relationship
            double neg_log10_pr = (acentric + 1) / (1 / 0.7 - 1)*(Tc / _T - 1);
            double ps_est = pc*pow(10.0, -neg_log10_pr);

            double ps;
            if (roots.size() == 2){
                double p0 = calc_pressure_nocache(_T, roots[0]);
                double p1 = calc_pressure_nocache(_T, roots[1]);
                if (p1 < p0){ std::swap(p0, p1); }
                //ps = CoolProp::BoundedSecant(resid, p0, p1, pc, -0.01*ps_est, 1e-5, 100);
                if (p0 > 0 && p1 < pc){
                    ps = CoolProp::Brent(resid, p0*1.0001, p1*0.9999, DBL_EPSILON, 1e-10, 100);
                }
                else{
                    ps = CoolProp::BoundedSecant(resid, ps_est, 1e-10, pc, -0.01*ps_est, 1e-5, 100);
                }
            }
            else{
                ps = CoolProp::BoundedSecant(resid, ps_est, 1e-10, pc, -0.01*ps_est, 1e-5, 100);
            }
            
            _p = ps;
            rhoL = resid.deltaL*cubic->T_r;
            rhoV = resid.deltaV*cubic->T_r;
            this->SatL->update(DmolarT_INPUTS, rhoL, _T);
            this->SatV->update(DmolarT_INPUTS, rhoV, _T);
        }
        else{
            HelmholtzEOSMixtureBackend::update(QT_INPUTS, _Q, _T);
            return;
        }
    }
    _rhomolar = 1/(_Q/rhoV+(1-_Q)/rhoL);
    _phase = iphase_twophase;
}

CoolPropDbl CoolProp::AbstractCubicBackend::solver_rho_Tp(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rho_guess){
    int Nsoln = 0;
    double rho0=0, rho1=0, rho2=0, rho = -1;
    rho_Tp_cubic(T, p, Nsoln, rho0, rho1, rho2); // Densities are sorted in increasing order
    if (Nsoln == 1){
        rho = rho0;
    }
    else if (Nsoln == 3){
        if (imposed_phase_index != iphase_not_imposed){
            // Use imposed phase to select root
            if (imposed_phase_index == iphase_gas || imposed_phase_index == iphase_supercritical_gas){
                rho = rho0;
            }
            else if (imposed_phase_index == iphase_liquid || imposed_phase_index == iphase_supercritical_liquid){
                rho = rho2;
            }
            else{
                throw ValueError("Specified phase is invalid");
            }
        }
        else{
            throw ValueError("Cubic has three roots, but phase not imposed and guess density not provided");
        }
    }
    else{
        throw ValueError("Obtained neither 1 nor three roots");
    }
    if (is_pure_or_pseudopure){
        // Set some variables at the end
        this->recalculate_singlephase_phase();
    }
    else{
        _phase = iphase_gas; // TODO: fix this
    }
    _Q = -1;
    return rho;
}

void CoolProp::AbstractCubicBackend::set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string &parameter, const double value){
    if (parameter == "kij" || parameter == "k_ij"){
        get_cubic()->set_kij(i,j,value);
        if (this->SatL.get() != NULL){ this->SatL->set_binary_interaction_double(i,j,"kij",value); }
        if (this->SatV.get() != NULL){ this->SatV->set_binary_interaction_double(i,j,"kij",value); }
    }
    else{
        throw ValueError(format("I don't know what to do with parameter [%s]", parameter.c_str()));
    }
};
double CoolProp::AbstractCubicBackend::get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string &parameter){
    if (parameter == "kij" || parameter == "k_ij"){
        return get_cubic()->get_kij(i,j);
    }
    else{
        throw ValueError(format("I don't know what to do with parameter [%s]", parameter.c_str()));
    }
};

void CoolProp::AbstractCubicBackend::copy_k(AbstractCubicBackend *donor){
    std::size_t N = mole_fractions.size();
    for (std::size_t i = 0; i < N; ++i){
        for (std::size_t j = i+1; j < N; ++j){
            double val = donor->get_binary_interaction_double(i, j, "kij");
            this->set_binary_interaction_double(i, j, "kij", val);
        }
    }
}

void CoolProp::AbstractCubicBackend::set_C_MC(double c1, double c2, double c3)
{
    get_cubic()->set_C_MC(c1, c2, c3);

    if (SatL.get() != NULL){ 
        /// Cast to this type so that we can access its internal methods
        CoolProp::AbstractCubicBackend* _SatL = static_cast<CoolProp::AbstractCubicBackend*>(this->SatL.get());
        _SatL->get_cubic()->set_C_MC(c1, c2, c3);
    }
    if (SatV.get() != NULL){ 
        /// Cast to this type so that we can access its internal methods
        CoolProp::AbstractCubicBackend* _SatV = static_cast<CoolProp::AbstractCubicBackend*>(this->SatV.get());
        _SatV->get_cubic()->set_C_MC(c1, c2, c3);
    }
}
