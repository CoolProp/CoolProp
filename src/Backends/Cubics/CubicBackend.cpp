#include "CubicBackend.h"
#include "Solvers.h"
#include "Backends/Helmholtz/VLERoutines.h"

void CoolProp::AbstractCubicBackend::setup(bool generate_SatL_and_SatV){
    // Set the pure fluid flag
    is_pure_or_pseudopure = cubic->get_Tc().size() == 1;
    // Resize the vector
    resize(cubic->get_Tc().size());
    N = cubic->get_Tc().size();
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
CoolPropDbl CoolProp::AbstractCubicBackend::calc_pressure_nocache(CoolPropDbl T, CoolPropDbl rhomolar){
    AbstractCubic *cubic = get_cubic().get();
    double tau = cubic->T_r / T;
    double delta = rhomolar / cubic->rho_r;
    return _rhomolar*gas_constant()*_T*(1+delta*cubic->alphar(tau, delta, this->get_mole_fractions_doubleref(), 0, 1));
}
void CoolProp::AbstractCubicBackend::update_DmolarT()
{
    // Only works for now when phase is specified
    if (this->imposed_phase_index != iphase_not_imposed){
        _p = calc_pressure_nocache(_T, _rhomolar);
        _Q = -1;
        _phase = imposed_phase_index;
    }
    else{
        // Pass call to parent class
        HelmholtzEOSMixtureBackend::update(DmolarT_INPUTS, this->_rhomolar, this->_T);
    }
};

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
            _rhomolar = value1; _T = value2; update_DmolarT(); break;
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
    double am = cubic->am_term(cubic->T_r/T, mole_fractions_double, 0);
    double bm = cubic->bm_term(mole_fractions);
    double cm = cubic->cm_term();

    // Introducing new variables to simplify the equation:
    double d1 = cm - bm;
    double d2 = cm + cubic->get_Delta_1()*bm;
    double d3 = cm + cubic->get_Delta_2()*bm;
    
    // Cubic coefficients:
    double crho0 = -p;
    double crho1 = R*T - p*(d1+d2+d3);
    double crho2 = R*T*(d2 + d3) - p*(d1*(d2 + d3) + d2*d3) - am;
    double crho3 = R*T*d2*d3 - p*d1*d2*d3 - d1*am;

    // Solving the cubic:
    solve_cubic(crho3, crho2, crho1, crho0, Nsolns, rho0, rho1, rho2);
    sort3(rho0, rho1, rho2);
    return;
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
            this->SatL->update(DmolarT_INPUTS, rhoL, _T);
            this->SatV->update(DmolarT_INPUTS, rhoV, _T);
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
            // ** Spinodal densities is disabled for now because it is VERY slow :(
            // std::vector<double> roots = spinodal_densities();
            std::vector<double> roots;

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
CoolPropDbl CoolProp::AbstractCubicBackend::solver_rho_Tp_global(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rhomolar_max)
{
    _rhomolar = solver_rho_Tp(T, p, 40000);
    return static_cast<double>(_rhomolar);
    
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
                if (rho0 > 0){
                    rho = rho0;
                }
                else if (rho1 > 0){
                    rho = rho1;
                }
                else{
                    throw CoolProp::ValueError(format("Unable to find gaseous density for T: %g K, p: %g Pa", T, p));
                }
            }
            else if (imposed_phase_index == iphase_liquid || imposed_phase_index == iphase_supercritical_liquid){
                rho = rho2;
            }
            else{
                throw ValueError("Specified phase is invalid");
            }
        }
        else{
            if (p < p_critical()){
                add_transient_pure_state();
                transient_pure_state->update(PQ_INPUTS, p, 0);
                if (T > transient_pure_state->T()){
                    double rhoV = transient_pure_state->saturated_vapor_keyed_output(iDmolar);
                    // Gas
                    if (rho0 > 0 && rho0 < rhoV){
                        rho = rho0;
                    }
                    else if( rho1 > 0 && rho1 < rhoV){
                        rho = rho1;
                    }
                    else{
                        throw CoolProp::ValueError(format("Unable to find gaseous density for T: %g K, p: %g Pa", T, p));
                    }
                }
                else{
                    // Liquid
                    rho = rho2;
                }
            }
            else{
                throw ValueError("Cubic has three roots, but phase not imposed and guess density not provided");
            }
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

CoolPropDbl CoolProp::AbstractCubicBackend::calc_molar_mass(void)
{
    double summer = 0;
    for (unsigned int i = 0; i < N; ++i)
    {
        summer += mole_fractions[i] * components[i].molemass;
    }
    return summer;
}

void CoolProp::AbstractCubicBackend::set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string &parameter, const double value){
    if (parameter == "kij" || parameter == "k_ij"){
        get_cubic()->set_kij(i, j, value);
    }
    else{
        throw ValueError(format("I don't know what to do with parameter [%s]", parameter.c_str()));
    }
    for (std::vector<shared_ptr<HelmholtzEOSMixtureBackend> >::iterator it = linked_states.begin(); it != linked_states.end(); ++it) {
        (*it)->set_binary_interaction_double(i,j,parameter,value);
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
    get_cubic()->set_kmat(donor->get_cubic()->get_kmat());
    for (std::vector<shared_ptr<HelmholtzEOSMixtureBackend> >::iterator it = linked_states.begin(); it != linked_states.end(); ++it) {
        AbstractCubicBackend *ACB = static_cast<AbstractCubicBackend *>(it->get());
        ACB->copy_k(this);
    }
}

void CoolProp::AbstractCubicBackend::set_C_MC(double c1, double c2, double c3)
{
    get_cubic()->set_C_MC(c1, c2, c3);
    for (std::vector<shared_ptr<HelmholtzEOSMixtureBackend> >::iterator it = linked_states.begin(); it != linked_states.end(); ++it) {
        AbstractCubicBackend *ACB = static_cast<AbstractCubicBackend *>(it->get());
        ACB->set_C_MC(c1,c2,c3);
    }
}

void CoolProp::AbstractCubicBackend::set_C_Twu(double L, double M, double N){
    get_cubic()->set_C_Twu(L, M, N);
    for (std::vector<shared_ptr<HelmholtzEOSMixtureBackend> >::iterator it = linked_states.begin(); it != linked_states.end(); ++it) {
        AbstractCubicBackend *ACB = static_cast<AbstractCubicBackend *>(it->get());
        ACB->set_C_Twu(L,M,N);
    }
}

void CoolProp::AbstractCubicBackend::set_fluid_parameter_double(const size_t i, const std::string parameter, const double value)
{
	// Set the volume translation parrameter, currently applied to the whole fluid, not to components.
	if (parameter == "c" || parameter == "cm" || parameter == "c_m") {
		get_cubic()->set_cm(value);
		if (this->SatL.get() != NULL) { this->SatL->set_fluid_parameter_double(i, "cm", value); }
		if (this->SatV.get() != NULL) { this->SatV->set_fluid_parameter_double(i, "cm", value); }
	}
	else {
		throw ValueError(format("I don't know what to do with parameter [%s]", parameter.c_str()));
	}
}
