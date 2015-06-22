

#include "Backends/Helmholtz/HelmholtzEOSBackend.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Backends/REFPROP/REFPROPMixtureBackend.h"
#include "Backends/REFPROP/REFPROPBackend.h"
#include <time.h>
#include "AbstractState.h"
#include "DataStructures.h"
#include <cstdio>
#include "CoolProp.h"
using namespace CoolProp;

#include "rapidjson/rapidjson_include.h"
#include "Backends/Helmholtz/Fluids/FluidLibrary.h"
#ifdef ENABLE_CATCH
#include "Tests.h"
#endif
#include "SpeedTest.h"
#include "HumidAirProp.h"
#include "time.h"
#include "Helmholtz.h"

#include "crossplatform_shared_ptr.h"

//#include <vld.h>
/*
void generate_melting_curve_data(const char* file_name, const char *fluid_name, double Tmin, double Tmax)
{

    FILE *fp;
    fp = fopen(file_name,"w");
    shared_ptr<AbstractState> State(AbstractState::factory(std::string("REFPROP"),std::string(fluid_name)));
    for (double T = Tmin; T < Tmax; T += 0.1)
    {
        try{
            double pp = State->calc_melt_p_T(T);
            State->update(PT_INPUTS,pp,T);
            double rho = State->rhomolar();
            State->update(DmolarT_INPUTS,rho,T);
            //double pp2 = State->p();
            //printf("%g,%g,%g\n",T,pp,rho);
            fprintf(fp, "%g,%g,%g\n",T,pp,rho);
        }
        catch(std::exception &e)
        {

            std::cout << fluid_name << " " << e.what() << std::endl;
            break;
        }
    }
    fclose(fp);
}
struct element
{
    double d,t,ld;
    int l;
};
*/
int main()
{
	#if 0
        double T0 = 295;
        double dT = 0.00000001;
        double pV = CoolProp::PropsSI("P","T",T0,"Q",1,"Water");    
        double hV = CoolProp::PropsSI("Hmolar","T",T0 + dT,"P",pV,"Water");
		std::cout << get_global_param_string("errstring");
        double TV = CoolProp::PropsSI("T","Hmolar",hV,"P",pV,"Water");
        std::cout << get_global_param_string("errstring");

        double ahV = CoolProp::PropsSI("P", "Q", 0, "T", 373.124, "Water");

		std::cout << get_global_param_string("parameter_list") << std::endl;
		
    #endif	
    #if 0
		
		double pm = CoolProp::PropsSI("pmax","T",-100,"Q",0.0000000000000000e+00,"REFPROP::Water");
		double Tm = PropsSI("Tmax","T",-100,"Q",0.0000000000000000e+00,"REFPROP::Water");
		double Tt = PropsSI("P","T",-100,"Q",0.0000000000000000e+00,"AceticAcid");
		double T = PropsSI("P","T",5.9020000000000005e+02,"Q",0.0000000000000000e+00,"AceticAcid");
		std::cout << get_global_param_string("errstring");
		double Tc = Props1SI("Water","T_triple");
		std::cout << get_global_param_string("errstring");
		double rhoc = Props1SI("Water","rhocrit");
		double pc = Props1SI("Water","pcrit");
		std::cout << Tc << rhoc << std::endl;
	
	#endif
    #if 0
        double T11 = PropsSI("T","P",101325,"Q",0, "Propane");
        double h3 = PropsSI("Hmolar","T",T11,"Q",0.5, "Propane");
        double T3 = PropsSI("T","P",101325,"Hmolar",h3, "Propane");
        double h4 = PropsSI("Hmolar","T",T3,"P",101325, "Propane");
        
        CoolProp::set_reference_stateS("Propane","NBP");
        double h0 = PropsSI("H","P",101325,"T",300, "Propane");
        double T0 = PropsSI("T","P",101325,"H",h0, "Propane");
        double h1 = PropsSI("H","T",T0,"P",101325, "Propane");
        
        double s0 = PropsSI("S","P",101325,"T",300, "Propane");
        double _T0 = PropsSI("T","P",101325,"S",s0, "Propane");
        double s1 = PropsSI("S","T",T0,"P",101325, "Propane");
        int r = 0;
    #endif
    #if 0
        std::string eos = get_BibTeXKey(std::string("R152A"),std::string("EOS"));
        std::cout << eos << std::endl;
        
        double p0 = PropsSI("C","T",300+273.15,"D",1e-10, "REFPROP::R152A");
        double pa = PropsSI("C","T",300+273.15,"D",1e-10, "HEOS::R152A");
        double p1 = PropsSI("P","T",273.15,"D",1,"REFPROP::R123");
        double p2 = PropsSI("P","T",273.15,"D",1,"HEOS::R123");
        std::cout << p1 << std::endl;
        std::cout << p2 << std::endl;
        double errd = p1 - p2;
        
        double l1 = PropsSI("D","T",368,"Q",1,"REFPROP::Propane");
        double l2 = PropsSI("D","T",368,"Q",1,"HEOS::Propane");
        double errl = l1 - l2;
        
        int err = 0;
    #endif
    #if 0
        shared_ptr<CoolProp::AbstractState> AS(AbstractState::factory("HEOS","Water"));
        AS->update(CoolProp::QT_INPUTS, 0, 300);
        double h1 = AS->hmass();
        shared_ptr<CoolProp::AbstractState> AR(AbstractState::factory("REFPROP","Air"));
        AR->update(CoolProp::QT_INPUTS, 0, 100);
        double h2 = AR->hmass();

        int rr = 0;
    #endif
    //if (1){
    //    shared_ptr<CoolProp::AbstractState> AS(AbstractState::factory("INCOMP","ExamplePure"));
    //    AS->update(CoolProp::PT_INPUTS, 101325, 373);
    //    double mu = AS->conductivity();
    //    int rr =0;
    //}
//    if (1)
//    {
//        std::string s = get_csv_parameter_list();
//        std::vector<std::string> keys = strsplit(s,',');
//        for (std::vector<std::string>::iterator it = keys.begin(); it != keys.end(); ++it){
//            std::string pp = (*it);
//            int key = CoolProp::get_parameter_index((*it));
//            std::string IO = get_parameter_information(key, "IO");
//            std::string units = get_parameter_information(key, "units");
//            std::string longg = get_parameter_information(key, "long");
//            std::cout << format("%s\t%s\t%s\t%s",(*it).c_str(), IO.c_str(), units.c_str(), longg.c_str()) << std::endl;
//        }
//        std::cout << s << std::endl;
//        int rr =0;
//    }
    CoolProp::set_debug_level(0);

    #if 0
		std::cout << get_global_param_string("incompressible_list_pure");
		std::cout << get_global_param_string("incompressible_list_solution");
		double rr2 = CoolProp::PropsSI("T","P",101325,"Q",0,"WATER");
		double rra = CoolProp::PropsSI("C","T",300,"D",1e-10,"HEOS::R1234ze(E)");
		double rrv = CoolProp::PropsSI("C","T",300,"D",1e-10,"REFPROP::R1234ze");
		
		std::cout << get_global_param_string("errstring");
        double rr0 = HumidAir::HAPropsSI("S","T",473.15,"W",0,"P",1e6);
        //CoolProp::set_reference_stateS("Air","RESET");
        double rr1 = HumidAir::HAPropsSI("B","T",473.15,"W",0.5,"P",101325);
        int r = 1;
    #endif
    #if 0
        // First type (slowest, most string processing, exposed in DLL)
        double r0A = PropsSI("Dmolar","T",298,"P",1e5,"Propane[0.5]&Ethane[0.5]"); // Default backend is HEOS
        double r0B = PropsSI("Dmolar","T",298,"P",1e5,"HEOS::Propane[0.5]&Ethane[0.5]");
        double r0C = PropsSI("Dmolar","T",298,"P",1e5,"REFPROP::Propane[0.5]&Ethane[0.5]");

        std::vector<double> z(2,0.5);
        // Second type (C++ only, a bit faster)
        double r1A = PropsSI("Dmolar","T",298,"P",1e5,"Propane&Ethane", z);
        double r1B = PropsSI("Dmolar","T",298,"P",1e5,"HEOS::Propane&Ethane", z);
        double r1C = PropsSI("Dmolar","T",298,"P",1e5,"REFPROP::Propane&Ethane", z);

        //const double *pz = &(z[0]);
        //int n = z.size();
        //// Third type (DLL)
        //double r2A = PropsSIZ("Dmolar","T",298,"P",1e5,"Propane&Ethane", pz, n);
        //double r2B = PropsSIZ("Dmolar","T",298,"P",1e5,"HEOS::Propane&Ethane", pz, n);
        //double r2C = PropsSIZ("Dmolar","T",298,"P",1e5,"REFPROP::Propane&Ethane", pz, n);

        double tt = 0;
    #endif

    #if 0
        shared_ptr<CoolProp::AbstractState> AS(AbstractState::factory("HEOS","ETHANE&PROPANE"));
        std::vector<double>x(2,0.5);
        AS->set_mole_fractions(x);
        //AS->build_phase_envelope("");
        AS->update(PQ_INPUTS,648000,0);
        for (int Q = 0; Q <= 1; Q++)
        {
            std::vector<double> TT,PP,RR;
            for (double p = 101325; ;p*=1.00005)
            {
                try{
                AS->update(PQ_INPUTS,p,Q);
                double T = AS->T();
                double rho = AS->rhomolar();
                double pp = AS->p();
                TT.push_back(T);
                PP.push_back(pp);
                RR.push_back(rho);
                printf("%g, %g, %g\n", T, rho, pp);
                }
                catch(std::exception &e){std::cout << e.what() << std::endl; break;}
                catch(...){break;}
            }
            rapidjson::Document d;
            d.SetObject();
            cpjson::set_double_array("T",TT,d,d);
            cpjson::set_double_array("P",PP,d,d);
            cpjson::set_double_array("rho",RR,d,d);
            std::string fname = format("%d.json",Q);
            FILE *fp = fopen(fname.c_str(),"w");
            fprintf(fp,"%s",cpjson::json2string(d).c_str());
            fclose(fp);
        }

        double rr = 0;
        return 0;
    #endif
    #if 0
    {
        std::string NBP_refs[] = {"D5","D6","MD2M","MDM","Benzene","Helium","Ethylene","Ethanol","n-Dodecane","Benzene","n-Undecane","Neon","Fluorine","Methanol","Acetone","Methane","Ethane","n-Pentane","n-Hexane","n-Heptane","n-Octane","CycloHexane","MD3M","MM","D4","MethylPalmitate","MethylStearate","MethylOleate","MethylLinoleate","MethylLinolenate","m-Xylene"};
        std::string IIR_refs[] = {"SES36","R143a","CycloPropane","Propylene","R227EA","R365MFC","R161","HFE143m","SulfurHexafluoride","CarbonDioxide","R1234ze(E)","R22","R124","Propyne","R507A","R152A","R123","R11","n-Butane","IsoButane","RC318","R21","R114","R13","R12","R113","R1233zd(E)","R41"};
        for (std::size_t i = 0; i < sizeof(NBP_refs)/sizeof(NBP_refs[0]); ++i)
        {
            try{
                //set_reference_stateS(NBP_refs[i],"RESET");
                std::vector<std::string> comps(1,NBP_refs[i]);
                HelmholtzEOSMixtureBackend HEOS(comps);
                HEOS.update(PQ_INPUTS, 101325, 0);
                double delta_a1 = HEOS.smass()/(HEOS.gas_constant()/HEOS.molar_mass());
                double delta_a2 = -HEOS.hmass()/(HEOS.gas_constant()/HEOS.molar_mass()*HEOS.get_reducing().T);
                std::cout << format("%s,%s,%16.15g,%16.15g\n",NBP_refs[i].c_str(),"NBP",delta_a1, delta_a2);
            }
            catch(const std::exception &e)
            {
                std::cout << "ERROR FOR " << NBP_refs[i] << std::endl;
            }
        }
        for (std::size_t i = 0; i < sizeof(IIR_refs)/sizeof(IIR_refs[0]); ++i)
        {
            try{
                set_reference_stateS(IIR_refs[i],"RESET");
                std::vector<std::string> comps(1,IIR_refs[i]);
                HelmholtzEOSMixtureBackend HEOS(comps);
                HEOS.update(QT_INPUTS, 0, 273.15);
                double delta_a1 = (HEOS.smass()-1000)/(HEOS.gas_constant()/HEOS.molar_mass());
                double delta_a2 = -(HEOS.hmass()-200000)/(HEOS.gas_constant()/HEOS.molar_mass()*HEOS.get_reducing().T);
                std::cout << format("%s,%s,%16.15g,%16.15g\n",IIR_refs[i].c_str(),"IIR",delta_a1, delta_a2);
            }
            catch(const std::exception &e)
            {
                std::cout << "ERROR FOR " << IIR_refs[i] << std::endl;
            }
        }
        std::string OTH_refs[] = {"Ammonia","Argon","R14"};
        for (std::size_t i = 0; i < sizeof(OTH_refs)/sizeof(OTH_refs[0]); ++i)
        {
            try{
                set_reference_stateS(OTH_refs[i],"RESET");
                std::vector<std::string> comps(1,OTH_refs[i]);
                HelmholtzEOSMixtureBackend HEOS(comps);
                REFPROPBackend REOS(OTH_refs[i]);
                HEOS.update(PQ_INPUTS, 101325, 0);
                REOS.update(PQ_INPUTS, 101325, 0);
                double delta_a1 = (HEOS.smass()-REOS.smass())/(HEOS.gas_constant()/HEOS.molar_mass());
                double delta_a2 = -(HEOS.hmass()-REOS.hmass())/(HEOS.gas_constant()/HEOS.molar_mass()*HEOS.get_reducing().T);
                std::cout << format("%s,%s,%16.15g,%16.15g\n",OTH_refs[i].c_str(),"OTH",delta_a1, delta_a2);
            }
            catch(const std::exception &e)
            {
                std::cout << "ERROR FOR " << OTH_refs[i] << ": " << e.what() << std::endl;
            }
        }

        std::string OTH2_refs[] = {"Air"};
        for (std::size_t i = 0; i < sizeof(OTH2_refs)/sizeof(OTH2_refs[0]); ++i)
        {
            try{
                //set_reference_stateS(OTH2_refs[i],"RESET");
                std::vector<std::string> comps(1,OTH2_refs[i]);
                HelmholtzEOSMixtureBackend HEOS(comps);
                REFPROPBackend REOS(OTH2_refs[i]);
                HEOS.update(QT_INPUTS, 0, 100);
                REOS.update(QT_INPUTS, 0, 100);
                double delta_a1 = (HEOS.smass()-REOS.smass())/(HEOS.gas_constant()/HEOS.molar_mass());
                double delta_a2 = -(HEOS.hmass()-REOS.hmass())/(HEOS.gas_constant()/HEOS.molar_mass()*HEOS.get_reducing().T);
                std::cout << format("%s,%s,%16.15g,%16.15g\n",OTH2_refs[i].c_str(),"OTH",delta_a1, delta_a2);
            }
            catch(const std::exception &e)
            {
                std::cout << "ERROR FOR " << OTH2_refs[i] << ": " << e.what() << std::endl;
            }
        }
        double rr = 0;
    }
	#endif
    #if 0
    {
        generate_melting_curve_data("Ethylene-I.mlt","ethylene",103.989,110.369);
        generate_melting_curve_data("Ethylene-II.mlt","ethylene",110.369,450);
        generate_melting_curve_data("Propylene-I.mlt","propylen",87.953,109.6);
        generate_melting_curve_data("Propylene-II.mlt","propylen",109.6,575);
        generate_melting_curve_data("ParaHydrogen-I.mlt","parahyd",13.8033,22);
        generate_melting_curve_data("ParaHydrogen-II.mlt","parahyd",22,2000);

        generate_melting_curve_data("n-Propane.mlt","propane",85.53,2000);
        generate_melting_curve_data("n-Butane.mlt","butane",134.9,2000);
        generate_melting_curve_data("n-Pentane.mlt","pentane",143.5,2000);
        generate_melting_curve_data("IsoButane.mlt","isobutan",113.73,2000);
        generate_melting_curve_data("Isopentane.mlt","ipentane",112.66,2000);
        generate_melting_curve_data("Argon.mlt","argon",83.8058,2000);
        generate_melting_curve_data("Ethane.mlt","ethane",90.37,2000);
        generate_melting_curve_data("Nitrogen.mlt","nitrogen",63.151,2000);
        generate_melting_curve_data("Fluorine.mlt","fluorine",53.4811,2000);
        generate_melting_curve_data("Methane.mlt","methane",90.70,2000);
        generate_melting_curve_data("Methanol.mlt","methanol",175.61,2000);
        generate_melting_curve_data("Krypton.mlt","krypton",115.775,2000);
        generate_melting_curve_data("Xenon.mlt","xenon",161.405,2000);
        generate_melting_curve_data("CarbonMonoxide.mlt","co",68.16,2000);
        generate_melting_curve_data("Oxygen.mlt","oxygen",54.361,2000);
        generate_melting_curve_data("CycloHexane.mlt","cyclohex",279.7,2000);
        generate_melting_curve_data("CarbonDioxide.mlt","CO2",217,2000);
    }
	#endif
    #if 0
    {
        std::cout << "Water DmolarT at 1e-3,300 \n-----------------\n";
        CoolProp::compare_REFPROP_and_CoolProp("Water",DmolarT_INPUTS, 1e-3, 300, 10000, 0, 1e-8);
        std::cout << "Water PT at 101325,300 \n-----------------\n";
        CoolProp::compare_REFPROP_and_CoolProp("Water",PT_INPUTS, 101325, 300, 10000, 0, 1e-8);
        std::cout << "Water QT at 350 K\n-----------------\n";
        CoolProp::compare_REFPROP_and_CoolProp("Water",QT_INPUTS, 1, 350, 10000, 0, 1e-8);
        std::cout << "Water PQ at 101325 Pa\n-----------------\n";
        CoolProp::compare_REFPROP_and_CoolProp("Water",PQ_INPUTS, 101325, 1, 10000, 1e-2, 0);

        std::cout << "R134a DmolarT at 1e-3,300 \n-----------------\n";
        CoolProp::compare_REFPROP_and_CoolProp("R134a",DmolarT_INPUTS, 1e-3, 300, 10000, 0, 1e-8);
        std::cout << "R134a PT at 101325,300 \n-----------------\n";
        CoolProp::compare_REFPROP_and_CoolProp("R134a",PT_INPUTS, 101325, 300, 10000, 0, 1e-8);
        std::cout << "R134a QT at 350 K\n-----------------\n";
        CoolProp::compare_REFPROP_and_CoolProp("R134a",QT_INPUTS, 1, 350, 10000, 0, 1e-8);
        std::cout << "R134a PQ at 101325 Pa\n-----------------\n";
        CoolProp::compare_REFPROP_and_CoolProp("R134a",PQ_INPUTS, 101325, 1, 10000, 1e-2, 0);
    }
	#endif
    #if 0
    {
        std::vector<std::string> ss = strsplit(get_global_param_string("FluidsList"),',');

        for (std::vector<std::string>::iterator it = ss.begin(); it != ss.end(); ++it)
        {
            AbstractState *S = AbstractState::factory("HEOS", (*it));
            S->update(QT_INPUTS, 0, S->Ttriple());
            std::cout << format("%s %17.15g\n", S->name().c_str(), S->p());
        }
    }
	#endif
	#if 0
	{
        shared_ptr<CoolProp::AbstractState> ASR(CoolProp::AbstractState::factory("HEOS","H2S"));
        ASR->update(PT_INPUTS, 1000e6, 200);
        double v  = ASR->viscosity();
        int rr =0;
    }
	#endif
    #if 0
    {
        shared_ptr<CoolProp::AbstractState> ASR(CoolProp::AbstractState::factory("REFPROP","CO2"));
        double p1 = ASR->calc_melt_p_T(250);

        shared_ptr<CoolProp::AbstractState> ASC(CoolProp::AbstractState::factory("HEOS","CO2"));
        double p2 = ASC->calc_melt_p_T(250);

        double rrr1 = PropsSI("D","T",200,"P",300,"Water");
        double rrr0 = PropsSI("Cvmolar","T",200,"Dmolar",14000,"REFPROP::R125");
        double rrr2 = PropsSI("speed_of_sound","T",300,"Dmolar",700,"R125");
        double rrr =0 ;
    }
	#endif
    #if 0
    {
        shared_ptr<AbstractState> ASR(AbstractState::factory("REFPROP","CO2"));
        ASR->update(QT_INPUTS, 1, 304);
        double muR0 = ASR->conductivity();

        shared_ptr<AbstractState> ASC(AbstractState::factory("HEOS","CO2"));
        ASC->update(QT_INPUTS, 1, 304);
        double muC = ASC->conductivity();
        double rr = 4;
    }
	#endif
    #if 0
    {
        double Tc = Props1SI("Water","Tcrit");
        double pc = Props1SI("Water","pcrit");
        double p = pc*2;
        double T = Tc*0.5;
        char ykey[] = "H";
        double y = PropsSI(ykey,"P",p,"T",T,"Water");
        double TT = PropsSI("T","P",p,ykey,y,"Water");
        int rr = 0;
    }
    #endif
    #if 0
    {
        double h = HumidAir::HAPropsSI("H","T",303.15,"R",1.0000000000000000e+00,"P",1.0132500000000000e+05 );
        double T = HumidAir::HAPropsSI("T","H",h,"R",1.0000000000000000e+00,"P",1.0132500000000000e+05 );
        double hh = HumidAir::HAPropsSI("H","T",T,"R",1.0000000000000000e+00,"P",1.0132500000000000e+05 );
        double s = HumidAir::HAPropsSI("S","T",2.1814999999999998e+02,"R",1.0000000000000000e+00,"P",1.0132500000000000e+05);
        int r = 3;
    }
    #endif
    #if 0
    {
        ::set_debug_level(0);
        std::vector<std::string> tags;
        tags.push_back("[mixture_derivs]");
        run_user_defined_tests(tags);
        char c;
        std::cin >> c;
    }
    #endif
    #if 0
    {
		run_tests();
		char c;
		std::cin >> c;
	}
	#endif
    #if 0
    {
        double TTT0 = PropsSI("T","Q",1,"P",3e6,"REFPROP::R32[0.3]&R125[0.7]");
        double TTT1 = PropsSI("T","Q",1,"P",3e6,"HEOS::R125[0.7]&R32[0.3]");
        int rr =0;
    }
    #endif
	#if 1
    double T1 = HumidAir::HAPropsSI("T", "P", 101325, "H", 202520, "W", 0.013);
    double T2 = HumidAir::HAPropsSI("T", "P", 101325, "H", 202520, "W", 0.007);
    double T3 = HumidAir::HAPropsSI("T", "P", 101325, "H", 202520, "W", 0.0000001);
	shared_ptr<CoolProp::AbstractState> CP(CoolProp::AbstractState::factory("HEOS", "Water"));
	shared_ptr<CoolProp::AbstractState> RP(CoolProp::AbstractState::factory("REFPROP", "Water"));
	std::vector<std::string> fluids = RP->fluid_names();
	double tt = RP->rhomolar_critical();
	int rr = 3l;
	#endif
    #if 0
    { 
		//double B1 = HumidAir::HAPropsSI("B","T",283.73,"W",0,"P",101325);// = 283,73
		//double B2 = HumidAir::HAPropsSI("B","T",193.92,"W",0,"P",101325);// = 193,92
		//double B3 = HumidAir::HAPropsSI("B","T",194.21,"W",0,"P",101325);// = 193,92
		shared_ptr<CoolProp::AbstractState> CP(CoolProp::AbstractState::factory("HEOS", "Water"));
		shared_ptr<CoolProp::AbstractState> RP(CoolProp::AbstractState::factory("REFPROP", "Water"));
		RP->update(PT_INPUTS, 101325, 300);
		CP->update(PT_INPUTS, 101325, 300);

		double are1 = RP->molar_mass();
		double are14 = CP->molar_mass();
		
		
		//double ar1 = RP->first_partial_deriv(iHmass,iT,iP);
		//double ar2 = CP->first_partial_deriv(iHmass,iT,iP);
		//double a01 = RP->cpmass();
		//double a02 = CP->cpmass();
		
		std::cout << format("alphar %Lg %Lg\n", RP->alphar(),               CP->alphar());
		std::cout << format("dalphar_dDelta %Lg %Lg\n", RP->dalphar_dDelta(),       CP->dalphar_dDelta());
		std::cout << format("dalphar_dTau %Lg %Lg\n", RP->dalphar_dTau(),         CP->dalphar_dTau());
		std::cout << format("d2alphar_dDelta2 %Lg %Lg\n", RP->d2alphar_dDelta2(),     CP->d2alphar_dDelta2());
		std::cout << format("d2alphar_dDelta_dTau %Lg %Lg\n", RP->d2alphar_dDelta_dTau(), CP->d2alphar_dDelta_dTau());
		std::cout << format("d2alphar_dTau2 %Lg %Lg\n", RP->d2alphar_dTau2(),       CP->d2alphar_dTau2());
		
		std::cout << format("alpha0 %Lg %Lg\n", RP->alpha0(),               CP->alpha0());
		std::cout << format("dalpha0_dDelta %Lg %Lg\n", RP->dalpha0_dDelta(),       CP->dalpha0_dDelta());
		std::cout << format("dalpha0_dTau %Lg %Lg\n", RP->dalpha0_dTau(),         CP->dalpha0_dTau());
		std::cout << format("d2alpha0_dDelta2 %Lg %Lg\n", RP->d2alpha0_dDelta2(),     CP->d2alpha0_dDelta2());
		std::cout << format("d2alpha0_dDelta_dTau %Lg %Lg\n", RP->d2alpha0_dDelta_dTau(), CP->d2alpha0_dDelta_dTau());
		std::cout << format("d2alpha0_dTau2 %Lg %Lg\n", RP->d2alpha0_dTau2(),       CP->d2alpha0_dTau2());
		
		double T1 = HumidAir::HAPropsSI("T", "V", 0.83, "R", 1, "P", 101325);// = 283,73
		double T32 = HumidAir::HAPropsSI("T","B",273.15,"W",0,"P",101325);// = 193,92
		double T2 = HumidAir::HAPropsSI("T","B",273.16,"W",0,"P",101325);// = 193,92
		double T3 = HumidAir::HAPropsSI("T","B",273.40,"W",0,"P",101325);// = 193,92
		double dd8 = PropsSI("D","T",8.5525000000000006e+01,"Q",0,"REFPROP::Propane");
		
        ::set_debug_level(0);
		::set_debug_level(0);
        
        shared_ptr<AbstractState> HEOS(AbstractState::factory("HEOS","Methane&Ethane"));
        std::vector<long double> z(2, 0.8); z[1] = 1-z[0];
        //shared_ptr<AbstractState> HEOS(AbstractState::factory("HEOS","Methane&Propane&Ethane&n-Butane"));
        //std::vector<long double> z(4, 0.1); z[1] = 0.35; z[2] = 0.35, z[3] = 0.2;
        HEOS->set_mole_fractions(z);
        
        time_t t1, t2;
        t1 = clock();
        try{
            HEOS->build_phase_envelope("dummy");
        }
        catch(std::exception &e){
            std::cout << get_global_param_string("errstring") << std::endl;
        }
        t2 = clock();
        HEOS->update(PSmolar_INPUTS, 4e6, 79.1048486373);
        long double TT = HEOS->T();
        HEOS->update(PQ_INPUTS, 1.3e5, 1);
        double ssat = HEOS->smolar();
        double hsat = HEOS->hmolar();
        double dsat = HEOS->rhomolar();
        for (long double s = ssat + 60; s > ssat; s -= 5){
            HEOS->update(PSmolar_INPUTS, 1.3e5, s);
            std::cout << s << " " << HEOS->rhomolar() << " " << dsat << std::endl;
        }
        std::cout << format("time: %g s/call\n", ((double)(t2-t1))/CLOCKS_PER_SEC);
        exit(EXIT_SUCCESS);
        
        std::cout << get_global_param_string("errstring") << std::endl;
        exit(EXIT_FAILURE);
        //double refretrte = PropsSI("P","Dmolar",107.9839357,"T",116.5360225,"Methane[0.5]&Propane[0.5]");
        for (double p = 101325; p < 9e6; p *= 1.05){
            std::cout << p << " " << PropsSI("T","P",p,"Q",1,"Methane[0.5]&Propane[0.5]") << std::endl;
            //std::cout << get_global_param_string("errstring") << std::endl;
        }
        int rr =1;
    }
    #endif
    #if 0
    {        
        ::set_debug_level(6);
        double dd6 = PropsSI("P","T",300,"Q",0,"TTSE&HEOS::n-Propane");
    }
    #endif
	#if 0
	{
        
        ::set_debug_level(0);
        
        shared_ptr<AbstractState> Water(AbstractState::factory("HEOS","water"));
        Water->update(PT_INPUTS, 800,  300);
        
        Water->update(HmolarP_INPUTS, 45960.1, 800);
        CoolProp::phases phase = Water->phase();
        double tt = Water->T();
        
        double ee6 = PropsSI("T","P",101325,"Q",0,"Water");
		char ykey[] = "H";
		double Ts, y, T2, dT = -1;
		double dd0 = CoolProp::Props1SI("Tmax","n-Propane");
        double dd1 = CoolProp::Props1SI("n-Propane","Tmax");
        double Tc = CoolProp::Props1SI("n-Propane","Tcrit");
        
        double dd6 = PropsSI("P","T",Tc-1e-4,"Q",0,"TTSE&HEOS::n-Propane");
        
        double dd7 = PropsSI("P","T",Tc-1e-4,"Q",0,"n-Propane");
        std::cout << get_global_param_string("errstring");
        
        double dd8 = PropsSI("D","T",-8.5525000000000006e+01,"Q",0,"REFPROP::Propane");
        std::cout << get_global_param_string("errstring");
        double dd9 = PropsSI("T","U",3.7175480877288617e+05,"P",1.7372110031440207e-04,"n-Propane");
        std::cout << get_global_param_string("errstring");
        
		
        Water->update(PT_INPUTS, 101325, 0);
		double ptt = Water->melting_line(iT, iP, 138.268e6);
        
        Water->update(QT_INPUTS, 0.5, 300);
        double hmolar = Water->hmolar();
        double p = Water->p();
        Water->update(HmolarP_INPUTS, hmolar, p);
        double T = Water->T();
		
		std::cout << get_global_param_string("errstring");
		y = PropsSI(ykey,"T",Ts+dT,"P",101325,"n-Propane");
		T2 = PropsSI("T",ykey,y,"P",101325,"n-Propane");
		std::cout << get_global_param_string("errstring");
		
        #if ENABLE_CATCH
            std::vector<std::string> tags;
            tags.push_back("[flash]");
            run_user_defined_tests(tags);
            double rr = 0;
			char c;
		std::cin >> c;
        #endif
    }
	#endif
    #if 0
    {
        time_t t1,t2;
        long N = 100000;
        double ss = 0;
        std::vector<std::string> names(1,"Propane");
        
        shared_ptr<HelmholtzEOSMixtureBackend> Water(new HelmholtzEOSMixtureBackend(names));
        Water->set_mole_fractions(std::vector<long double>(1,1));
        ResidualHelmholtzGeneralizedExponential GenExp = Water->get_components()[0]->pEOS->alphar.GenExp;
        
        HelmholtzDerivatives derivs1, derivs2;
        double tau = 0.8, delta = 1.1;
        
        ss = 0;
        t1 = clock();
        for (long i = 0; i < N; ++i){
            derivs1.reset();
            GenExp.all(tau, delta+i*1e-10, derivs1);
            ss += derivs1.alphar;
        }
        t2 = clock();
        std::cout << format("value(all): %0.13g, %0.13g, %g us/call\n", ss, derivs1.alphar, ((double)(t2-t1))/CLOCKS_PER_SEC/double(N)*1e6);
        
        ss = 0;
        t1 = clock(); 
        for (long i = 0; i < N; ++i){
            derivs2.reset();
            GenExp.allEigen(tau, delta+i*1e-10, derivs2);
            ss += derivs2.alphar;
        }
        t2 = clock();
        std::cout << format("value(allEigen): %0.13g, %0.13g, %g us/call\n", ss, derivs2.alphar, ((double)(t2-t1))/CLOCKS_PER_SEC/double(N)*1e6);

        int r44 =0;
    }
    #endif
    #if 0
    {
        
        t1 = clock();
        for (long i = 0; i < N; ++i){
            Water->update(PQ_INPUTS, 10132+i, 0);
            ss += Water->p();
        }
        t2 = clock();
        std::cout << format("value: %0.13g, %g us/call\n", ss, ((double)(t2-t1))/CLOCKS_PER_SEC/double(N)*1e6);
        ss = 0;
        {
            shared_ptr<REFPROPMixtureBackend> Water(new REFPROPMixtureBackend(names));
            Water->set_mole_fractions(std::vector<long double>(1,1));
            
            t1 = clock();
            for (long i = 0; i < N; ++i){
                Water->update(PQ_INPUTS, 10132+i, 0);
                ss += Water->p();
            }
            t2 = clock();
            std::cout << format("value: %0.13g, %g us/call\n", ss, ((double)(t2-t1))/CLOCKS_PER_SEC/double(N)*1e6);
        }
    }
    #endif
    #if 0
    {
        std::vector<std::string> names(1,"Water");
        shared_ptr<HelmholtzEOSMixtureBackend> Water(new HelmholtzEOSMixtureBackend(names));
        Water->set_mole_fractions(std::vector<long double>(1,1));
        Water->update(PQ_INPUTS, 101325, 0);
        
        HelmholtzDerivatives derivs;
        time_t t1,t2;
        long N = 100000;
        double ss = 0;
        
        std::vector<CoolPropFluid*> components = Water->get_components();
        ResidualHelmholtzGeneralizedExponential GenExp = components[0]->pEOS->alphar.GenExp;
        ResidualHelmholtzNonAnalytic NonAnal = components[0]->pEOS->alphar.NonAnalytic;
        
        long double tau = 0.8, delta = 2.2;
        derivs.reset();
        GenExp.all(tau, delta, derivs);
        
        t1 = clock();
        for (long i = 0; i < N; ++i){
            derivs.reset();
            GenExp.all(tau, delta+i*1e-10, derivs);
            ss += derivs.alphar+derivs.dalphar_ddelta+derivs.dalphar_dtau+derivs.d2alphar_ddelta2+derivs.d2alphar_ddelta_dtau+derivs.d2alphar_dtau2;
        }
        t2 = clock();
        std::cout << format("value: %0.13g, %g us/call\n", ss, ((double)(t2-t1))/CLOCKS_PER_SEC/double(N)*1e6);
        
        tau = 0.99; delta = 1.01;
        derivs.reset();
        NonAnal.all(tau, delta, derivs);
        long double a00 = NonAnal.base(tau, delta);
        long double a10 = NonAnal.dDelta(tau, delta);
        long double a01 = NonAnal.dTau(tau, delta);
        long double a20 = NonAnal.dDelta2(tau, delta);
        long double a11 = NonAnal.dDelta_dTau(tau, delta);
        long double a02 = NonAnal.dTau2(tau, delta);
        long double a03 = NonAnal.dTau3(tau, delta);
        long double a12 = NonAnal.dDelta_dTau2(tau, delta);
        long double a21 = NonAnal.dDelta2_dTau(tau, delta);
        long double a30 = NonAnal.dDelta3(tau, delta);
        
        exit(0);
    }
    #endif
    #if 0
    {
        /*double h1 = PropsSI("S","P",101325,"Q",0,"n-Pentane");
        std::string er = get_global_param_string("errstring");
        set_reference_stateS("n-Propane","NBP");
        double h2 = PropsSI("H","P",101325,"Q",0,"n-Propane");*/

        //std::string RPname = get_fluid_param_string("Water", "REFPROPname");
        //std::string s = get_BibTeXKey("n-Propane", "rr");



        #if ENABLE_CATCH
        run_tests();
        #endif

        std::string fl = get_global_param_string("FluidsList");
        double rr = PropsSI("D", "P", 3e5, "T", 300, "Nitrogen");
        shared_ptr<AbstractState> AS(AbstractState::factory("HEOS","Nitrogen"));

        AS->update(DmolarT_INPUTS, 40, 300);
        double p1 = AS->umolar();
        double d1 = AS->rhomolar();
        double T1 = AS->delta();
        double dpdT_constrho = AS->first_partial_deriv(iUmolar, iDelta, iTau);
        AS->update(DmolarT_INPUTS, 40+1e-6, 300);
        double p2 = AS->umolar();
        double d2 = AS->rhomolar();
        double T2 = AS->delta();

        double dpdT_constrho2 = (p2-p1)/(T2-T1);

        AS->update(PT_INPUTS, 101000, 300);

        std::vector<double> T(2,300), P(2,101325), o, z(1,1);
        std::string in1 = "Dmass", in2 = "T", in3 = "P", Ref = "Nitrogen";
        T[1] = 400;
        o = PropsSI(in1,in2,T,in3,P,Ref,z);
        double tr = 0;
    }
	#endif
    #if 0
    {
        // First type (slowest, most string processing, exposed in DLL)
        double r0A = PropsSI("Dmolar","T",298,"P",1e5,"Propane[0.5]&Ethane[0.5]"); // Default backend is HEOS
        double r0B = PropsSI("Dmolar","T",298,"P",1e5,"HEOS::Propane[0.5]&Ethane[0.5]");
        double r0C = PropsSI("Dmolar","T",298,"P",1e5,"REFPROP::Propane[0.5]&Ethane[0.5]");

        std::vector<double> z(2,0.5);
        // Second type (C++ only, a bit faster)
        double r1A = PropsSI("Dmolar","T",298,"P",1e5,"Propane&Ethane", z);
        double r1B = PropsSI("Dmolar","T",298,"P",1e5,"HEOS::Propane&Ethane", z);
        double r1C = PropsSI("Dmolar","T",298,"P",1e5,"REFPROP::Propane&Ethane", z);

        //const double *pz = &(z[0]);
        //int n = z.size();
        //// Third type (DLL)
        //double r2A = PropsSIZ("Dmolar","T",298,"P",1e5,"Propane&Ethane", pz, n);
        //double r2B = PropsSIZ("Dmolar","T",298,"P",1e5,"HEOS::Propane&Ethane", pz, n);
        //double r2C = PropsSIZ("Dmolar","T",298,"P",1e5,"REFPROP::Propane&Ethane", pz, n);

        double tt = 0;
    }
	#endif
    #if 0
    {
        #ifdef ENABLE_CATCH
        run_tests();
        #endif
    }
	#endif
    #if 0
    {


        typedef double dbltype;
        dbltype n[] = {0.0125335479355233,                        7.8957634722828,                        -8.7803203303561,                        0.31802509345418,                        -0.26145533859358,-0.0078199751687981,0.0088089493102134,                        -0.66856572307965,                        0.20433810950965,                        -6.621260503968699e-005,                        -0.19232721156002,                        -0.25709043003438,                        0.16074868486251,                        -0.040092828925807,                        3.9343422603254e-007,                        -7.5941377088144e-006,                        0.00056250979351888,                        -1.5608652257135e-005,                        1.1537996422951e-009,                        3.6582165144204e-007,                        -1.3251180074668e-012,                        -6.2639586912454e-010,                        -0.10793600908932,                        0.017611491008752,                        0.22132295167546,                        -0.40247669763528,                        0.58083399985759,                        0.0049969146990806,                        -0.031358700712549,                        -0.74315929710341,                        0.4780732991548,                        0.020527940895948,                        -0.13636435110343,                        0.014180634400617,                        0.008332650488071301,                        -0.029052336009585,                        0.038615085574206,                        -0.020393486513704,                        -0.0016554050063734,                        0.0019955571979541,                        0.00015870308324157,                        -1.638856834253e-005,                        0.043613615723811,                        0.034994005463765,-0.076788197844621,0.022446277332006,-6.2689710414685e-005,-5.5711118565645e-010,-0.19905718354408,0.31777497330738,-0.11841182425981};
        dbltype d[] = {1,                         1,                        1,                        2,                        2,                        3,                        4,                        1,                        1,                        1,                        2,                        2,                        3,                        4,                        4,                        5,                        7,                        9,                        10,                        11,                        13,                        15,                        1,                        2,                        2,                        2,                        3,                        4,                        4,                        4,                        5,                        6,                        6,                        7,                        9,                        9,                        9,                        9,                        9,                        10,                        10,                        12,                        3,                        4,                        4,                        5,                        14,                        3,                        6,                        6,                        6                    };
        dbltype t[] = {-0.5,                        0.875,                        1,                        0.5,                        0.75,                        0.375,                        1,                        4,                        6,                        12,                        1,                        5,                        4,                        2,                        13,                        9,                        3,                        4,                        11,                        4,                        13,                        1,                        7,                        1,                        9,                        10,                        10,                        3,                        7,                        10,                        10,                        6,                        10,                        10,                        1,                        2,                        3,                        4,                        8,                        6,                        9,                        8,                        16,                        22,                        23,                        23,                        10,                        50,                        44,                        46,                        50                    };
        dbltype l[] = {0,                        0,                        0,                        0,                        0,                        0,                        0,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        3,                        3,                        3,                        3,                        4,                        6,                        6,                        6,                        6                    };
        dbltype summer = 0;
        std::vector<element> elements;
        for (std::size_t i = 0; i < 51; ++i)
        {
            element el;
            el.d = d[i];
            el.t = t[i];
            el.l = (int)l[i];
            el.ld = (dbltype)l[i];
            elements.push_back(el);
        }

        long N = 1000000;
        dbltype pow_delta_li, di, ti, lid;
        std::vector<dbltype> s(51);
        dbltype delta, log_tau, log_delta, tau, gamma;
        double t1 = clock();

        for (std::size_t ii = 0; ii < N; ++ii)
        {
            delta = 1.3, tau = 0.7;
            log_tau  = log(tau), log_delta = log(delta);
            pow_delta_li = pow(delta, 5);
            for (int jj = 0; jj < 51; ++jj)
            {
                //di = d[jj]; ti = t[jj]; lid = l[jj]; gamma = 1;
                //int li = (int)lid;

                element &el = elements[jj];
                ti = el.t;
                //int li = (int)lid;

                summer += exp(ti*log_tau);//(di-gamma*lid*pow_delta_li)*pow(tau,ti)*pow_delta_li*exp(-gamma*pow_delta_li);
                //summer += __ieee754_exp(ti*log_tau);//(di-gamma*lid*pow_delta_li)*pow(tau,ti)*pow_delta_li*exp(-gamma*pow_delta_li);
//                if (li > 0){
//
//
//                }
//                else{
//                    s[jj] += di*exp(ti*log_tau+(di-1)*log_delta);
//                }
            }
        }
        //summer = std::accumulate(s.begin(), s.end(), (dbltype)(0));
        double t2 = clock();
        double elap = (t2-t1)/CLOCKS_PER_SEC/((dbltype)N)*1e6;
        printf("%g %g\n",elap, summer);
        int rr =5;
    }
	#endif
    #if 0
    {
        shared_ptr<AbstractState> MixRP(AbstractState::factory(std::string("REFPROP"),std::string("propane")));
        MixRP->update(QT_INPUTS, 0, 330);
        long double s1 = MixRP->surface_tension();

        shared_ptr<AbstractState> Mix(AbstractState::factory(std::string("HEOS"), std::string("propane")));
        Mix->update(QT_INPUTS, 0, 330);
        long double s2 = Mix->surface_tension();
    }
	#endif
    #if 0
    {

        double T = 300;

        shared_ptr<AbstractState> MixRP(AbstractState::factory(std::string("REFPROP"),std::string("propane")));
        {
            long N = 100000;
            double t1 = clock(), summer = 0;
            for (std::size_t ii = 0; ii < N; ++ii)
            {
                MixRP->update(QT_INPUTS, 0.0, T+10/((double)N)*ii);
                summer += MixRP->p();
            }
            double t2 = clock();
            double elap = (t2-t1)/CLOCKS_PER_SEC/((double)N)*1e6;
            printf("%g %g\n",elap, summer);
        }
        double p2 = MixRP->p();
        double cv2 = MixRP->cvmolar();
        double cp2 = MixRP->cpmolar();
        double T2 = MixRP->T();

        shared_ptr<AbstractState> Mix(AbstractState::factory(std::string("HEOS"), std::string("propane")));
        {
            long N = 100000;
            double t1 = clock(), summer = 0;
            for (std::size_t ii = 0; ii < N; ++ii)
            {
                Mix->update(QT_INPUTS, 0.0, T+10/((double)N)*ii);
                summer += Mix->p();
            }
            double t2 = clock();
            double elap = (t2-t1)/CLOCKS_PER_SEC/((double)N)*1e6;
            printf("%g %g\n",elap, summer);
        }
        double p1 = Mix->p();
        double cv1 = Mix->cvmolar();
        double cp1 = Mix->cpmolar();
        double T1 = Mix->T();

        double rr = 0;
    }
	#endif
    #if 0
    {
        int N = 2;
        std::vector<long double> z(N, 1.0/N);
        double Q = 1, T = 250, p = 300000;

        int inputs = PQ_INPUTS; double val1 = p, val2 = Q;

        shared_ptr<AbstractState> MixRP(AbstractState::factory(std::string("REFPROP"), std::string("Ethane,propane")));
        MixRP->set_mole_fractions(z);
        MixRP->update(inputs, val1, val2);
        double p2 = MixRP->p();
        double rho2 = MixRP->rhomolar();
        double cv2 = MixRP->cvmolar();
        double cp2 = MixRP->cpmolar();
        double w2 = MixRP->speed_sound();
        double h2 = MixRP->hmolar();
        double s2 = MixRP->smolar();
        double phi20 = MixRP->fugacity_coefficient(0);
        double phi21 = MixRP->fugacity_coefficient(1);

        shared_ptr<AbstractState> Mix(AbstractState::factory(std::string("HEOS"), std::string("Ethane,propane")));
        Mix->set_mole_fractions(z);
        Mix->update(inputs, val1, val2);
        double p1 = Mix->p();
        double cv1 = Mix->cvmolar();
        double cp1 = Mix->cpmolar();
        double w1 = Mix->speed_sound();
        double h1 = Mix->hmolar();
        double s1 = Mix->smolar();
        double phi10 = Mix->fugacity_coefficient(0);
        double phi11 = Mix->fugacity_coefficient(1);

        double rr = 0;
    }
	#endif
    #if 0
    {
        int N = 2;
        std::vector<long double> z(N, 1.0/N);
        double Q = 0, T = 250, p = 300000;

        shared_ptr<AbstractState> Mix(AbstractState::factory(std::string("HEOS"), std::string("Ethane,n-Propane")));

        Mix->set_mole_fractions(z);

        for (double T = 210; ;T += 0.1)
        {
            Mix->update(QT_INPUTS, Q, T);
            std::cout << format(" %g %g\n",Mix->p(),Mix->T());
        }
    }
	#endif
    #if 0
    {
        time_t t1,t2;

        std::size_t N = 1000000;
        shared_ptr<AbstractState> State(AbstractState::factory(std::string("HEOS"), std::string("Water")));

        double p = State->p();
        double summer = 0;
        t1 = clock();
        for (std::size_t ii = 0; ii < N; ++ii)
        {
            //AbstractState *State = new REFPROPBackend("Methane");
            //summer += EOS->dalphar_dDelta(0.7,1.3);
            /*for (int i = 0; i < 50; i++)
            {
                summer += exp(1.3+1e-10*ii+1e-10*i);
            }
            summer += log(0.7-1e-10*ii);*/
            //summer += log(1.3);
            State->update(PT_INPUTS,101325,300);
            summer += State->p();
        }
        t2 = clock();
        double elap = ((double)(t2-t1))/CLOCKS_PER_SEC/((double)N)*1e6;
        printf("%g %g\n",elap, summer/((double)N));
        double eee = 0;
        return 0;
    }
	#endif

    #if 0
    {
        shared_ptr<AbstractState> State(AbstractState::factory(std::string("REFPROP"), std::string("Methane|Ethane")));

        std::vector<long double> x(2,0.5);
        State->set_mole_fractions(x);
        State->update(DmassT_INPUTS,1,250);
        double hh = State->hmolar();
        double mu = State->viscosity();
        double sigma = State->surface_tension();
    }
	#endif
    #if 0
    {
        time_t t1,t2;
        t1 = clock();
        long N = 100000;
        for (long ii = 0; ii < N; ii++)
        {
            shared_ptr<AbstractState> State(AbstractState::factory(std::string("REFPROP"), std::string("Methane")));
            //AbstractState *State = new REFPROPBackend("Methane");
        }
        t2 = clock();
        double elap = ((double)(t2-t1))/CLOCKS_PER_SEC/((double)N)*1e6;
        printf("%g\n",elap);
    }
	#endif
    #if 0
    {
        shared_ptr<AbstractState> State(AbstractState::factory(std::string("REFPROP"), std::string("Methane")));

        State->update(DmassT_INPUTS,1,300);
        double hh = State->hmolar();
        double mu = State->viscosity();

        time_t t1,t2;
        t1 = clock();
        for (long ii = 0; ii < 1000000; ii++)
        {
            State->update(PQ_INPUTS,300000,1-ii*1e-6);
            //State->update(DmassT_INPUTS,1-ii*1e-10,180);
            //double hh1 = State->hmolar();
            //double mu2 = State->viscosity();
        }
        t2 = clock();
        double elap = ((double)(t2-t1))/CLOCKS_PER_SEC;
        printf("%g\n",elap);
    }
	#endif
}
