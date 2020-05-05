using SMath.Math.Symbolic;

namespace coolprop_wrapper
{
    class Unit
    {
        static MItem unit(string Format, params string[] units)
        {
            var lUnits = new System.Collections.Generic.List<string>(units.Length);
            foreach (var unit in units)
                lUnits.Add(SMath.Manager.UnitsManager.GetCurrentUnitName(unit));
            return Converter.ToMItem(string.Format(Format, lUnits.ToArray()));
        }

        public static MItem K { get { return unit("{0}", "'K"); } }
        // [Pa] = [kg/{m*s^2}]
        public static MItem Pa { get { return unit("{0}", "'Pa"); } }
        // [J] = [kg*m^2/s^2]
        public static MItem J_kg { get { return unit("{0}/{1}", "'J", "'kg"); } } // J/kg
        public static MItem mol_m3 { get { return unit("{0}/{{{1}^3}}", "'mol","'m"); } } // mol/m^3
        public static MItem kg_m3 { get { return unit("{0}/{{{1}^3}}", "'kg", "'m"); } } // kg/m^3
        public static MItem J_mol { get { return unit("{0}/{1}", "'J", "'mol"); } } // J/mol
        public static MItem J_mol_K { get { return unit("{0}/{{{1}*{2}}}", "'J", "'mol", "'K"); } } // J/mol/K
        public static MItem J_kg_K { get { return unit("{0}/{{{1}*{2}}}", "'J", "'kg", "'K"); } } // J/kg/K
        public static MItem m_s { get { return unit("{0}/{1}", "'m", "'s"); } } // m/s
        // [W]=[{kg*m^2}/{s^3}]
        public static MItem W_m_K { get { return unit("{0}/{{{1}*{2}}}", "'W", "'m", "'K"); } } // W/m/K
        public static MItem _K { get { return unit("1/{0}", "'K"); } } // 1/K
        public static MItem _Pa { get { return unit("1/{0}", "'Pa"); } } // 1/Pa
        // [N]=[{kg*m}/{s^2}]
        public static MItem N_m { get { return unit("{0}/{1}", "'N", "'m"); } } // N/m
        public static MItem kg_mol { get { return unit("{0}/{1}", "'kg", "'mol"); } } // kg/mol
        public static MItem Pa_s { get { return unit("{0}*{1}", "'Pa", "'s"); } } // Pa*s
        public static MItem m3_kg { get { return unit("{{{0}^3}}/{1}", "'m", "'kg"); } } // m^3/kg

        public static MItem unitless { get { return unit("1"); } }

        static void AddUnits(ref System.Collections.Generic.Dictionary<string, MItem> dic, MItem unit, params string[] names)
        {
            foreach (var name in names)
                dic.Add(name, unit);
        }

        static System.Collections.Generic.Dictionary<string, MItem> InitUnitsDictionary()
        {
            var dic = new System.Collections.Generic.Dictionary<string, MItem>(160);
            AddUnits(ref dic, unitless,
                     "DELTA", "Delta",                                             //           IO Reduced density (rho/rhoc)
                     "Q",                                                          // [mol/mol] IO Mass vapor quality
                     "TAU", "Tau",                                                 //           IO Reciprocal reduced temperature (Tc/T)
                     "ACENTRIC", "acentric",                                       //           O  Acentric factor
                     "ALPHA0", "alpha0",                                           //           O  Ideal Helmholtz energy
                     "ALPHAR", "alphar",                                           //           O  Residual Helmholtz energy
                     "BVIRIAL", "Bvirial",                                         //           O  Second virial coefficient
                     "CVIRIAL", "Cvirial",                                         //           O  Third virial coefficient
                     "DALPHA0_DDELTA_CONSTTAU", "dalpha0_ddelta_consttau",         //           O  Derivative of ideal Helmholtz energy with delta
                     "DALPHA0_DTAU_CONSTDELTA", "dalpha0_dtau_constdelta",         //           O  Derivative of ideal Helmholtz energy with tau
                     "DALPHAR_DDELTA_CONSTTAU", "dalphar_ddelta_consttau",         //           O  Derivative of residual Helmholtz energy with delta
                     "DALPHAR_DTAU_CONSTDELTA", "dalphar_dtau_constdelta",         //           O  Derivative of residual Helmholtz energy with tau
                     "DBVIRIAL_DT", "dBvirial_dT",                                 //           O  Derivative of second virial coefficient with respect to T
                     "DCVIRIAL_DT", "dCvirial_dT",                                 //           O  Derivative of third virial coefficient with respect to T
                     "FH",                                                         //           O  Flammability hazard
                     "FRACTION_MAX", "fraction_max",                               //           O  Fraction (mole, mass, volume) maximum value for incompressible solutions
                     "FRACTION_MIN", "fraction_min",                               //           O  Fraction (mole, mass, volume) minimum value for incompressible solutions
                     "FUNDAMENTAL_DERIVATIVE_OF_GAS_DYNAMICS", "fundamental_derivative_of_gas_dynamics", // O Fundamental_derivative_of_gas_dynamics
                     "GWP100",                                                     //           O  100-year global warming potential
                     "GWP20",                                                      //           O  20-year global warming potential
                     "GWP500",                                                     //           O  500-year global warming potential
                     "HH",                                                         //           O  Health hazard
                     "ODP",                                                        //           O  Ozone depletion potential
                     "PHASE", "Phase",                                             //           O  Phase index as a float
                     "PH",                                                         //           O  Physical hazard
                     "PRANDTL", "Prandtl",                                         //           O  Prandtl number
                     "Z");                                                         //           O  Compressibility factor
            AddUnits(ref dic, mol_m3,
                     "DMOLAR", "Dmolar",                                           // [mol/m^3] IO Molar density
                     "RHOMOLAR_CRITICAL", "rhomolar_critical",                     // [mol/m^3] O Molar density at critical point
                     "RHOMOLAR_REDUCING", "rhomolar_reducing");                    // [mol/m^3] O Molar density at reducing point
            AddUnits(ref dic, kg_m3,
                     "D", "DMASS", "Dmass",                                        // [kg/m^3]  IO Mass density
                     "RHOCRIT", "RHOMASS_CRITICAL", "rhocrit", "rhomass_critical", // [kg/m^3]  O  Mass density at critical point
                     "RHOMASS_REDUCING", "rhomass_reducing");                      // [kg/m^3]  O  Mass density at reducing point
            AddUnits(ref dic, J_mol,
                     "HMOLAR", "Hmolar",                                           // [J/mol]   IO Molar specific enthalpy
                     "UMOLAR", "Umolar",                                           // [J/mol]   IO Molar specific internal energy
                     "GMOLAR", "Gmolar",                               // [J/mol/K] O  Molar gas constant
                     "HMOLAR_RESIDUAL", "Hmolar_residual",                               // [J/mol/K] O  Molar gas constant
                     "GMOLAR_RESIDUAL", "Gmolar_residual");                                          // [J/mol]   O  Molar specific Gibbs energy
            AddUnits(ref dic, J_kg,
                     "H", "HMASS", "Hmass",                                        // [J/kg]    IO Mass specific enthalpy
                     "U", "UMASS", "Umass",                                        // [J/kg]    IO Mass specific internal energy
                     "G", "GMASS", "Gmass");                                       // [J/kg]    O  Mass specific Gibbs energy
            AddUnits(ref dic, Pa,
                     "P",                                                          // [Pa]      IO Pressure
                     "PCRIT", "P_CRITICAL", "Pcrit", "p_critical", "pcrit",        // [Pa]      O  Pressure at the critical point
                     "PMAX", "P_MAX", "P_max", "pmax",                             // [Pa]      O  Maximum pressure limit
                     "PMIN", "P_MIN", "P_min", "pmin",                             // [Pa]      O  Minimum pressure limit
                     "PTRIPLE", "P_TRIPLE", "p_triple", "ptriple",                 // [Pa]      O  Pressure at the triple point (pure only)
                     "P_REDUCING", "p_reducing");                                  // [Pa]      O  Pressure at the reducing point
            AddUnits(ref dic, J_mol_K,
                     "SMOLAR", "Smolar",                                           // [J/mol/K] IO Molar specific entropy
                     "CP0MOLAR", "Cp0molar",                                       // [J/mol/K] O  Ideal gas molar specific constant presssure specific heat
                     "CPMOLAR", "Cpmolar",                                         // [J/mol/K] O  Molar specific constant presssure specific heat
                     "CVMOLAR", "Cvmolar",                                         // [J/mol/K] O  Molar specific constant volume specific heat
                     "GAS_CONSTANT", "gas_constant",                               // [J/mol/K] O  Molar gas constant
                     "SMOLAR_RESIDUAL", "Smolar_residual");                        // [J/mol/K] O  Residual molar entropy (sr/R = tau*dar_dtau-ar)
            AddUnits(ref dic, J_kg_K,
                     "S", "SMASS", "Smass",                                        // [J/kg/K]  IO Mass specific entropy
                     "CVMASS", "Cvmass", "O",                                      // [J/kg/K]  O  Mass specific constant volume specific heat
                     "CP0MASS", "Cp0mass",                                         // [J/kg/K]  O  Ideal gas mass specific constant presssure specific heat
                     "C", "CPMASS", "Cpmass");                                     // [J/kg/K]  O  Mass specific constant presssure specific heat
            AddUnits(ref dic, K,
                     "T",                                                          // [K]       IO Temperature
                     "TCRIT", "T_CRITICAL", "T_critical", "Tcrit",                 // [K]       O  Temperature at the critical point
                     "TMAX", "T_MAX", "T_max", "Tmax",                             // [K]       O  Maximum temperature limit
                     "TMIN", "T_MIN", "T_min", "Tmin",                             // [K]       O  Minimum temperature limit
                     "TTRIPLE", "T_TRIPLE", "T_triple", "Ttriple",                 // [K]       O  Temperature at the triple point
                     "T_FREEZE", "T_freeze",                                       // [K]       O  Freezing temperature for incompressible solutions
                     "T_REDUCING", "T_reducing");                                  // [K]       O  Temperature at the reducing point
            AddUnits(ref dic, m_s,
                     "A", "SPEED_OF_SOUND", "speed_of_sound");                     // [m/s]     O  Speed of sound
            AddUnits(ref dic, W_m_K,
                     "CONDUCTIVITY", "L", "conductivity");                         // [W/m/K]   O  Thermal conductivity
            AddUnits(ref dic, _K,
                     "ISOBARIC_EXPANSION_COEFFICIENT", "isobaric_expansion_coefficient"); // [1/K] O Isobaric expansion coefficient
            AddUnits(ref dic, _Pa,
                     "ISOTHERMAL_COMPRESSIBILITY", "isothermal_compressibility");  // [1/Pa]    O  Isothermal compressibility
            AddUnits(ref dic, N_m,
                     "I", "SURFACE_TENSION", "surface_tension");                   // [N/m]     O  Surface tension
            AddUnits(ref dic, kg_mol,
                     "M", "MOLARMASS", "MOLAR_MASS", "MOLEMASS", "molar_mass", "molarmass", "molemass"); // kg/mol O Molar mass
            AddUnits(ref dic, Pa_s,
                     "V", "VISCOSITY", "viscosity");                               // [Pa s]    O  Viscosity

            return dic;
        }
        static System.Collections.Generic.Dictionary<string, MItem> dic = InitUnitsDictionary();

        public static MItem Find(string param)
        {
            if (!dic.ContainsKey(param))
                return unitless;

            // HACK: fix for SS-2414
            return SMath.Math.Decision.SymbolicCalculation(dic[param].ToEntry(), new SMath.Math.Store());
            //return dic[param];
        }

        static System.Collections.Generic.Dictionary<string, MItem> InitHAUnitsDictionary()
        {
            var dic = new System.Collections.Generic.Dictionary<string, MItem>(50);
            AddUnits(ref dic, unitless,
                     "Omega", "HumRat", "W",         // Humidity Ratio [kg water/kg dry air]
                     "RH", "RelHum", "R",            // Relative humidity in (0,1) [-]
                     "psi_w", "Y",                   // Mole fraction of water [mol_w/mol]
                     "Z");                           // Compressibility factor (Z=pv/(RT)) [-]
            AddUnits(ref dic, K,
                     "Tdb", "T_db", "T",             // Dry-Bulb Temperature [K]
                     "Tdp", "T_dp", "DewPoint", "D", // Dew-Point Temperature [K]
                     "Twb", "T_wb", "WetBulb", "B"); // Wet-Bulb Temperature [K]
            AddUnits(ref dic, J_kg,
                     "Enthalpy", "H", "Hda",         // Mixture enthalpy per dry air [J/kg dry air]
                     "Hha");                         // Mixture enthalpy per humid air [J/kg humid air]
            AddUnits(ref dic, J_kg_K,
                     "Entropy", "S", "Sda",          // Mixture entropy per unit dry air [J/kg dry air/K]
                     "Sha",                          // Mixture entropy per unit humid air [J/kg humid air/K]
                     "C", "cp",                      // Mixture specific heat per unit dry air [J/kg dry air/K]
                     "Cha", "cp_ha");                // Mixture specific heat per unit humid air [J/kg humid air/K]
            AddUnits(ref dic, Pa,
                     "P",                            // Pressure [Pa]
                     "P_w");                         // Partial pressure of water vapor [Pa]
            AddUnits(ref dic, m3_kg,
                     "V", "Vda",                     // Mixture volume per unit dry air [m^3/kg dry air]
                     "Vha");                         // Mixture volume per unit humid air [m^3/kg humid air]
            AddUnits(ref dic, Pa_s,
                     "mu", "Visc", "M");             // Mixture viscosity [Pa-s]
            AddUnits(ref dic, W_m_K,
                     "k", "Conductivity", "K");      // Mixture thermal conductivity [W/m/K]

            return dic;
        }
        static System.Collections.Generic.Dictionary<string, MItem> HAdic = InitHAUnitsDictionary();

        public static MItem FindHA(string param)
        {
            if (!HAdic.ContainsKey(param))
                return unitless;

            // HACK: fix for SS-2414
            return SMath.Math.Decision.SymbolicCalculation(HAdic[param].ToEntry(), new SMath.Math.Store());
            //return HAdic[param];
        }

        static void matchInternal(MItem unit, SMath.Math.Numeric.TNumber val, ref SMath.Math.Store context)
        {
            // First, reduce our units to basic representation to allow MathEqual()
            var t = SMath.Math.Decision.Preprocessing(unit.ToEntry(), context);
            var u = Expression.SimplifyEx(Converter.ToMItem(t), context);
            if (!val.obj.Units.MathEqual(u) && !val.obj.Units.MathEqual(unitless))
                throw new SMath.Manager.MathException(SMath.Manager.Errors.UnitsDontMatch);
        }
        public static void match(string param, SMath.Math.Numeric.TNumber val, SMath.Math.Store context)
        {
            matchInternal(Find(param), val, ref context);
        }
        public static void matchHA(string param, SMath.Math.Numeric.TNumber val, SMath.Math.Store context)
        {
            matchInternal(FindHA(param), val, ref context);
        }
    }
}
