import pylab
import numpy as np
import scipy.stats

class FluidClass(object):
    pass

def is_float(a):
    try:
        a = float(a)
        return True
    except ValueError:
        return False
     
def parse_Melinder_nonaqueous():
    
    lines = open('Melinder coefficients.csv').readlines()
    
    #find #Fluid# dividers
    dividers = [i for i in range(len(lines)) if lines[i].find('#Fluid#') > -1]
    dividers.append(len(lines))
    Fluids = []
    for i in range(len(dividers)-1):
        Fluid = FluidClass()
        L = dividers[i]
        R = dividers[i+1]
        
        l = lines[L:R]
        l.pop(0) # #Fluid# remove
        Fluid.name, Fluid.longname = l.pop(0).split(',')[0].split('::')
        Fluid.description = l.pop(0).split(',')[0]
        l.pop(0) #empty line
        l.pop(0) #header
        l.pop(0) #units
        
        Fluid.t,Fluid.rho,Fluid.cp,Fluid.k,Fluid.mu = [],[],[],[],[]
        for line in l:
            ls = line.split(',')
            if is_float(ls[0]):
                tf,t,rho,cp,k,mu =  [float(thing) for thing in line.split(',')]
            else:
                if not is_float(ls[1]):
                    break
                t,rho,cp,k,mu =  [float(thing) for thing in line.split(',')[1::]]
            Fluid.t.append(t)
            Fluid.rho.append(rho)
            Fluid.cp.append(cp)
            Fluid.k.append(k)
            Fluid.mu.append(mu)
        pylab.plot(Fluid.t, Fluid.rho)
        
        Fluid.t = np.array(Fluid.t)
        Fluid.rho = np.array(Fluid.rho)
        Fluid.k = np.array(Fluid.k)
        Fluid.mu = np.array(Fluid.mu)
        Fluid.cp = np.array(Fluid.cp)
        
        Fluids.append(Fluid)
        
    pylab.close()
    
    

    #Create the classes
    for Fluid in Fluids:
        Fluid.ClassName = Fluid.name.strip()+"LiquidClass"
        #String for density
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Fluid.t+273.15,Fluid.rho)
        Fluid.rhostr = "{slope:g}*T_K+{intercept:g}".format(slope = slope, intercept = intercept)
        #String for conductivity
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Fluid.t+273.15,Fluid.k)
        Fluid.condstr = "{slope:g}*T_K+{intercept:g}".format(slope = slope/1000, intercept = intercept/1000)
        #String for specific heat
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Fluid.t+273.15,Fluid.cp)
        Fluid.cpstr = "{slope:g}*T_K+{intercept:g}".format(slope = slope/1000, intercept = intercept/1000)
        #String for internal energy
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Fluid.t+273.15,Fluid.cp)
        Fluid.ustr = "{slope:g}*(T_K*T_K-298*298)/2.0+{intercept:g}*(T_K-298)".format(slope = slope/1000, intercept = intercept/1000)
        #String for entropy
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Fluid.t+273.15,Fluid.cp)
        Fluid.sstr = "{slope:g}*(T_K-298)/2.0+{intercept:g}*log(T_K/298)".format(slope = slope/1000, intercept = intercept/1000)
        #String for viscosity
        n = np.polyfit(Fluid.t+273.15,np.log(np.array(Fluid.mu)/1000),2)
        Fluid.viscstr = "exp({n0:g}*T_K*T_K+{n1:g}*T_K+{n2:g})".format(n0 = n[0], n1 = n[1], n2 = n[2])
    
        Fluid.description = Fluid.description.lstrip('"')
        import textwrap
        print textwrap.dedent(
        """
        class {ClassName:s} : public IncompressibleLiquid{{
            
        public:
            
            // Constructor
            {ClassName:s}(){{
                name = std::string("{name:s}");
                description = std::string("{description:s}");
                reference = std::string("{reference:s}");
            }};

            ///< Destructor.  No implementation
            ~{ClassName:s}(){{}};

            double rho(double T_K){{
                return {rhostr:s};
            }}
            double cp(double T_K){{
                return {cpstr:s};
            }}
            double u(double T_K){{
                return {ustr:s};
            }}
            double s(double T_K){{
                return {sstr:s};
            }}
            double visc(double T_K){{
                return {viscstr:s};
            }}
            double cond(double T_K){{
                return {condstr:s};
            }}
        }};
        """.format(ClassName = Fluid.ClassName,
                    rhostr = Fluid.rhostr,
                    cpstr = Fluid.cpstr,
                    ustr = Fluid.ustr,
                    viscstr = Fluid.viscstr,
                    sstr = Fluid.sstr,
                    condstr = Fluid.condstr,
                    name = Fluid.name.strip(),
                    description = Fluid.longname.strip(),
                    reference = Fluid.description+ ' - from Ake Melinder, 2010, \\"Properties of Secondary Working Fluids for Indirect Systems\\", IIR')
        )

parse_Melinder_nonaqueous()

#for i in range(len(lines)):
#    if lines[i].strip() == '#Fluid#':
#        print lines[i],lines[i+1],lines[i+2]