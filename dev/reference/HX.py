"""
Supplemental code for paper:
I. Bell et al., "A Generalized Moving-Boundary Algorithm to Predict the Heat Transfer Rate of
Counterflow Heat Exchangers for any Phase Configuration", Applied Thermal Engineering, 2014

================================ PORTING REFERENCE ============================
This is the ORIGINAL author's supplemental script for the moving-boundary HX
paper, recovered verbatim from the 2014 archive and kept here as the reference
oracle for the C++ port in `dev/demo_hx_speedup.cpp` (issue CoolProp-q00v).

Recovered from:
    /Volumes/Windows7_OS/Users/Belli/Documents/Papers/Completed/1_ 2014 BPHE/Figs/HX.py
    (mtime 2014-10-14; the cleaner of two copies — equation numbers wired into
     the comments rather than the earlier draft's "Equation XXX" placeholders.)

The NUMERICS are intentionally left untouched so this file can serve as a
bit-faithful oracle for the port. Only compatibility/correctness edits needed
to make it import and run on a modern toolchain were applied; every such change
is tagged inline with `# [PORT]` and summarized here:

  1. CoolProp.Plots import (`from CoolProp.Plots import Ph, Ts`) — the old
     functional `Ph`/`Ts` plot API was REMOVED in CoolProp >= 6.x (verified
     absent in 7.2.1dev). Because it was a top-level import, the whole module
     failed to import on a current install, so even the non-plotting solver
     path was unreachable. Made the import lazy + optional; `plot_ph_pair` /
     `plot_Ts_pair` now raise a clear message if the legacy API is missing.
     The solver, `calculate_cell_boundaries`, and `plot_cells` are unaffected.
  2. `plt.xlabel('$\\hat h$ [-]')` — the bare backslash is an invalid escape
     under Python 3.12+ (SyntaxWarning -> future SyntaxError). Made it a raw
     string. No visual change.

Behavioural caveats for the port (NOT changed here — preserve & match them):

  * Table-3 mass-flow TRANSPOSITION (real, in the original): in
    `PropaneEvaporatorPinching()` the HeatExchanger is constructed as
    `HeatExchanger('Water', mdot_c, ..., 'n-Propane', mdot_h, ...)`. With
    mdot_h=0.01 and mdot_c=0.1 *as defined locally*, the HOT fluid (Water) is
    actually given 0.1 kg/s and the COLD fluid (n-Propane) 0.01 kg/s — i.e. the
    opposite of Table 3 in the paper. This is the "apparent mass-flow
    transposition" the design spec flags; it is what produces the evaporating
    cell structure of Figs 8/9. The C++ port must reproduce this assignment to
    match, then document it.
  * Two-phase heat-transfer coefficient is hard-coded to alpha = 1000 W/m^2/K
    (liquid/vapor = 100) inside `objective_function`. The design spec lists the
    two-phase alpha as 2000. Reconcile this before asserting the correctness
    gate — the oracle here uses 1000.
  * Brent bracket is [1e-5, Qmax - 1e-10] with rtol=1e-14, xtol=1e-10.
===============================================================================
"""

from __future__ import division, print_function

import CoolProp
import CoolProp.CoolProp as CP

# [PORT] Legacy functional plot API (Ph/Ts) was removed in modern CoolProp.
# Import lazily so the solver path works even when it is unavailable; the two
# methods that use it degrade with a clear error instead of breaking import.
try:
    from CoolProp.Plots import Ph, Ts  # noqa: F401  (legacy API, may be absent)
    _HAVE_LEGACY_PLOTS = True
except ImportError:
    _HAVE_LEGACY_PLOTS = False

import matplotlib.pyplot as plt
import numpy as np
from math import log
import scipy.optimize

# Set to True to enable some debugging output to screen
debug = False

class struct(object): 
    """
    
    A dummy class that allows you to set variables like::
    
        S = struct()
        S.A = 'apple'
        S.N = 3
    """
    pass
    
class HeatExchanger(object):
    
    def __init__(self, Fluid_h, mdot_h, p_hi, h_hi, Fluid_c, mdot_c, p_ci, h_ci):
        """
        
        Parameters
        ----------
        
        """
        
        # Set variables in the class instance
        self.Fluid_h = Fluid_h
        self.mdot_h = mdot_h
        self.h_hi = h_hi
        self.p_hi = p_hi
        self.Fluid_c = Fluid_c
        self.mdot_c = mdot_c
        self.h_ci = h_ci
        self.p_ci = p_ci
        
        # Determine the inlet temperatures from the pressure/enthalpy pairs
        self.T_ci = CP.PropsSI('T', 'P', self.p_ci, 'H', self.h_ci, self.Fluid_c)
        self.T_hi = CP.PropsSI('T', 'P', self.p_hi, 'H', self.h_hi, self.Fluid_h)
        
        # Calculate the bubble and dew enthalpies for each stream
        self.T_cbubble = CP.PropsSI('T', 'P', self.p_ci, 'Q', 0, self.Fluid_c)
        self.T_cdew    = CP.PropsSI('T', 'P', self.p_ci, 'Q', 1, self.Fluid_c)
        self.T_hbubble = CP.PropsSI('T', 'P', self.p_hi, 'Q', 0, self.Fluid_h)
        self.T_hdew    = CP.PropsSI('T', 'P', self.p_hi, 'Q', 1, self.Fluid_h)
        self.h_cbubble = CP.PropsSI('H', 'T', self.T_cbubble, 'Q', 0, self.Fluid_c)
        self.h_cdew    = CP.PropsSI('H', 'T', self.T_cdew, 'Q', 1, self.Fluid_c)
        self.h_hbubble = CP.PropsSI('H', 'T', self.T_hbubble, 'Q', 0, self.Fluid_h)
        self.h_hdew    = CP.PropsSI('H', 'T', self.T_hdew, 'Q', 1, self.Fluid_h)
        
    
    def external_pinching(self):
        """ Determine the maximum heat transfer rate based on the external pinching analysis """

        # Equation 5
        self.h_ho = CP.PropsSI('H','T',self.T_ci,'P',self.p_hi,self.Fluid_h)        

        # Equation 4
        Qmaxh = self.mdot_h*(self.h_hi-self.h_ho)
        
        # Equation 7
        self.h_co = CP.PropsSI('H','T',self.T_hi,'P',self.p_ci,self.Fluid_c)
        
        # Equation 6
        Qmaxc = self.mdot_c*(self.h_co-self.h_ci)
        
        Qmax = min(Qmaxh, Qmaxc)
        
        if debug:
            print('Qmax (external pinching) is', Qmax)
        
        self.calculate_cell_boundaries(Qmax)
        
        return Qmax
        
    def calculate_cell_boundaries(self, Q):
        """ Calculate the cell boundaries for each fluid """
        
        # Re-calculate the outlet enthalpies of each stream
        self.h_co = self.h_ci + Q/self.mdot_c
        self.h_ho = self.h_hi - Q/self.mdot_h
        
        # Start with the external boundaries (sorted in increasing enthalpy)
        self.hvec_c = [self.h_ci, self.h_co]
        self.hvec_h = [self.h_ho, self.h_hi]
        
        # Add the bubble and dew enthalpies for the hot stream
        if self.h_hdew is not None and self.h_hi > self.h_hdew > self.h_ho:
            self.hvec_h.insert(-1, self.h_hdew)
        if self.h_hbubble is not None and self.h_hi > self.h_hbubble > self.h_ho:
            self.hvec_h.insert(1, self.h_hbubble)
        
        # Add the bubble and dew enthalpies for the cold stream
        if self.h_cdew is not None and self.h_ci < self.h_cdew < self.h_co:
            self.hvec_c.insert(-1, self.h_cdew)
        if self.h_cbubble is not None and self.h_ci < self.h_cbubble < self.h_co:
            self.hvec_c.insert(1, self.h_cbubble)
            
        if debug:
            print(self.hvec_c, self.hvec_h)
            
        # Fill in the complementary cell boundaries
        # Start at the first element in the vector
        k = 0
        while k < len(self.hvec_c)-1 or k < len(self.hvec_h)-1:
            if len(self.hvec_c) == 2 and len(self.hvec_h) == 2:
                break
                
            # Determine which stream is the limiting next cell boundary
            Qcell_hk = self.mdot_h*(self.hvec_h[k+1]-self.hvec_h[k])
            Qcell_ck = self.mdot_c*(self.hvec_c[k+1]-self.hvec_c[k])
            
            if abs(Qcell_hk/Qcell_ck - 1)< 1e-6:
                k +=1
                break
            elif Qcell_hk > Qcell_ck:
                # Hot stream needs a complementary cell boundary
                self.hvec_h.insert(k+1, self.hvec_h[k] + Qcell_ck/self.mdot_h)
            else:
                # Cold stream needs a complementary cell boundary
                self.hvec_c.insert(k+1, self.hvec_c[k] + Qcell_hk/self.mdot_c)
            
            if debug:
                print(k,len(self.hvec_c),len(self.hvec_h),Qcell_hk, Qcell_ck)
            
            if debug:
                # Calculate the temperature and entropy at each cell boundary
                self.Tvec_c = CP.PropsSI('T','H',self.hvec_c,'P',self.p_ci,self.Fluid_c)
                self.Tvec_h = CP.PropsSI('T','H',self.hvec_h,'P',self.p_hi,self.Fluid_h)
                self.svec_c = CP.PropsSI('S','H',self.hvec_c,'P',self.p_ci,self.Fluid_c)
                self.svec_h = CP.PropsSI('S','H',self.hvec_h,'P',self.p_hi,self.Fluid_h)
                self.plot_cells()
                plt.show()
            
            Qcell_hk = self.mdot_h*(self.hvec_h[k+1]-self.hvec_h[k])
            Qcell_ck = self.mdot_c*(self.hvec_c[k+1]-self.hvec_c[k])
            assert (abs(Qcell_hk/Qcell_ck-1) < 1e-6)
            
            # Increment index
            k += 1
        
        assert(len(self.hvec_h) == len(self.hvec_c))
        Qhs = np.array([self.mdot_h*(self.hvec_h[i+1]-self.hvec_h[i]) for i in range(len(self.hvec_h)-1)])
        Qcs = np.array([self.mdot_c*(self.hvec_c[i+1]-self.hvec_c[i]) for i in range(len(self.hvec_c)-1)])
        if debug:
            if np.max(np.abs(Qcs/Qhs))<1e-5:
                print(Qhs, Qcs)
        
        # Calculate the temperature and entropy at each cell boundary
        self.Tvec_c = CP.PropsSI('T','H',self.hvec_c,'P',self.p_ci,self.Fluid_c)
        self.Tvec_h = CP.PropsSI('T','H',self.hvec_h,'P',self.p_hi,self.Fluid_h)
        self.svec_c = CP.PropsSI('S','H',self.hvec_c,'P',self.p_ci,self.Fluid_c)
        self.svec_h = CP.PropsSI('S','H',self.hvec_h,'P',self.p_hi,self.Fluid_h)
        
        # Calculate the phase in each cell
        self.phases_h = []
        for i in range(len(self.hvec_h)-1):
            havg = (self.hvec_h[i] + self.hvec_h[i+1])/2.0
            if havg < self.h_hbubble:
                self.phases_h.append('liquid')
            elif havg > self.h_hdew:
                self.phases_h.append('vapor')
            else:
                self.phases_h.append('two-phase')
        
        self.phases_c = []
        for i in range(len(self.hvec_c)- 1):
            havg = (self.hvec_c[i] + self.hvec_c[i+1])/2.0
            if havg < self.h_cbubble:
                self.phases_c.append('liquid')
            elif havg > self.h_cdew:
                self.phases_c.append('vapor')
            else:
                self.phases_c.append('two-phase')
            
    def internal_pinching(self, stream):
        """
        Determine the maximum heat transfer rate based on the internal pinching analysis 
        """
        
        if stream == 'hot':
            
            # Try to find the dew point enthalpy as one of the cell boundaries
            # that is not the inlet or outlet
        
            # Check for the hot stream pinch point
            for i in range(1,len(self.hvec_h)-1):
                
                # Check if enthalpy is equal to the dewpoint enthalpy of hot
                # stream and hot stream is colder than cold stream (impossible)
                if (abs(self.hvec_h[i] - self.h_hdew) < 1e-6 
                            and self.Tvec_c[i] > self.Tvec_h[i]):
            
                    # Enthalpy of the cold stream at the pinch temperature
                    # Equation 10
                    h_c_pinch = CP.PropsSI('H','T',self.T_hdew,'P',self.p_ci, self.Fluid_c)
                    
                    # Heat transfer in the cell
                    # Equation 9
                    Qright = self.mdot_h*(self.h_hi-self.h_hdew)
                    
                    # New value for the limiting heat transfer rate
                    # Equation 12
                    Qmax = self.mdot_c*(h_c_pinch-self.h_ci) + Qright
                    
                    # Recalculate the cell boundaries
                    self.calculate_cell_boundaries(Qmax)

                    return Qmax
        
        elif stream == 'cold':
            # Check for the cold stream pinch point
            for i in range(1,len(self.hvec_c)-1):
                
                # Check if enthalpy is equal to the bubblepoint enthalpy of cold
                # stream and hot stream is colder than cold stream (impossible)
                if (abs(self.hvec_c[i] - self.h_cbubble) < 1e-6 
                            and self.Tvec_c[i] > self.Tvec_h[i]):
            
                    # Enthalpy of the cold stream at the pinch temperature
                    # Equation 14
                    h_h_pinch = CP.PropsSI('H','T',self.T_cbubble,'P',self.p_hi, self.Fluid_h)
                    
                    # Heat transfer in the cell
                    # Equation 13
                    Qleft = self.mdot_c*(self.h_cbubble-self.h_ci)
                    
                    # New value for the limiting heat transfer rate
                    # Equation 16
                    Qmax = Qleft + self.mdot_h*(self.h_hi-h_h_pinch)
                    
                    # Recalculate the cell boundaries
                    self.calculate_cell_boundaries(Qmax)

                    return Qmax
        else:
            raise ValueError
        
    def run(self, only_external = False, and_solve = False):
        # Check the external pinching & update cell boundaries  
        Qmax_ext = self.external_pinching()
        Qmax = Qmax_ext
    
        if not only_external:
            # Check the internal pinching
            for stream in ['hot','cold']:
                # Check stream internal pinching & update cell boundaries
                Qmax_int = self.internal_pinching(stream)
                if Qmax_int is not None:
                    Qmax = Qmax_int
                
        self.Qmax = Qmax
        
        if and_solve and not only_external:
            Q = self.solve()
            
        Qtotal = self.mdot_c*(self.hvec_c[-1]-self.hvec_c[0])
        
        # Build the normalized enthalpy vectors
        self.hnorm_h = self.mdot_h*(np.array(self.hvec_h)-self.hvec_h[0])/Qtotal
        self.hnorm_c = self.mdot_c*(np.array(self.hvec_c)-self.hvec_c[0])/Qtotal
        
        if and_solve:
            return Q
        
    def objective_function(self, Q):
        
        self.calculate_cell_boundaries(Q)

        w = []
        for k in range(len(self.hvec_c)-1):
            Thi = self.Tvec_h[k+1]
            Tci = self.Tvec_c[k]
            Tho = self.Tvec_h[k]
            Tco = self.Tvec_c[k+1]
            DTA = Thi - Tco
            DTB = Tho - Tci

            if DTA == DTB:
                LMTD = DTA
            else:
                try:
                    LMTD = (DTA-DTB)/log(abs(DTA/DTB))
                except ValueError as VE:
                    print(Q, DTA, DTB)
                    raise
            UA_req = self.mdot_h*(self.hvec_h[k+1]-self.hvec_h[k])/LMTD
            if self.phases_c[k] in ['liquid','vapor']:
                alpha_c = 100
            else:
                alpha_c = 1000
            if self.phases_h[k] in ['liquid','vapor']:
                alpha_h = 100
            else:
                alpha_h = 1000
            
            UA_avail = 1/(1/(alpha_h*self.A_h)+1/(alpha_c*self.A_c))
            w.append(UA_req/UA_avail)
            
        if debug:
            print(Q, 1-sum(w))
            
        return 1-sum(w)
    
    def solve(self):
        """ 
        Solve the objective function using Brent's method and the maximum heat transfer 
        rate calculated from the pinching analysis
        """
        self.Q = scipy.optimize.brentq(self.objective_function, 1e-5, self.Qmax-1e-10, rtol = 1e-14, xtol = 1e-10)
        return self.Q
        
    def plot_objective_function(self, N = 100):
        """ Plot the objective function """
        Q = np.linspace(1e-5,self.Qmax,N)
        r = np.array([self.objective_function(_Q) for _Q in Q])
        plt.plot(Q, r)
        plt.show()
        
    def plot_ph_pair(self):
        """ Plot p-h plots for the pair of working fluids """
        if not _HAVE_LEGACY_PLOTS:  # [PORT] legacy Ph() removed in modern CoolProp
            raise RuntimeError("plot_ph_pair requires the legacy CoolProp.Plots.Ph API, "
                               "which is unavailable in this CoolProp version.")
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        Ph(self.Fluid_h, axis = ax1)
        Ph(self.Fluid_c, axis = ax2)
        ax1.set_title('')
        ax1.plot([self.h_hi, self.h_ho],[self.p_hi,self.p_hi],'s-')
        ax2.set_title('')
        ax2.plot([self.h_ci, self.h_co],[self.p_ci,self.p_ci],'s-')
        plt.show()
    
    def plot_Ts_pair(self):
        """ Plot a T-s plot for the pair of working fluids """
        if not _HAVE_LEGACY_PLOTS:  # [PORT] legacy Ts() removed in modern CoolProp
            raise RuntimeError("plot_Ts_pair requires the legacy CoolProp.Plots.Ts API, "
                               "which is unavailable in this CoolProp version.")
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        Ts(self.Fluid_h, axis = ax1)
        Ts(self.Fluid_c, axis = ax2)
        ax1.set_title('')
        ax1.plot(self.svec_h, self.Tvec_h, 's-')
        ax2.set_title('')
        ax2.plot(self.svec_c, self.Tvec_c, 's-')
        plt.show()
        
    def plot_cells(self, fName = '', dpi = 400):
        """ Plot the cells of the heat exchanger """
        plt.figure(figsize = (2.4,2.4))
        plt.plot(self.hnorm_h, self.Tvec_h, 'rs-')
        plt.plot(self.hnorm_c, self.Tvec_c, 'bs-')
        plt.xlim(0,1)
        plt.ylabel('T [K]') 
        plt.xlabel(r'$\hat h$ [-]')  # [PORT] raw string: bare \h is an invalid py3.12+ escape
        plt.tight_layout(pad = 0.2) 
        if fName != '':
            plt.savefig(fName, dpi = dpi)
        
def PropaneEvaporatorPinching():
    p_Water = 101325
    h_Water = CP.PropsSI('H','T',330,'P',p_Water,'Water')
    mdot_h = 0.01
    
    p_ref = CP.PropsSI('P','T',300,'Q',1,'n-Propane')
    h_ref = CP.PropsSI('H','T',275,'P',p_ref,'n-Propane')
    mdot_c = 0.1
    
    HX = HeatExchanger('Water',mdot_c,p_Water,h_Water,'n-Propane',mdot_h,p_ref,h_ref) 
    
    HX.A_h = HX.A_c = 4
    #Actually run the HX code
    HX.run(and_solve = True)
    HX.plot_cells('full.png')
    
if __name__=='__main__':
    # If the script is run directly, this code will be executed.
    PropaneEvaporatorPinching()
    