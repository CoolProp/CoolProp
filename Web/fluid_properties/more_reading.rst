*************************
Advanced Fluid Properties
*************************

Use of Extended Corresponding States for Transport Properties
-------------------------------------------------------------

For a limited selection of fluids, correlations are provided for the viscosity and the thermal conductivity.  But for many fluids, no correlations are available, and therefore other methods must be employed.  The extended corresponding states is a method of estimating the transport properties of a fluid by analogy with the transport properties of a fluid that are well defined.

Implementing the ECS method is quite a challenge, but CoolProp is one of the only fluid property databases that properly implements it.  And the onlyopen-source package that does.  A multi-step method is required, which is hopefully clearly laid out here.

To begin with, the reference fluid must be selected that the fluid of interest will be compared with.  Ideally the shape of the molecules should be similar, but in practice, most fluids use R134a as the reference fluid since its thermodynamic and transport properties are well quantified with reference-quality correlations.

Once the reference fluid has been selected, the conformal state of the reference fluid must be determined.  The conformal state is the state at which the transport properties of the reference fluid and the fluid of interest are (in theory) the same.  In practice, at low densities the shape factors are assumed to be unity, and the conformal temperature and molar density are obtained from 

.. math::

    T_0 = T\frac{T_0^{c}}{T_j^c}
    
.. math::

    \overline{\rho_0} = \overline{\rho}\frac{\overline{\rho_0}^c}{\overline{\rho_j}^c}

Exact solution for the conformal temperature

If you have Helmholtz EOS for both the fluid and the reference fluid, you need to find a conformal temperature for the reference fluid that will yield the same compressibility factor and the residual Helmholtz energy

.. math::

    Z_j(T_j,\rho_j) = Z_0(T_0,\rho_0)

.. math::

    \alpha_j^r(T_j,\rho_j) = \alpha_0^r(T_0,\rho_0)

where "j" is for the fluid of interest, and the subscript "0" is for the reference fluid.  The left side of each equation is already known from the equation of state.  Literature suggests that solving for :math:`T_0` and :math:`\rho_0` directly is quite challenging.  See McLinden 2000 or Klein 1997.

Alternatively, if the shape factors :math:`\theta` and :math:`\phi` are known, either from correlation or otherwise, the conformal temperature and density can be calculated directly.

.. math::

    T_0 = \frac{T}{f} = T\frac{T_0^{c}}{T_j^c\theta(T_j,\rho_j)}
    
.. math::

    \rho_0 = \rho h = \rho\frac{\rho_0^c}{\rho_j^c}\phi(T_j,\rho_j)


Conversion from ideal gas term to Helmholtz energy term
-------------------------------------------------------

Much of the time the coefficients for the ideal-gas part of the Helmholtz energy are given directly, but sometimes only the gas-specific heat is provided.  Therefore you need to be able to go from specific heat to ideal-gas Helmholtz Energy.  The ideal-gas Helmholtz energy is given by Equation 23 from Lemmon, 2004, Equations of State for Mixtures of R-32, R-125, R-134a, R-143a, and R-152a, J. Phys. Chem. Ref. Data, Vol. 33, No. 2, 2004 or

.. math::

    a_0 = -RT+RT\ln\frac{\rho T}{\rho_0T_0}+h_0^0-Ts_0^0+\int_{T_0}^T c_p^0(T)dT-T\int_{T_0}^T \frac{c_p^0(T)}{T}dT
    
non-dimensionalizing

.. math::

    \alpha_0 =\frac{a_0}{RT}= -1+\ln\frac{\rho T}{\rho_0T_0}+\frac{h_0^0}{RT}-\frac{s_0^0}{R}+\frac{1}{RT}\int_{T_0}^T c_p^0(T)dT-\frac{1}{R}\int_{T_0}^T \frac{c_p^0(T)}{T}dT
    
Now we might want to do a change of variable in the integrals.  If so, do a u-substitution in the integrals.
    
.. math::

    T=\frac{T_c}{\tau}

where

.. math::

    dT=-\frac{T_c}{\tau^2}d\tau
    
.. math::

    \alpha_0 = -1+\ln\frac{\rho T}{\rho_0T_0}+\frac{h_0^0}{RT}-\frac{s_0^0}{R}+\frac{1}{RT}\int_{\tau_0}^{\tau} c_p^0(T)(-\frac{T_c}{\tau^2}d\tau)-\frac{1}{R}\int_{\tau_0}^{\tau} \frac{c_p^0(\tau)}{T}(-\frac{T_c}{\tau^2}d\tau)
    
Simplifying and factoring the :math:`\tau` term yields

.. math::

    \alpha_0 = -1+\ln\frac{\rho T}{\rho_0T_0}+\frac{h_0^0}{RT}-\frac{s_0^0}{R}-\frac{\tau}{R}\int_{\tau_0}^{\tau} \frac{c_p^0(\tau)}{\tau^2}d\tau+\frac{1}{R}\int_{\tau_0}^{\tau} \frac{c_p^0(\tau)}{\tau}d\tau
        
which finally yields the solution as of Equation 3 from Lemmon, 2003 (and others)

The specific-heat contribution can then be taken as a sum of the contributions 

for a term of the form

.. math::

    \frac{c_p^0}{R}=\frac{(B/T)^2\exp(B/T)}{(\exp(B/T)-1)^2}

the contribution is found from 

.. math::

    \frac{1}{T}\int_{T_0}^T \frac{(B/T)^2\exp(B/T)}{(\exp(B/T)-1)^2} dT-\int_{T_0}^T \frac{(B/T)^2\exp(B/T)}{(\exp(B/T)-1)^2}\frac{1}{T}dT
    
.. math::

    \frac{1}{T} \left[ \frac{B}{\exp(B/T)-1 }\right|_{T_0}^T - \left[ \frac{B}{T}\left(\frac{1}{\exp(B/T)-1}+1\right) - \log[\exp(B/T)-1] \right|_{T_0}^T dT

Factor out a B, First two terms cancel, leaving

.. math::

    - \left[ \frac{B}{T} - \log[\exp(B/T)-1] \right|_{T_0}^T dT
    
.. math::

    \left[\log[\exp(B/T)-1] - \frac{B}{T} \right|_{T_0}^T dT
    
.. math::

    \log[\exp(B/T)-1] - \frac{B}{T} -(\log[\exp(B/T_0)-1] - \frac{B}{T_0})
    
or in terms of :math:`\tau`

.. math::

    \log[\exp(B\tau/Tc)-1] - \frac{B\tau}{Tc} -(\log[\exp(B\tau_0/T_c)-1] - \frac{B\tau_0}{T_c})
    
for a term of the form

.. math::

    \frac{c_p^0}{R}=c

the contribution is found from 

.. math::

    \frac{1}{T}\int_{T_0}^T c dT-\int_{T_0}^T \frac{c}{T}dT
    
.. math::

    \frac{c}{T}(T-T_0)-c\log(T/T_0)
    
or in terms of :math:`\tau`

.. math::

    c-\frac{cT_0\tau}{T_c}+c\log(\tau/\tau_0)
    
    
for a term of the form

.. math::

    \frac{c_p^0}{R}=cT^t, t \neq 0

the contribution is found from 

.. math::

    \frac{1}{T}\int_{T_0}^T c T^t dT-\int_{T_0}^T \frac{c T^t}{T}dT
    
.. math::

    \frac{c}{T}\left(\frac{T^{t+1}}{t+1}-\frac{T_0^{t+1}}{t+1}\right)-c\left(\frac{T^{t}}{t}-\frac{T_0^{t}}{t}\right)

.. math::

    cT^{t}\left(\frac{1}{t+1}-\frac{1}{t}\right)-c\frac{T_0^{t+1}}{T(t+1)}+c\frac{T_0^t}{t}

or in terms of :math:`\tau`

.. math::

    cT_c^{t}\tau^{-t}\left(\frac{1}{t+1}-\frac{1}{t}\right)-c\frac{T_0^{t+1}\tau}{T_c(t+1)}+c\frac{T_0^t}{t}
    
..
    .. math::
        
        \int\limits_{{\tau _0}}^\tau  {\left[ {aT_c^t{\tau ^{ - t - 1}}} \right]d\tau }  - \tau \int\limits_{{\tau _0}}^\tau  {\left[ {aT_c^t{\tau ^{ - t - 2}}} \right]d\tau } \\

    .. math::

        aT_c^t\left( {\int\limits_{{\tau _0}}^\tau  {{\tau ^{ - t - 1}}d\tau }  - \tau \int\limits_{{\tau _0}}^\tau  {{\tau ^{ - t - 2}}d\tau } } \right)\\

    if :math:`t=0`

    .. math::

        a\left( {\int\limits_{{\tau _0}}^\tau  {\frac{1}{\tau }d\tau }  - \tau \int\limits_{{\tau _0}}^\tau  {{\tau ^{ - 2}}d\tau } } \right)

    .. math::

        a\left( {\left[ {\ln \left( \tau  \right)} \right]_{{\tau _0}}^\tau  - \tau \left[ {\frac{{{\tau ^{ - 1}}}}{{ - 1}}} \right]_{{\tau _0}}^\tau } \right)
        
    .. math::
        a\left( \ln \left( \tau  \right) - \ln \left( {{\tau _0}} \right) \right)


    if :math:`t\neq0`:

    .. math::
        
        aT_c^t\left( {\left[ {\frac{{{\tau ^{ - t}}}}{{ - t}}} \right]_{{\tau _0}}^\tau  - \tau \left[ {\frac{{{\tau ^{ - t - 1}}}}{{ - t - 1}}} \right]_{{\tau _0}}^\tau } \right)\\

    .. math::

        aT_c^t\left( {\frac{{{\tau ^{ - t}}}}{{ - t}} - \frac{{\tau _0^{ - t}}}{{ - t}} - \tau \left[ {\frac{{{\tau ^{ - t - 1}}}}{{ - t - 1}} - \frac{{\tau _0^{ - t - 1}}}{{ - t - 1}}} \right]} \right)\\
     
    .. math::
     
        - aT_c^t\left( {\frac{{{\tau ^{ - t}}}}{t} - \frac{{\tau _0^{ - t}}}{t} - \left[ {\frac{{{\tau ^{ - t}}}}{{t + 1}} - \frac{{\tau _0^{ - t}}}{{t + 1}}} \right]} \right)
        
    .. math::

        - aT_c^t\frac{{{\tau ^{ - t}}}}{t} + aT_c^t\frac{{{\tau ^{ - t}}}}{{t + 1}} + aT_c^t\frac{{\tau _0^{ - t}}}{t} - aT_c^t\frac{{\tau _0^{ - t}}}{{t + 1}}

    .. math::

        - aT_c^t\frac{{{\tau ^{ - t}}}}{t} + aT_c^t\frac{{{\tau ^{ - t}}}}{{t + 1}} + aT_c^t\tau _0^{ - t}\left[ {\frac{1}{t} - \frac{1}{{t + 1}}} \right]
        
    .. math::

        aT_c^t{\tau ^{ - t}}\left[ {\frac{1}{{t + 1}} - \frac{1}{t}} \right] + aT_c^t\tau _0^{ - t}\left[ {\frac{1}{t} - \frac{1}{{t + 1}}} \right]\\
        
    if :math:`t = 1`

    .. math::
     
        - \frac{{a{T_c}{\tau ^{ - 1}}}}{2} + \frac{{a{T_c}\tau _0^{ - 1}}}{2}
    
These terms can be summarized by the following table:

.. math::

    \begin{array}{*{20}{c}}
    {\dfrac{{c_p^0}}{R}{\rm{ Term}}}&{{\alpha ^0}{\rm{ Term}}}&{{\rm{Class Name}}}&{}&{}&{}&{}&{}\\
    {{a_k}\dfrac{{{{\left( {{b_k}/T} \right)}^2}\exp \left( {{b_k}/T} \right)}}{{{{\left( {\exp \left( {{b_k}/T} \right) - 1} \right)}^2}}}}&{{a_k}\ln \left[ {1 - \exp \left( {\frac{{ - {b_k}\tau }}{{{T_c}}}} \right)} \right]}&{{\rm{phi0\_Planck\_Einstein}}(a,b/Tc,[iStart,iEnd])}&{}&{}&{}&{}&{}\\
    {ac\frac{{{{\left( {b/T} \right)}^2}\exp \left( { - b/T} \right)}}{{{{\left( {c\exp \left( { - b/T} \right) + 1} \right)}^2}}}}&{a\ln \left[ {c + \exp \left( {\frac{{b\tau }}{{{T_c}}}} \right)} \right]}&{{\rm{phi0\_Planck\_Einstein2}}(a,b/Tc,c)}&{}&{}&{}&{}&{}\\
    {yuck}&{{a_k}{\tau ^{{b_k}}}}&{{\rm{phi0\_power}}\left( {a,b,[iStart,iEnd]} \right)}&{}&{}&{}&{}&{}\\
    a&{a - a\frac{\tau }{{{\tau _0}}} + a\ln \left( {\frac{\tau }{{{\tau _0}}}} \right)}&{{\rm{phi0\_cp0\_constant}}(a,Tc,T0)}&{}&{}&{}&{}&{}\\
    {{a_1} + {a_2}{{\left( {\frac{{{a_3}/T}}{{\sinh \left( {{a_3}/T} \right)}}} \right)}^2} + {a_4}{{\left( {\frac{{{a_5}/T}}{{\cosh \left( {{a_5}/T} \right)}}} \right)}^2}}&{yuck}&{{\rm{phi0\_cp0\_AlyLee}}(a,Tc,T0,R)}&{}&{}&{}&{}&{}\\
    {{\rm{n/a}}}&{\log (\delta ) + {a_1} + {a_2}\tau }&{{\rm{phi0\_lead(}}a1,{\rm{ }}a2{\rm{)}}}&{}&{}&{}&{}&{}\\
    {{\rm{n/a}}}&{a\log \tau }&{{\rm{phi0\_logtau}}(a)}&{}&{}&{}&{}&{}
    \end{array}

If the reference enthalpy is known, you can determine the constants from 

.. math::

    \frac{h_0}{RT}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right] +\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}+1
    
.. math::

    \left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta} = \frac{1}{\tau}\left(\frac{h_0}{RT}-\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}-1\right)- \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta}
    
For the specific heat
The two integral terms are

.. math::
    
    - \frac{\tau }{R}\int_{{\tau _0}}^\tau  {\frac{{c_p^0}}{{{\tau ^2}}}d\tau }  + \frac{1}{R}\int_{{\tau _0}}^\tau  {\frac{{c_p^0}}{\tau }d\tau }

First derivative

.. math::

    \frac{d}{{d\tau }}\left[ { - \frac{\tau }{R}\int_{{\tau _0}}^\tau  {\frac{{c_p^0}}{{{\tau ^2}}}d\tau }  + \frac{1}{R}\int_{{\tau _0}}^\tau  {\frac{{c_p^0}}{\tau }d\tau } } \right] =  - \frac{{c_p^0}}{{\tau R}} - \frac{1}{R}\int_{{\tau _0}}^\tau  {\frac{{c_p^0}}{{{\tau ^2}}}d\tau }  + \frac{{c_p^0}}{{\tau R}} =  - \frac{1}{R}\int_{{\tau _0}}^\tau  {\frac{{c_p^0}}{{{\tau ^2}}}d\tau }

Second Derivative

.. math::

    \frac{{{d^2}}}{{d{\tau ^2}}}\left[ { - \frac{\tau }{R}\int_{{\tau _0}}^\tau  {\frac{{c_p^0}}{{{\tau ^2}}}d\tau }  + \frac{1}{R}\int_{{\tau _0}}^\tau  {\frac{{c_p^0}}{\tau }d\tau } } \right] = \frac{d}{{d\tau }}\left[ { - \frac{1}{R}\int_{{\tau _0}}^\tau  {\frac{{c_p^0}}{{{\tau ^2}}}d\tau } } \right] =  - \frac{{c_p^0}}{{{\tau ^2}R}}


Converting Bender and mBWR EOS
------------------------------

If the EOS is of the form

.. math::

    \frac{p}{{\rho RT}} = Z\left( {T,\rho } \right) = 1 + \sum\limits_i {{n_i}{T^{{s_i}}}{\rho ^{{r_i}}}}  + \sum\limits_i {{n_i}{T^{{s_i}}}{\rho ^{{r_i}}}} \exp \left( { - {\gamma _i}{{\left( {\frac{\rho }{{{\rho _c}}}} \right)}^2}} \right)

To convert to standard power form in CoolProp, use

.. math::

    \delta \sum\limits_i {{d_i}{a_i}{\tau ^{{t_i}}}{\delta ^{{d_i} - 1}}}  = \sum\limits_i {{n_i}{T^{{s_i}}}{\rho ^{{r_i}}}}  = \sum\limits_i {{n_i}{{\left( {\frac{{{T_c}}}{\tau }} \right)}^{{s_i}}}{{\left( {{\rho _c}\delta } \right)}^{{r_i}}}}  = \sum\limits_i {{n_i}T_c^{{s_i}}\rho _c^{{r_i}}{\tau ^{ - {s_i}}}{\delta ^{{r_i}}}}

The left-hand-side is the derivative of the residual Helmholtz energy with respect 
to delta times the reduced density since 

.. math::

    \frac{p}{\rho RT}=1+\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}

where

.. math::

    \delta : {d_i} - 1 + 1 = {r_i} \Rightarrow {d_i} = {r_i}
    
.. math::

    \tau : {t_i} =  - {s_i}

.. math::

    c : {d_i}{a_i} = {n_i}T_c^{{s_i}}\rho _c^{{r_i}}

.. math::

    p = \rho RT + \sum\limits_i {{n_i}{T^{{s_i}}}{\rho ^{{r_i}}}}  + \sum\limits_i {{n_i}{T^{{s_i}}}{\rho ^{{r_i}}}} \exp \left( { - {\gamma _i}{{\left( {\frac{\rho }{{{\rho _c}}}} \right)}^2}} \right){\rm{   (Eq 3}}{\rm{.28)}}
    
.. math::

    \frac{p}{{\rho RT}} = 1 + \sum\limits_i {\frac{{{n_i}}}{R}{T^{{s_i} - 1}}{\rho ^{{r_i} - 1}}}  + \sum\limits_i {\frac{{{n_i}}}{R}{T^{{s_i} - 1}}{\rho ^{{r_i} - 1}}} \exp \left( { - {\gamma _i}{{\left( {\frac{\rho }{{{\rho _c}}}} \right)}^2}} \right)
    
.. math::
    
    \delta \sum\limits_i {{d_i}{a_i}{\tau ^{{t_i}}}{\delta ^{{d_i} - 1}}}  = \sum\limits_i {\frac{{{n_i}}}{R}{{\left( {\frac{{{T_c}}}{\tau }} \right)}^{{s_i} - 1}}{{\left( {{\rho _c}\delta } \right)}^{{r_i} - 1}}}  = \sum\limits_i {\frac{{{n_i}}}{R}T_c^{{s_i} - 1}\rho _c^{{r_i} - 1}{\tau ^{ - ({s_i} - 1)}}{\delta ^{{r_i} - 1}}} 

.. math::

    \delta :1 + {d_i} - 1 = {r_i} - 1

.. math::

    \tau :{t_i} =  - \left( {s_i - 1} \right)
    
.. math::
    
    c:{d_i}{a_i} = \frac{{{n_i}}}{R}T_c^{{s_i} - 1}\rho _c^{{r_i} - 1}
    
In the Bender EOS, for the exponential part you have terms that can be converted to reduced form

.. math::
    
    a_i\delta^{d_i}\tau^{t_i}\exp(-\gamma \delta^2)
    
which yields the terms in the following table (from Span, 2000)

.. math::

    \begin{array}{*{4}{c}|*{4}{c}}
    \multicolumn{4}{c}{\mbox{From Bender}} & \multicolumn{4}{c}{\mbox{Power term}}\\
    {i}&{d_i}&{t_i}&{\gamma_i}&{n_i}&{d_i}&{t_i}&{\gamma_i}\\\hline
    {14}&2&3&\gamma &{{n_{14}}/(2\gamma)  + {n_{17}}/(2{\gamma ^2})}&0&3&0\\
    {15}&2&4&\gamma &{{n_{15}}/(2\gamma)  + {n_{17}}/(2{\gamma ^2})}&0&4&0\\
    {16}&2&5&\gamma &{{n_{16}}/(2\gamma)  + {n_{17}}/(2{\gamma ^2})}&0&5&0\\
    {17}&4&3&\gamma &{ - {n_{14}}/(2\gamma)  - {n_{17}}/(2{\gamma ^2})}&0&3&\gamma \\
    {18}&4&4&\gamma &{ - {n_{15}}/(2\gamma)  - {n_{18}}/(2{\gamma ^2})}&0&4&\gamma \\
    {19}&4&5&\gamma &{ - {n_{16}}/(2\gamma)  - {n_{19}}/(2{\gamma ^2})}&0&5&\gamma \\
    {20}&{}&{}&{}&{ - {n_{17}}/(2{\gamma})}&2&3&\gamma \\
    {21}&{}&{}&{}&{ - {n_{18}}/(2{\gamma})}&2&4&\gamma \\
    {22}&{}&{}&{}&{ - {n_{19}}/(2{\gamma})}&2&5&\gamma 
    \end{array}
    
.. warning::

    If the terms in the EOS are in terms of :math:`T` and :math:`\rho` rather than :math:`\tau` and :math:`\delta`, make sure to multiply appropriately by the critical densities in the exponential term.  For instance in Polt paper, the first constant should be :math:`n_{14}\rho_c^2/(2\gamma)+n_{17}\rho_c^4/(2\gamma^2)/T_c^3` Be careful!

