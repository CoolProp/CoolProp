from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropsPlot
import matplotlib.pylab as pl
from numpy import *
%matplotlib inline

eta_e_s = 0.88
eta_p_s = 0.5
T_max = 550 + 273.15
p_max = 13000.e3 
p_cd = 5.e3

def SimpleRankineCycle(T3, p3, p1, epsilon_e, epsilon_p, fluid='water'):
    h1 = PropsSI('H', 'P', p1, 'Q', 0., fluid)
    s1 = PropsSI('S', 'P', p1, 'Q', 0., fluid)
    T1 = PropsSI('T', 'P', p1, 'Q', 0., fluid)
    
    p2 = p3
    h2 = h1 + (PropsSI('H', 'P', p2, 'S', s1, fluid) - h1) / epsilon_p
    s2 = PropsSI('S', 'H', h2, 'P', p2, fluid)
    T2 = PropsSI('T', 'H', h2, 'P', p2, fluid)
    
    h3 = PropsSI('H', 'P', p3, 'T', T3, fluid)
    s3 = PropsSI('S', 'H', h3, 'P', p3, fluid)
    
    p4 = p1
    h4 = h3 - epsilon_e * (h3 - PropsSI('H', 'P', p4, 'S', s3, fluid))
    s4 = PropsSI('S', 'H', h4, 'P', p4, fluid)
    T4 = PropsSI('T', 'H', h4, 'P', p4, fluid)
    
    w_net = h3 - h4
    q_boiler = h3 - h2
    eta_c = w_net / q_boiler
    
    Ts = PropsPlot(fluid, 'Ts')
    Ts.draw_isolines('P', [p1, p3], num=10)
    Ts.set_axis_limits([0., 12., 200., 900.])
    
    ax = Ts.axis
    ax.text(s1/1000, T1,' 1', fontsize=10, rotation=0, color='r')
    ax.text(s2/1000, T2,' 2', fontsize=10, rotation=0, color='r')
    ax.text(s3/1000, T3,' 3', fontsize=10, rotation=0, color='r')
    ax.text(s4/1000, T4,' 4', fontsize=10, rotation=0, color='r')
    ax.text(8., 850., "Efficiency: %.1f%%" %(eta_c*100.))
    ax.text(8., 800., "Net work: %d kJ/kg" %(w_net/1000))
    ax.text(8., 750., "Heat input: %d kJ/kg" %(q_boiler/1000))
    return Ts

Ts = SimpleRankineCycle(T_max, p_max, p_cd, eta_e_s, eta_p_s, fluid="water")