import matplotlib,numpy

from CoolProp.CoolProp import PropsSI
from scipy.optimize import newton

def SimpleCycle(Ref,Te,Tc,DTsh,DTsc,eta_a,Ts_Ph='Ph',skipPlot=False,axis=None):
    """
    This function plots a simple four-component cycle, on the current axis, or that given by the optional parameter *axis*

    Required parameters:

    * Ref : A string for the refrigerant
    * Te : Evap Temperature in K
    * Tc : Condensing Temperature in K
    * DTsh : Evaporator outlet superheat in K
    * DTsc : Condenser outlet subcooling in K
    * eta_a : Adiabatic efficiency of compressor (no units) in range [0,1]

    Optional parameters:

    * Ts_Ph : 'Ts' for a Temperature-Entropy plot, 'Ph' for a Pressure-Enthalpy
    * axis : An axis to use instead of the active axis
    * skipPlot : If True, won't actually plot anything, just print COP

    """
    T=numpy.zeros((6))
    h=numpy.zeros_like(T)
    p=numpy.zeros_like(T)
    s=numpy.zeros_like(T)
    T[1]=Te+DTsh
    pe=PropsSI('P','T',Te,'Q',1.0,Ref)
    pc=PropsSI('P','T',Tc,'Q',1.0,Ref)
    h[1]=PropsSI('H','T',T[1],'P',pe,Ref)
    s[1]=PropsSI('S','T',T[1],'P',pe,Ref)
    T2s=newton(lambda T: PropsSI('S','T',T,'P',pc,Ref)-s[1],T[1]+30)
    h2s=PropsSI('H','T',T2s,'P',pc,Ref)
    h[2]=(h2s-h[1])/eta_a+h[1]
    T[2]=PropsSI('T','H',h[2],'P',pc,Ref)
    s[2]=PropsSI('S','T',T[2],'P',pc,Ref)

    sbubble_c=PropsSI('S','P',pc,'Q',0,Ref)
    sdew_c=PropsSI('S','P',pc,'Q',1,Ref)
    sbubble_e=PropsSI('S','P',pe,'Q',0,Ref)
    sdew_e=PropsSI('S','P',pe,'Q',1,Ref)
    T[3]=Tc-DTsc
    h[3]=PropsSI('H','T',T[3],'P',pc,Ref)
    s[3]=PropsSI('S','T',T[3],'P',pc,Ref)
    h[4]=h[3]
    h[5]=h[1]
    s[5]=s[1]
    T[5]=T[1]
    p=[numpy.nan,pe,pc,pc,pe,pe]
    COP=(h[1]-h[4])/(h[2]-h[1])
    COPH=(h[2]-h[3])/(h[2]-h[1])

    hsatL=PropsSI('H','T',Te,'Q',0,Ref)
    hsatV=PropsSI('H','T',Te,'Q',1,Ref)
    ssatL=PropsSI('S','T',Te,'Q',0,Ref)
    ssatV=PropsSI('S','T',Te,'Q',1,Ref)
    vsatL=1/PropsSI('D','T',Te,'Q',0,Ref)
    vsatV=1/PropsSI('D','T',Te,'Q',1,Ref)
    x=(h[4]-hsatL)/(hsatV-hsatL)
    s[4]=x*ssatV+(1-x)*ssatL
    T[4]=x*Te+(1-x)*Te

    print(COP,COPH)
    if skipPlot==False:
        if axis==None:
            ax=matplotlib.pyplot.gca()
        if Ts_Ph in ['ph','Ph']:
            ax.plot(h,p)
        elif Ts_Ph in ['Ts','ts']:
            s=list(s)
            T=list(T)
            s.insert(5,sdew_e)
            T.insert(5,Te)
            s.insert(3,sbubble_c)
            T.insert(3,Tc)
            s.insert(3,sdew_c)
            T.insert(3,Tc)
            ax.plot(s[1::],T[1::],'b')
        else:
            raise TypeError('Type of Ts_Ph invalid')

def TwoStage(Ref,Q,Te,Tc,DTsh,DTsc,eta_oi,f_p,Tsat_ic,DTsh_ic,Ts_Ph='Ph',prints=False,skipPlot=False,axis=None,**kwargs):
    """
    This function plots a two-stage cycle, on the current axis, or that given by the optional parameter *axis*

    Required parameters:

    * Ref : Refrigerant [string]
    * Q : Cooling capacity [W]
    * Te : Evap Temperature [K]
    * Tc : Condensing Temperature [K]
    * DTsh : Evaporator outlet superheat [K]
    * DTsc : Condenser outlet subcooling [K]
    * eta_oi : Adiabatic efficiency of compressor (no units) in range [0,1]
    * f_p : fraction of compressor power lost as ambient heat transfer in range [0,1]
    * Tsat_ic : Saturation temperature corresponding to intermediate pressure [K]
    * DTsh_ic : Superheating at outlet of intermediate stage [K]

    Optional parameters:

    * Ts_Ph : 'Ts' for a Temperature-Entropy plot, 'Ph' for a Pressure-Enthalpy
    * prints : True to print out some values
    * axis : An axis to use instead of the active axis
    * skipPlot : If True, won't actually plot anything, just print COP

    """

    T=numpy.zeros((8))
    h=numpy.zeros_like(T)
    p=numpy.zeros_like(T)
    s=numpy.zeros_like(T)
    rho=numpy.zeros_like(T)
    T[0]=numpy.NAN
    s[0]=numpy.NAN
    T[1]=Te+DTsh
    pe=PropsSI('P','T',Te,'Q',1.0,Ref)
    pc=PropsSI('P','T',Tc,'Q',1.0,Ref)
    pic=PropsSI('P','T',Tsat_ic,'Q',1.0,Ref)
    Tbubble_c=PropsSI('T','P',pc,'Q',0,Ref)
    Tbubble_e=PropsSI('T','P',pe,'Q',0,Ref)

    h[1]=PropsSI('H','T',T[1],'P',pe,Ref)
    s[1]=PropsSI('S','T',T[1],'P',pe,Ref)
    rho[1]=PropsSI('D','T',T[1],'P',pe,Ref)
    T[5]=Tbubble_c-DTsc
    h[5]=PropsSI('H','T',T[5],'P',pc,Ref)
    s[5]=PropsSI('S','T',T[5],'P',pc,Ref)
    rho[5]=PropsSI('D','T',T[5],'P',pc,Ref)
    mdot=Q/(h[1]-h[5])

    rho1=PropsSI('D','T',T[1],'P',pe,Ref)
    h2s=PropsSI('H','S',s[1],'P',pic,Ref)
    Wdot1=mdot*(h2s-h[1])/eta_oi
    h[2]=h[1]+(1-f_p)*Wdot1/mdot
    T[2]=PropsSI('T','H',h[2],'P',pic,Ref)
    s[2]=PropsSI('S','T',T[2],'P',pic,Ref)
    rho[2]=PropsSI('D','T',T[2],'P',pic,Ref)
    T[3]=288
    p[3]=pic
    h[3]=PropsSI('H','T',T[3],'P',pic,Ref)
    s[3]=PropsSI('S','T',T[3],'P',pic,Ref)
    rho[3]=PropsSI('D','T',T[3],'P',pic,Ref)
    rho3=PropsSI('D','T',T[3],'P',pic,Ref)
    h4s=PropsSI('H','T',s[3],'P',pc,Ref)
    Wdot2=mdot*(h4s-h[3])/eta_oi
    h[4]=h[3]+(1-f_p)*Wdot2/mdot
    T[4]=PropsSI('T','H',h[4],'P',pc,Ref)
    s[4]=PropsSI('S','T',T[4],'P',pc,Ref)
    rho[4]=PropsSI('D','T',T[4],'P',pc,Ref)

    sbubble_e=PropsSI('S','T',Tbubble_e,'Q',0,Ref)
    sbubble_c=PropsSI('S','T',Tbubble_c,'Q',0,Ref)
    sdew_e=PropsSI('S','T',Te,'Q',1,Ref)
    sdew_c=PropsSI('S','T',Tc,'Q',1,Ref)

    hsatL=PropsSI('H','T',Tbubble_e,'Q',0,Ref)
    hsatV=PropsSI('H','T',Te,'Q',1,Ref)
    ssatL=PropsSI('S','T',Tbubble_e,'Q',0,Ref)
    ssatV=PropsSI('S','T',Te,'Q',1,Ref)
    vsatL=1/PropsSI('D','T',Tbubble_e,'Q',0,Ref)
    vsatV=1/PropsSI('D','T',Te,'Q',1,Ref)
    x=(h[5]-hsatL)/(hsatV-hsatL)
    s[6]=x*ssatV+(1-x)*ssatL
    T[6]=x*Te+(1-x)*Tbubble_e
    rho[6]=1.0/(x*vsatV+(1-x)*vsatL)

    h[6]=h[5]
    h[7]=h[1]
    s[7]=s[1]
    T[7]=T[1]
    p=[numpy.nan,pe,pic,pic,pc,pc,pe,pe]
    COP=Q/(Wdot1+Wdot2)
    RE=h[1]-h[6]

    if prints==True:
        print('x5:',x)
        print('COP:', COP)
        print('COPH', (Q+Wdot1+Wdot2)/(Wdot1+Wdot2))
        print(T[2]-273.15,T[4]-273.15,p[2]/p[1],p[4]/p[3])
        print(mdot,mdot*(h[4]-h[5]),pic)
        print('Vdot1',mdot/rho1,'Vdisp',mdot/rho1/(3500/60.)*1e6/0.7)
        print('Vdot2',mdot/rho3,'Vdisp',mdot/rho3/(3500/60.)*1e6/0.7)
        print(mdot*(h[4]-h[5]),Tc-273.15)
        for i in range(1,len(T)-1):
            print('%d & %g & %g & %g & %g & %g \\\\' %(i,T[i]-273.15,p[i],h[i],s[i],rho[i]))
    else:
        print(Tsat_ic,COP)

    if skipPlot==False:
        if axis==None:
            ax=matplotlib.pyplot.gca()
        else:
            ax=axis
        if Ts_Ph in ['ph','Ph']:
            ax.plot(h,p)
        elif Ts_Ph in ['Ts','ts']:
            s_copy=s.copy()
            T_copy=T.copy()
            for i in range(1,len(s)-1):
                ax.plot(s[i],T[i],'bo',mfc='b',mec='b')
                dT=[0,-5,5,-20,5,5,5]
                ds=[0,0.05,0,0,0,0,0]
                ax.text(s[i]+ds[i],T[i]+dT[i],str(i))

            s=list(s)
            T=list(T)
            s.insert(7,sdew_e)
            T.insert(7,Te)
            s.insert(5,sbubble_c)
            T.insert(5,Tbubble_c)
            s.insert(5,sdew_c)
            T.insert(5,Tc)

            ax.plot(s,T)
            s=s_copy
            T=T_copy
        else:
            raise TypeError('Type of Ts_Ph invalid')
    return COP

def EconomizedCycle(Ref,Qin,Te,Tc,DTsh,DTsc,eta_oi,f_p,Ti,Ts_Ph='Ts',skipPlot=False,axis=None,**kwargs):
    """
    This function plots an economized cycle, on the current axis, or that given by the optional parameter *axis*

    Required parameters:

    * Ref : Refrigerant [string]
    * Qin : Cooling capacity [W]
    * Te : Evap Temperature [K]
    * Tc : Condensing Temperature [K]
    * DTsh : Evaporator outlet superheat [K]
    * DTsc : Condenser outlet subcooling [K]
    * eta_oi : Adiabatic efficiency of compressor (no units) in range [0,1]
    * f_p : fraction of compressor power lost as ambient heat transfer in range [0,1]
    * Ti : Saturation temperature corresponding to intermediate pressure [K]

    Optional parameters:

    * Ts_Ph : 'Ts' for a Temperature-Entropy plot, 'Ph' for a Pressure-Enthalpy
    * axis : An axis to use instead of the active axis
    * skipPlot : If True, won't actually plot anything, just print COP

    """

    m=1

    T=numpy.zeros((11))
    h=numpy.zeros_like(T)
    p=numpy.zeros_like(T)
    s=numpy.zeros_like(T)
    rho=numpy.zeros_like(T)

    T[0]=numpy.NAN
    s[0]=numpy.NAN
    T[1]=Te+DTsh
    pe=PropsSI('P','T',Te,'Q',1.0,Ref)
    pc=PropsSI('P','T',Tc,'Q',1.0,Ref)
    pi=PropsSI('P','T',Ti,'Q',1.0,Ref)
    p[1]=pe
    h[1]=PropsSI('H','T',T[1],'P',pe,Ref)
    s[1]=PropsSI('S','T',T[1],'P',pe,Ref)
    rho[1]=PropsSI('D','T',T[1],'P',pe,Ref)
    h2s=PropsSI('H','S',s[1],'P',pi,Ref)
    wdot1=(h2s-h[1])/eta_oi
    h[2]=h[1]+(1-f_p[0])*wdot1
    p[2]=pi
    T[2]=T_hp(Ref,h[2],pi,T2s)
    s[2]=PropsSI('S','T',T[2],'P',pi,Ref)
    rho[2]=PropsSI('D','T',T[2],'P',pi,Ref)

    T[5]=Tc-DTsc
    h[5]=PropsSI('H','T',T[5],'P',pc,Ref)
    s[5]=PropsSI('S','T',T[5],'P',pc,Ref)
    rho[5]=PropsSI('D','T',T[5],'P',pc,Ref)

    p[5]=pc
    p[6]=pi
    h[6]=h[5]

    p[7]=pi
    p[8]=pi
    p[6]=pi
    T[7]=Ti
    h[7]=PropsSI('H','T',Ti,'Q',1,Ref)
    s[7]=PropsSI('S','T',Ti,'Q',1,Ref)
    rho[7]=PropsSI('D','T',Ti,'Q',1,Ref)
    T[8]=Ti
    h[8]=PropsSI('H','T',Ti,'Q',0,Ref)
    s[8]=PropsSI('S','T',Ti,'Q',0,Ref)
    rho[8]=PropsSI('D','T',Ti,'Q',0,Ref)
    x6=(h[6]-h[8])/(h[7]-h[8]) #Vapor Quality
    s[6]=s[7]*x6+s[8]*(1-x6)
    rho[6]=1.0/(x6/rho[7]+(1-x6)/rho[8])
    T[6]=Ti

    #Injection mass flow rate
    x=m*(h[6]-h[8])/(h[7]-h[6])


    p[3]=pi
    h[3]=(m*h[2]+x*h[7])/(m+x)
    T[3]=T_hp(Ref,h[3],pi,T[2])
    s[3]=PropsSI('S','T',T[3],'P',pi,Ref)
    rho[3]=PropsSI('D','T',T[3],'P',pi,Ref)
    T4s=newton(lambda T: PropsSI('S','T',T,'P',pc,Ref)-s[3],T[2]+30)
    h4s=PropsSI('H','T',T4s,'P',pc,Ref)
    p[4]=pc
    wdot2=(h4s-h[3])/eta_oi
    h[4]=h[3]+(1-f_p[1])*wdot2
    T[4]=T_hp(Ref,h[4],pc,T4s)
    s[4]=PropsSI('S','T',T[4],'P',pc,Ref)
    rho[4]=PropsSI('D','T',T[4],'P',pc,Ref)

    p[9]=pe
    h[9]=h[8]
    T[9]=Te
    hsatL_e=PropsSI('H','T',Te,'Q',0,Ref)
    hsatV_e=PropsSI('H','T',Te,'Q',1,Ref)
    ssatL_e=PropsSI('S','T',Te,'Q',0,Ref)
    ssatV_e=PropsSI('S','T',Te,'Q',1,Ref)
    vsatL_e=1/PropsSI('D','T',Te,'Q',0,Ref)
    vsatV_e=1/PropsSI('D','T',Te,'Q',1,Ref)
    x9=(h[9]-hsatL_e)/(hsatV_e-hsatL_e) #Vapor Quality
    s[9]=ssatV_e*x9+ssatL_e*(1-x9)
    rho[9]=1.0/(x9*vsatV_e+(1-x9)*vsatL_e)

    s[10]=s[1]
    T[10]=T[1]
    h[10]=h[1]
    p[10]=p[1]

    Tbubble_e=Te
    Tbubble_c=Tc
    sbubble_e=PropsSI('S','T',Tbubble_e,'Q',0,Ref)
    sbubble_c=PropsSI('S','T',Tbubble_c,'Q',0,Ref)
    sdew_e=PropsSI('S','T',Te,'Q',1,Ref)
    sdew_c=PropsSI('S','T',Tc,'Q',1,Ref)

    Wdot1=m*wdot1
    Wdot2=(m+x)*wdot2
    if skipPlot==False:
        if axis==None:
            ax=matplotlib.pyplot.gca()
        else:
            ax=axis
        if Ts_Ph in ['ph','Ph']:
            ax.plot(h,p)
            ax.set_yscale('log')
        elif Ts_Ph in ['Ts','ts']:
            ax.plot(numpy.r_[s[7],s[3]],numpy.r_[T[7],T[3]],'b')
            s_copy=s.copy()
            T_copy=T.copy()
            dT=[0,-5,5,-12,5,12,-12,0,0,0]
            ds=[0,0.05,0.05,0,0.05,0,0.0,0.05,-0.05,-0.05]
            for i in range(1,len(s)-1):
                ax.plot(s[i],T[i],'bo',mfc='b',mec='b')
                ax.text(s[i]+ds[i],T[i]+dT[i],str(i),ha='center',va='center')

            s=list(s)
            T=list(T)
            s.insert(10,sdew_e)
            T.insert(10,Te)
            s.insert(5,sbubble_c)
            T.insert(5,Tbubble_c)
            s.insert(5,sdew_c)
            T.insert(5,Tc)
            ax.plot(s,T,'b')

            s=s_copy
            T=T_copy
        else:
            raise TypeError('Type of Ts_Ph invalid')

    COP=m*(h[1]-h[9])/(m*(h[2]-h[1])+(m+x)*(h[4]-h[3]))
    for i in range(1,len(T)-1):
            print('%d & %g & %g & %g & %g & %g \\\\' %(i,T[i]-273.15,p[i],h[i],s[i],rho[i]))
    print(x,m*(h[1]-h[9]),(m*(h[2]-h[1])+(m+x)*(h[4]-h[3])),COP)
    mdot=Qin/(h[1]-h[9])
    mdot_inj=x*mdot
    print('x9',x9,)
    print('Qcond',(mdot+mdot_inj)*(h[4]-h[5]),'T4',T[4]-273.15)
    print(mdot,mdot+mdot_inj)
    f=3500/60.
    eta_v=0.7
    print('Vdisp1: ',mdot/(rho[1]*f*eta_v)*1e6,'cm^3')
    print('Vdisp2: ',(mdot+mdot_inj)/(rho[1]*f*eta_v)*1e6,'cm^3')
    return COP

if __name__=='__main__':
    from CoolProp.Plots import Ph,Ts

    Ref='R290'
    fig=matplotlib.pyplot.figure(figsize=(4,3))
    ax=fig.add_axes((0.15,0.15,0.8,0.8))
    Ph(Ref,Tmin=273.15-30,hbounds=[0,600],axis=ax)
    COP=TwoStage('Propane',10000,273.15-5,273.15+43.3,5,7,0.7,0.3,15+273.15,3,prints = True)
    matplotlib.pyplot.show()

    Ref='R290'
    fig=matplotlib.pyplot.figure(figsize=(4,3))
    ax=fig.add_axes((0.15,0.15,0.8,0.8))
    Ph(Ref,Tmin=273.15-30,hbounds=[0,600],axis=ax)
    COP=SimpleCycle(Ref,273.15-5,273.15+45,5,7,0.7,Ts_Ph='Ph')
    matplotlib.pyplot.show()

    Ref='R410A'
    fig=matplotlib.pyplot.figure(figsize=(4,3))
    ax=fig.add_axes((0.15,0.15,0.8,0.8))
    Ts(Ref,Tmin=273.15-100,sbounds=[0,600],axis=ax)
    COP=SimpleCycle(Ref,273.15-5,273.15+45,5,7,0.7,Ts_Ph='Ts')
    matplotlib.pyplot.show()




##     for x in numpy.linspace(0,1):
##         Ref='REFPROP-MIX:R152A[%g]&R32[%g]' %(x,1-x)
##         COP=SimpleCycle(273.15+8,273.15+44,5,7,0.7,skipPlot=True,Ts_Ph='Ph')
##     matplotlib.pyplot.show()
