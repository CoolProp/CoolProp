import numpy as np
from CPIncomp.CoefficientObjects import CoefficientData
from CPIncomp.BaseObjects import IncompressibleData




# TODO: Convert mass to volume polynomial

class SecCoolExample(CoefficientData):
    """ 
    Ethanol-Water mixture according to Melinder book
    Source: SecCool Software
    """ 
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "ExampleSecCool"
        self.description = "Methanol solution"
        #self.reference = "SecCool software"
        self.Tmax =  20 + 273.15
        self.Tmin = -50 + 273.15
        self.xmax = 0.5
        self.xmin = 0.0
        self.TminPsat =  20 + 273.15
    
        self.Tbase =  -4.48 + 273.15
        self.xbase =  31.57 / 100.0
       
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = self.convertSecCoolArray(np.array([
           960.24665800, 
          -1.2903839100, 
          -0.0161042520, 
          -0.0001969888, 
           1.131559E-05, 
           9.181999E-08,
          -0.4020348270,
          -0.0162463989,
           0.0001623301,
           4.367343E-06,
           1.199000E-08,
          -0.0025204776,
           0.0001101514,
          -2.320217E-07,
           7.794999E-08,
           9.937483E-06, 
          -1.346886E-06,
           4.141999E-08]))
        
        
        
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = self.convertSecCoolArray(np.array([
           3822.9712300,
          -23.122409500,
           0.0678775826,
           0.0022413893,
          -0.0003045332,
          -4.758000E-06, 
           2.3501449500, 
           0.1788839410, 
           0.0006828000, 
           0.0002101166, 
          -9.812000E-06, 
          -0.0004724176, 
          -0.0003317949, 
           0.0001002032, 
          -5.306000E-06, 
           4.242194E-05, 
           2.347190E-05, 
          -1.894000E-06]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = self.convertSecCoolArray(np.array([
           0.4082066700, 
          -0.0039816870, 
           1.583368E-05, 
          -3.552049E-07, 
          -9.884176E-10, 
           4.460000E-10, 
           0.0006629321, 
          -2.686475E-05, 
           9.039150E-07, 
          -2.128257E-08, 
          -5.562000E-10, 
           3.685975E-07, 
           7.188416E-08, 
          -1.041773E-08, 
           2.278001E-10, 
           4.703395E-08, 
           7.612361E-11, 
          -2.734000E-10]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_POLYNOMIAL
        self.viscosity.coeffs = self.convertSecCoolArray(np.array([
           1.4725525500, 
           0.0022218998, 
          -0.0004406139, 
           6.047984E-06, 
          -1.954730E-07, 
          -2.372000E-09, 
          -0.0411841566, 
           0.0001784479, 
          -3.564413E-06, 
           4.064671E-08, 
           1.915000E-08, 
           0.0002572862, 
          -9.226343E-07, 
          -2.178577E-08, 
          -9.529999E-10, 
          -1.699844E-06, 
          -1.023552E-07, 
           4.482000E-09]))

        self.T_freeze.type = self.T_freeze.INCOMPRESSIBLE_POLYOFFSET
        self.T_freeze.coeffs = self.convertSecCoolTfreeze(np.array([
           27.755555600,
          -22.973221700, 
          -1.1040507200, 
          -0.0120762281, 
          -9.343458E-05]))
        


class ZitrecAC(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self) 

        self.name = "ZiAC"
        self.description = "ZitrecAC in water (corrosion inhibitor)"
        self.reference = "SecCool Software"

        self.Tmin     =   0 + 273.15
        self.Tmax     = 100 + 273.15
        self.TminPsat = Tmax

        self.xmin     = 0.05
        self.xmax     = 0.50

        self.Tbase    = 50.00 + 273.15
        self.xbase    = 22.75 / 100.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = self.convertSecCoolArray(np.array([
          1003.4314200
          0.6164672840
         -0.0075340011
          8.227043E-05
          1.416356E-05
         -3.611589E-07
         -0.4849253070
         -0.0015163769
         -3.076387E-05
         -4.631673E-07
          3.958918E-08
         -0.0053161289
         -3.675038E-06
          2.245807E-06
         -5.921160E-08
         -2.222469E-05
         -1.016950E-06
          1.391098E-08
        cRho.clear(
        cRho = makeMatrix(tmpVector

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = self.convertSecCoolArray(np.array([
          4129.8221400
         -2.4602756000
          0.0024608374
         -0.0018878988
         -7.525318E-05
          3.812299E-06
          0.3207990110
         -0.0059380735
         -0.0007210516
          1.852057E-05
          1.280127E-07
          0.0099912322
          0.0001857846
         -2.266658E-06
         -5.937414E-09
         -0.0001641590
         -1.708329E-06
          8.577111E-08
        cHeat.clear(
        cHeat = makeMatrix(tmpVector

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = self.convertSecCoolArray(np.array([
          0.5997434470
         -0.0014518043
          7.240125E-06
         -6.458395E-07
         -1.035513E-08
          6.862018E-10
          0.0011059284
         -1.686714E-06
         -1.542610E-08
         -3.730822E-09
          1.286049E-11
         -5.567860E-06
          5.331112E-08
         -1.738845E-09
          2.641378E-11
         -1.195582E-08
         -3.329696E-10
          2.343560E-11
        cCond.clear(
        cCond = makeMatrix(tmpVector

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        self.viscosity.coeffs = self.convertSecCoolArray(np.array([
         -0.2833999350
          0.0113457408
         -0.0002173449
          9.645169E-06
          6.299118E-07
         -2.190184E-08
         -0.0185725092
         -8.290652E-05
         -2.186512E-06
          1.126354E-07
         -2.387450E-09
          0.0001256061
          1.147598E-06
          1.128849E-08
         -1.723001E-09
          7.142151E-07
         -5.140398E-09
          4.535194E-10
        cVisc.clear(
        cVisc = makeMatrix(tmpVector

        self.T_freeze.type = self.T_freeze.INCOMPRESSIBLE_POLYOFFSET
        self.T_freeze.coeffs = self.convertSecCoolTfreeze(np.array([
        cTfreeze.push_back( 22.750000000 // reference concentration in per cent
        cTfreeze.push_back(-2.2469093100
        cTfreeze.push_back(-0.0942887708
        cTfreeze.push_back( 0.0002636562
        cTfreeze.push_back( 9.008030E-07

    }

}


class IceSlurryEA(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "IceEA"
        self.description = "Ethanol-water mixture with slurry ice"
        self.reference = "SecCool Software"

        self.Tmin     = -35 + 273.15
        self.Tmax     = -10 + 273.15
        self.TminPsat = Tmax

        self.xmin     = 0.05
        self.xmax     = 0.35

        self.Tbase    = -22.5 + 273.15
        self.xbase    =  20.0 / 100.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = self.convertSecCoolArray(np.array([
          959.65328700
         -0.6063295290
         -1.249095E-05
         -2.175929E-05
          1.414141E-06
          8.888899E-08
          0.4964468440
         -0.0048152589
         -3.703183E-05
          7.619048E-07
          1.454545E-07
         -0.0105000000
          0.0001660431
          1.020408E-07
         -1.587301E-07
          0.0003601411
         -4.814815E-06
          1.763668E-08
        cRho.clear(
        cRho = makeMatrix(tmpVector

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = self.convertSecCoolArray(np.array([
          8.241062E+04
          4133.4319000
         -0.0525342683
         -0.0332405251
         -0.0007070707
          -9.555626E-05
         -175.73225900
         -19.750427100
         -0.5552331250
         -0.0024634921
          0.0001094372
         -20.031632700
         -1.1448185800
         -0.0049591837
          0.0003301587
         -1.3004938300
         -0.0241058201
          0.0012176367
        cHeat.clear(
        cHeat = makeMatrix(tmpVector

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = self.convertSecCoolArray(np.array([
          0.5466734710
          0.0098512311
          6.916942E-05
          4.423616E-07
          1.313131E-09
          1.333317E-10
          0.0064403362
          6.261462E-05
         -1.361660E-07
         -4.190476E-09
          9.350649E-11
          1.914626E-05
         -4.656462E-07
         -5.544218E-09
         -4.761923E-11
          1.093474E-08
         -1.322751E-09
          1.763668E-11
        cCond.clear(
        cCond = makeMatrix(tmpVector

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        self.viscosity.coeffs = self.convertSecCoolArray(np.array([
          3.9714277000
          0.0217778490
         -0.0002520703
          3.339905E-06
         -6.179120E-08
          2.654128E-09
         -0.0820149233
          3.060191E-06
         -1.789494E-07
         -3.463500E-09
         -5.507847E-11
         -0.0010965961
          5.631692E-08
          1.027114E-08
         -8.462427E-10
         1.855647E-05
          6.616363E-09
          1.185315E-10
        cVisc.clear(
        cVisc = makeMatrix(tmpVector

        cTfreeze.clear(
    }
    /// Define freezing point calculations
    double Tfreeze(double p, double x){
        return Tmin
    }
}


class IceSlurryPG(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "IcePG"
        self.description = "Propylene glycol-water mixture with slurry ice"
        self.reference = "SecCool Software"

        self.Tmin     = -45 + 273.15
        self.Tmax     = -10 + 273.15
        self.TminPsat = Tmax

        self.xmin     = 0.05
        self.xmax     = 0.35

        self.Tbase    = -27.5 + 273.15
        self.xbase    =  20.0 / 100.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = self.convertSecCoolArray(np.array([
         1026.0807100
         -1.5901706300
          0.0023632306
         -5.555564E-05
          3.787877E-07
          1.500003E-07
         -0.8565823200
          0.0158439454
         -4.885161E-05
          5.820106E-07
         -1.298701E-08
         -0.0205963719
          0.0002962207
         -3.287982E-07
          2.645504E-08
         -0.0002976431
          3.477633E-06
          1.827802E-08
        cRho.clear(
        cRho = makeMatrix(tmpVector

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = self.convertSecCoolArray(np.array([
          9.202807E+04
          4259.8599100
         -13.425892900
          0.1972224350
          0.0016666667
          2.333263E-05
         -1347.4964300
         -40.622180100
          1.2466366000
         -0.0089735450
         -0.0002542569
          9.1439909300
          0.3821466440
         -0.0142154195
         -0.0004074074
          0.4962000960
          0.0010591631
         -0.0011599808
        cHeat.clear(
        cHeat = makeMatrix(tmpVector

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = self.convertSecCoolArray(np.array([
          0.5308857010
          0.0102779335
          6.990287E-05
          4.527783E-07
          3.030303E-09
          9.999831E-11
          0.0037803301
          3.128644E-05
         -1.033550E-07
         -3.809524E-09
          4.617604E-11
          8.588549E-05
          3.319350E-07
         -5.147392E-09
         -1.164022E-10
          2.443001E-06
          8.542569E-09
         -1.558442E-10
        cCond.clear(
        cCond = makeMatrix(tmpVector

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        self.viscosity.coeffs = self.convertSecCoolArray(np.array([
          5.4382616400
          0.0178821732
         -0.0001873851
          2.395101E-06
         -3.751781E-08
          5.875223E-10
         -0.1513106190
         -6.271586E-06
          3.964457E-07
         -1.180828E-08
          9.672800E-11
         -0.0002684929
         -1.099919E-07
          6.150740E-09
         -1.902903E-10
         -4.725022E-06
         -1.027292E-10
          8.515353E-11
        cVisc.clear(
        cVisc = makeMatrix(tmpVector

        cTfreeze.clear(
    }
    /// Define freezing point calculations
    double Tfreeze(double p, double x){
        return Tmin
    }
}


class IceSlurryNA(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "IceNA"
        self.description = "Sodium chloride-water mixture with slurry ice"
        self.reference = "SecCool Software"

        self.Tmin     = -20 + 273.15
        self.Tmax     =  -5 + 273.15
        self.TminPsat = Tmax

        self.xmin     = 0.05
        self.xmax     = 0.35

        self.Tbase    = -12.5 + 273.15
        self.xbase    =  20.0 / 100.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = self.convertSecCoolArray(np.array([
          1081.6353100
         -2.4559523700
          0.0058152057
         -7.500013E-05
         -7.575759E-07
          1.666671E-07
         -5.6609963900
          0.1002726190
         -0.0004797330
          1.333333E-06
          3.636364E-08
         -0.0852857143
          0.0007904762
          1.428571E-06
          6.666668E-07
         -0.0037650794
          3.333333E-05
          6.984127E-07
        cRho.clear(
        cRho = makeMatrix(tmpVector

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = self.convertSecCoolArray(np.array([
          7.434384E+04
          3669.8467100
         -2.0844426400
          0.0312501929
         -0.0002727273
         -8.333396E-05
         -794.24689800
         -33.895515900
          0.3610772000
         -0.0016888889
          0.0001406061
          12.209523800
          0.3702381290
         -0.0099523810
          0.0001333331
         -0.1358730160
          0.0145714286
         -0.0014412698
        cHeat.clear(
        cHeat = makeMatrix(tmpVector

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = self.convertSecCoolArray(np.array([
          0.7579141770
          0.0124563700
          5.749080E-05
          2.263889E-07
         -7.575758E-10
          1.333333E-10
          0.0009894098
         -5.386429E-05
          2.049928E-07
          1.333333E-09
         -3.757576E-10
          8.761905E-06
         -9.531746E-07
          3.809524E-09
          2.222222E-10
         -9.777778E-07
         -5.904762E-08
         -1.269841E-10
        cCond.clear(
        cCond = makeMatrix(tmpVector

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        self.viscosity.coeffs = self.convertSecCoolArray(np.array([
          1.9270346900
          0.0216118832
         -0.0002820062
          4.476720E-06
         -9.032034E-08
          2.020892E-09
         -0.0713269985
         -2.729779E-05
          2.115509E-06
         -7.777973E-08
          2.272317E-09
          0.0003930707
         -7.171429E-07
          4.470017E-08
         -1.194034E-09
          2.156903E-05
         -1.127279E-08
          3.086759E-09
        cVisc.clear(
        cVisc = makeMatrix(tmpVector

        cTfreeze.clear(
    }
    /// Define freezing point calculations
    double Tfreeze(double p, double x){
        return Tmin
    }
}


class PK2000(DigitalData):
    """ 
    Lithium Bromide solution from CoolProp 4
    Source: Patek et al.
    """ 
    def __init__(self):
        DigitalData.__init__(self) 

        self.name = "PK2000"
        self.description = "Pekasol 2000 in water (Potassium acetate and formate)"
        self.reference = "SecCool Software"
        
        self.Tmin     = -62 + 273.15
        self.Tmax     = 100 + 273.15
        self.TminPsat = Tmax

        self.xminVol  = 0.36
        self.xmaxVol  = 1.00
        self.xbaseVol = 67.60 / 100.0  #volume percent!
        
        
        
        self.m2Vcoeffs= np.array([63.693536900, 1.0864767400, 0.0033972173, 1.986130E-05])

        self.Tbase    = 33.31 + 273.15
        self.xbase    = 67.60 / 100.0  #volume percent!
        
        self.temperature.data         = self.getTrange()
        self.concentration.data       = self.getxrange()
        
#        data = [self.density.data,self.specific_heat.data,self.saturation_pressure.data]
#        keys = ["D",              "C",                    "Psat"]
#        
#        import os
#        for i in range(len(keys)):
#            def func(T,x):
#                return CP.PropsSI(keys[i],'T',T,'P',1e8,self.name+"-{0:.4f}%".format(x*100.0))
#            #if os.path.isfile(self.getFile(key)): os.remove(self.getFile(key))
#            data[i] = self.getArray(func,keys[i])

        key = 'D'
        def funcD(T,x):
            return CP.PropsSI(key,'T',T,'P',1e8,self.name+"-{0:.4f}%".format(x*100.0))
        self.density.data = self.getArray(funcD,key)
        
        key = 'C'
        def funcC(T,x):
            return CP.PropsSI(key,'T',T,'P',1e8,self.name+"-{0:.4f}%".format(x*100.0))
        self.specific_heat.data = self.getArray(funcC,key)
        
        key = 'Psat'
        def funcP(T,x):
            return CP.PropsSI(key,'T',T,'P',1e8,self.name+"-{0:.4f}%".format(x*100.0))
        self.saturation_pressure.data = self.getArray(funcP,key)


    def setConcentrationConversion(self):
        xMass = np.array([36.0, 41.4, 46.7, 51.8, 56.8, 61.7, 66.5, 71.1, 75.6, 80.0, 84.2, 88.4, 92.4, 96.2, 100.0])/100.0
        xVolu = np.array([30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0])/100.0
        mass = IncompressibleData()
        mass.data = xMass
        volu = IncompressibleData()
        volu.data = xVolu
        
        std_coeffs = np.zeros((1,4))
        
        mass.coeffs = np.copy(std_coeffs)
        mass.type   = mass.INCOMPRESSIBLE_POLYNOMIAL
        mass.fitCoeffs(0.0, volu.data, 0.0, 0.0)
        
        volu.coeffs = np.copy(std_coeffs)
        volu.type   = volu.INCOMPRESSIBLE_POLYNOMIAL
        volu.fitCoeffs(0.0, mass.data, 0.0, 0.0)
        
        
        


    def massToVolume(self, xMass):
        coeffs = np.array([63.693536900, 1.0864767400, 0.0033972173, 1.986130E-05])
        # Redefine x_m to convert from mass to volume fraction
        # M1 = M-M_m
        # Vol% = A[1] + A[2]*M1 + A[3]*M1^2 + A[4]*M1^3
        xMassBase = 69.92
        xCalc     = xMass * 100.0 - xMassBase
        return np.polynomial.polynomial.polyval2d(xCalc, 0.0, c)
    
    def volumeToMass(self,xVolume):
        
        # The residual function            
        def fun(xArray,yArray,coefficients):
            """ 
            This function takes the coefficient array and solves
            for x using y data. It evaluates the function and returns
            the sum of the squared residuals if yArray is not
            equal to None
            """
            calculated = np.polynomial.polynomial.polyval(xArray, coefficients)
            if yArray==None: return calculated
            data       = yArray
            res        = np.sum(np.square(calculated-data))
            return res
        
        arguments  = (xData,yData)
        res = minimize(fun, coeffs_start, method=fit[counter], args=arguments, tol=tol)
        
        
    
        
        
        
        
        

class PK2000(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "PK2000"
        self.description = "Pekasol 2000 in water (Potassium acetate and formate)"
        self.reference = "SecCool Software"

        self.Tmin     = -62 + 273.15
        self.Tmax     = 100 + 273.15
        self.TminPsat = Tmax

        self.xmin     = 0.36
        self.xmax     = 1.00

        self.Tbase    = 33.31 + 273.15
        self.xbase    = 67.60 / 100.0 // volume percent!

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = self.convertSecCoolArray(np.array([
          1197.8306300
          2.7580390000
         -0.0046328716
         -2.118894E-05
          5.174717E-08
          2.265693E-09
         -0.5038076360
         -0.0020754530
          7.670317E-06
          1.443587E-08
         -1.451878E-09
         -0.0015712351
          2.546509E-05
         -8.010506E-08
         -4.948979E-10
          2.947414E-08
         -4.747692E-09
          8.927986E-11
        cRho.clear(
        cRho = makeMatrix(tmpVector

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = self.convertSecCoolArray(np.array([
          3012.5363200
         -11.345089400
          0.0512475571
         -0.0004261743
          5.582778E-06
          2.339332E-08
          1.8771264300
         -0.0024687047
         -0.0001136660
         -3.281309E-06
         -2.598223E-08
         -0.0117207757
          0.0001728598
          2.221588E-06
         -1.099247E-07
          6.566134E-05
          3.261243E-06
         -5.841138E-08
        cHeat.clear(
        cHeat = makeMatrix(tmpVector

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = self.convertSecCoolArray(np.array([
          0.5060902210
         -0.0015058953
          4.296707E-06
         -1.171421E-09
         -5.116965E-12
          8.612329E-14
          0.0009077322
         -4.052813E-06
          8.942927E-09
         -3.126326E-11
          1.729361E-12
          2.222018E-09
         -4.276078E-10
         -3.160239E-11
          1.220910E-12
          3.600639E-10
         -3.839049E-11
          7.540478E-13
        cCond.clear(
        cCond = makeMatrix(tmpVector

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        self.viscosity.coeffs = self.convertSecCoolArray(np.array([
          0.5631031820
          0.0147645079
         -1.024523E-05
         -7.780506E-07
          1.243805E-08
          7.494764E-12
         -0.0199399025
         -1.895556E-05
         -8.956373E-07
         -1.780877E-08
          9.735698E-11
          0.0001153086
          2.823599E-07
          7.608202E-09
         -3.413081E-10
         -1.580556E-06
         -1.693264E-08
         -1.648529E-10
        cVisc.clear(
        cVisc = makeMatrix(tmpVector

        self.T_freeze.type = self.T_freeze.INCOMPRESSIBLE_POLYOFFSET
        self.T_freeze.coeffs = self.convertSecCoolTfreeze(np.array([
        cTfreeze.push_back( 65.000000000 // reference concentration in per cent
        cTfreeze.push_back(-28.549321300
        cTfreeze.push_back(-0.7312947390
        cTfreeze.push_back(-0.0061085973
        cTfreeze.push_back(-8.714598E-06

        cMaVo.clear(
        cMaVo.push_back( 63.693536900
        cMaVo.push_back( 1.0864767400
        cMaVo.push_back( 0.0033972173
        cMaVo.push_back( 1.986130E-05

    }

    // Redefine x_m to convert from mass to volume fraction
    // M1 = M-M_m
    // Vol% = A[1] + A[2]*M1 + A[3]*M1^2 + A[4]*M1^3
    double getxInput(double curxValue){
        double xVolume = polyval(cMaVo, (curxValue*100.-69.92)
        return xVolume-xbase*100.0
    }


}        

