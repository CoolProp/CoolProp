import unittest
from CoolProp.CoolProp import Props
import CoolProp
import numpy as np
       
def test_input_types():
    for Fluid in ['Water']:
        for Tvals in [0.5*Props(Fluid,'Tmin')+0.5*Props(Fluid,'Tcrit'),
                      [Props(Fluid,'Tmin')+1e-5,Props(Fluid,'Tcrit')-1e-5],
                      np.linspace(Props(Fluid,'Tmin')+1e-5, Props(Fluid,'Tcrit')-1e-5,30)
                      ]:
            yield check_type, Fluid, Tvals

def check_type(fluid, Tvals):
    Props('P','T',Tvals,'Q',0,fluid)
    
class PropsFailures(unittest.TestCase):
    def testUnmatchedLengths(self):
        self.assertRaises(TypeError,Props,'P','T',[280,290,300],'Q',[0,1],'R134a')
    def testMatrix(self):
        self.assertRaises(TypeError,Props,'P','T',np.array([280,290,300,280,290,300]).reshape(2,3),'Q',np.array([0,0.5,1,0.0,0.5,1]).reshape(2,3),'R134a')

if __name__=='__main__':
    import nose
    nose.runmodule()