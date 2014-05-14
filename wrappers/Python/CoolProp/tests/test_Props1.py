import unittest
from CoolProp.CoolProp import Props
import CoolProp

class Props1BadInputParameters(unittest.TestCase):
    """ All fluids, all parameters """
    def testEmptyFluid(self):
        self.assertRaises(ValueError,Props,'','Tcrit')
    def testIntegerFluid(self):
        self.assertRaises(TypeError,Props,1,'Tcrit')
    def testFloatFluid(self):
        self.assertRaises(TypeError,Props,1.0,'Tcrit')
    
    def testEmptyParam(self):
        self.assertRaises(ValueError,Props,'R134a','')
    def testBadParam(self):
        self.assertRaises(ValueError,Props,'R134a','R134a')
        
def testAllCoolPropPairs():
    for fluid in CoolProp.__fluids__:
        for param in ["Ttriple","Tcrit","pcrit","ptriple","Tmin",
                      "molemass","rhocrit","accentric","PHASE_LIQUID",
                      "PHASE_GAS","PHASE_SUPERCRITICAL","PHASE_TWOPHASE",
                      "ODP","GWP100"]:
            yield check_param,fluid,param

def check_param(fluid, param):
    val = Props(fluid,param)
    assert val > -2
    assert val < 1e9

if __name__=='__main__':
    import nose
    nose.runmodule()