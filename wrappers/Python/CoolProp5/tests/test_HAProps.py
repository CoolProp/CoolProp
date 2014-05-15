import unittest
import CoolProp
from CoolProp.HumidAirProp import HAProps
import numpy as np

def test_TRP():
    for R in np.linspace(0, 1, 11):
        for p in [101.325]:#np.linspace(0.1, 1000, 10):
            for T in np.linspace(220,370.15,10):
                for o in ['W','H','S','V']:
                    yield check_HAProps,o,'T',T,'R',R,'P',p

def check_HAProps(*args):
    val = HAProps(*args)
        
def test_input_types():

    pairs = [
             [300, 0.003],
              [[300,305],0.003],
              [np.linspace(300,305,6), 0.003],
              [300, [0.003, 0.0034]],
              [300, np.linspace(0.003,0.0034)]
             ]
    
    for T,w in pairs:
        yield check_type, T, w

def check_type(Tvals, wvals):
    HAProps('H','T',Tvals,'P',101.325,'W', wvals)
    
class PropsFailures(unittest.TestCase):
    def testUnmatchedLengths(self):
        self.assertRaises(TypeError,HAProps,'H','T',[300,301,302],'P',101.325,'W', [0.003,0.004])
    
if __name__=='__main__':
    import nose
    nose.runmodule()