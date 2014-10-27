
from CoolProp import CoolProp as cp
from CoolProp import __version__ as cpv
from CoolProp import __fluids__ as cpf
from sys import version

def cpversion():
	return cpv

def props(cp1,cp2,cp3,cp4,cp5,cp6):
	return cp.Props(cp1,cp2,cp3,cp4,cp5,cp6)

"""if __name__=='__main__':
	print props('P','T',300,'Q',1,'Propane')

 """
