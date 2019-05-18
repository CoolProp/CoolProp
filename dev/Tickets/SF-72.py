import matplotlib
matplotlib.use('Qt4Agg')
from CoolProp.Plots.Plots import Ph, drawIsoLines
from CoolProp.Plots.SimpleCycles import SimpleCycle
from CoolProp.CoolProp import Props
import matplotlib.pyplot as plt

Ref = 'Propane'
ax = Ph(Ref)
SimpleCycle(Ref, 260, 320, 5, 5, 0.7)
# ax.set_xlim([200,900])
# ax.set_ylim([300,530])
quality = drawIsoLines(Ref, 'ph', 'Q', [0.2, 0.4, 0.6, 0.8], axis=ax)
isobars = drawIsoLines(Ref, 'ph', 'T', [250.0, 300.0], axis=ax)
#isochores  = drawIsoLines(Ref, 'Ts', 'D', [2,    600]    , num=7, axis=ax)
plt.show()
