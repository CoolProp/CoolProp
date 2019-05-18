import matplotlib.pyplot as plt
from CoolProp.Plots import Ph
from CoolProp.Plots.SimpleCycles import SimpleCycle

fig = plt.figure(figsize=(6, 1))
ax = fig.add_axes((0, 0, 1, 1))
ax.set_xlim(0, 6)
ax.set_ylim(0, 2)
plt.text(3, 1, 'CoolProp', size=50, name='Times', ha='center', va='center')
ax.axis('off')

ax2 = fig.add_axes((0, 0, 0.3, 1))
Ph('R410A', lw=2, color='white')
SimpleCycle(Ref='R410A', Te=280, Tc=310, DTsh=5, DTsc=5, eta_a=0.7, Ts_Ph='Ph', axis=ax2, lw=2, color='b')
ax2.set_xlim(-6, 600)
ax2.set_ylim(-1000, 6000)
ax2.axis('off')
plt.draw()
plt.savefig('_static/header.png', transparent=True)
plt.show()
