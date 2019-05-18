import matplotlib
matplotlib.use('Qt4Agg')

import matplotlib.pyplot
from CoolProp.Plots.Plots import drawIsoLines, getIsoLines, drawLines

fluid = "n-Pentane"

fig, ((ax1, ax2)) = matplotlib.pyplot.subplots(1, 2, sharey='row')

drawIsoLines(fluid, 'Ts', 'Q', iValues=[0.0, 1.0], axis=ax1)  # for predefined styles
drawIsoLines(fluid, 'Ts', 'Q', iValues=[0.3, 0.7], axis=ax1)  # for predefined styles

# Get the data points
saturation = getIsoLines(fluid, 'Ts', 'Q', [0.0, 1.0], axis=ax2)
quality = getIsoLines(fluid, 'Ts', 'Q', [0.3, 0.7], axis=ax2)
# define custom styles
plt_kwargs = {"color": "green", "linewidth": 1.5}
drawLines(fluid, saturation, ax2, plt_kwargs=plt_kwargs)
plt_kwargs = {"color": "red", "linewidth": 0.75}
drawLines(fluid, quality, ax2, plt_kwargs=plt_kwargs)

matplotlib.pyplot.show()
