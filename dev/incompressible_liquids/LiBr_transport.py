import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy.ma as ma
from numpy.random import uniform, seed
# make up some randomly distributed data
# seed(1234)
#npts = 200
#x = uniform(-2,2,npts)
#y = uniform(-2,2,npts)
#z = x*np.exp(-x**2-y**2)
# define grid.
#xi = np.linspace(-2.1,2.1,100)
#yi = np.linspace(-2.1,2.1,100)
# grid the data.

Ti = np.linspace(300, 500)
xi = np.linspace(0.4, 0.7)

Tv = [312.9, 333.0, 353.2, 373.9, 393.3, 413.2, 433.0, 453.2, 472.5]
xv = [0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45]
zv = [1.919, 1.363, 1.054, 0.854, 0.725, 0.632, 0.564, 0.516, 0.480]

Tv += [314.9, 333.2, 353.7, 373.2, 393.3, 412.3, 432.6, 453.2, 472.6]
xv += [0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50]
zv += [2.446, 1.822, 1.400, 1.131, 0.942, 0.813, 0.712, 0.650, 0.589]

Tv += [314.2, 333.6, 353.1, 373.1, 393.5, 412.7, 433.0, 453.2, 472.6]
xv += [0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55]
zv += [3.418, 2.476, 1.895, 1.503, 1.228, 1.037, 0.900, 0.796, 0.720]

Tv += [316.0, 333.5, 353.6, 372.9, 394.3, 412.7, 433.0, 453.2, 472.9]
xv += [0.599, 0.599, 0.599, 0.599, 0.599, 0.599, 0.599, 0.599, 0.599]
zv += [4.952, 3.603, 2.657, 2.066, 1.634, 1.372, 1.164, 1.018, 0.886]

Tv += [333.1, 353.2, 373.6, 393.8, 413.2, 432.9, 452.9, 472.8]
xv += [0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63]
zv += [4.964, 3.537, 2.665, 2.090, 1.720, 1.463, 1.245, 1.099]

Tv += [333.4, 352.4, 372.3, 393.6, 413.2, 433.4, 452.6, 472.7]
xv += [0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65]
zv += [5.680, 4.158, 3.095, 2.372, 1.925, 1.633, 1.405, 1.240]

zi = griddata((Tv, xv), zv, (Ti[None, :], xi[:, None]), method='linear')  # , method='cubic')

#Tv += [000.0, 000.0, 000.0, 000.0, 000.0, 000.0, 000.0, 000.0, 000.0]
#xv += [ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00]

x = Tv
y = xv
z = zv

yi = xi
xi = Ti

zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')
# contour the gridded data, plotting dots at the randomly spaced data points.
CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.jet)
plt.colorbar()  # draw colorbar
# plot data points.
plt.scatter(x, y, marker='o', c='b', s=5)
# plt.xlim(-2,2)
# plt.ylim(-2,2)
#plt.title('griddata test (%d points)' % npts)
plt.show()
