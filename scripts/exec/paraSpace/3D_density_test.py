#	https://stackoverflow.com/questions/25286811/how-to-plot-a-3d-density-map-in-python-with-matplotlib
#
#	http://docs.enthought.com/mayavi/mayavi/
#
#	Mayavi: 3D scientific data visualization and plotting in Python


import numpy as np
from scipy import stats
from mayavi import mlab

mu, sigma = 0, 0.1 
x = 10*np.random.normal(mu, sigma, 5000)
y = 10*np.random.normal(mu, sigma, 5000)
z = 10*np.random.normal(mu, sigma, 5000)

xyz = np.vstack([x,y,z])
kde = stats.gaussian_kde(xyz)
density = kde(xyz)


print('x length= ',len(x))
print('y length= ',len(y))
print('z length= ',len(z))
print('density length=',len(density))


# Plot scatter with mayavi
figure = mlab.figure('DensityPlot')
pts = mlab.points3d(x, y, z, density, scale_mode='none', scale_factor=0.07)
mlab.axes()
mlab.show()
