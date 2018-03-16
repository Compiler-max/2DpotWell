"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import matplotlib as mpl

from util_2dPW_Interf import read_Niu_shift

out_file = 'basic_animation.mp4'

plotState = 1

aUtoAngstrm 	= 0.52917721092
aX 				= 10.0*aUtoAngstrm
aY 				= 5.00*aUtoAngstrm



xlim=[-np.pi/aX, np.pi/aX]
ylim=[-np.pi/aY, np.pi/aY]



frame_delay_in_ms = 20


#sequenital colormaps
cmap 			= mpl.cm.viridis
#cmap 			= mpl.cm.plasma
#cmap 			= mpl.cm.inferno
#cmap 			= mpl.cm.magma

#diverging colormaps
#cmap			= mpl.cm.coolwarm

#well known
#cmap			= mpl.cm.jet


# First set up the figure, the axis, and the plot element we want to animate
#fig = plt.figure()
#ax = plt.axes(xlim=xlim, ylim=ylim)
#line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,


# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(0, 2, 100)
    y = np.sin(2 * np.pi * (x - 0.01 * i))
    line.set_data(x, y)
    return line,


## call the animator.  blit=True means only re-draw the parts that have changed.
#anim = animation.FuncAnimation(fig, animate, init_func=init,
#                               frames=200, interval=frame_delay_in_ms, blit=True)
#
## save the animation as an mp4.  This requires ffmpeg or mencoder to be
## installed.  The extra_args ensure that the x264 codec is used, so that
## the video can be embedded in html5.  You may need to adjust this for
## your system: for more information, see
## http://matplotlib.sourceforge.net/api/animation_api.html
##anim.save(out_file, fps=30, extra_args=['-vcodec', 'libx264'])
#
#plt.show()
#


#now try contour animation
fig = plt.figure()
ax = plt.axes(xlim=xlim, ylim=ylim)
#cont, = ax.contourf([], [], [], 500)



def my_contour_init(): 
	cont.set_data([],[],[])
	return cont,

def my_contour_test(i):
	do_f2, f2_xi, f2_yi, f2_plot_x, f2_plot_y, f2_min, f2_max,f2_bz= read_Niu_shift(plotState, "./a0.0/f2response.txt")

	zmin = f2_min
	zmax = f2_max

	print('try to plot:',f2_plot_x)
	plt.title(i)
	cont	= plt.contourf(f2_xi, f2_yi, f2_plot_x, 15, cmap=cmap,vmax=zmax, vmin=zmin)
	
	
	return cont





anim = animation.FuncAnimation(fig, my_contour_test,
                               frames=5, interval=frame_delay_in_ms, blit=True)

plt.show()

