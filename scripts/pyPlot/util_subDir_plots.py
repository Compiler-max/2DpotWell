import os
import os.path
import matplotlib.pyplot as plt
import matplotlib as mpl
#my stuff
from plot_bandStruct 	import plotBands
from plot_niu_2dcolor 	import plotNiuColor


#BANDSTRUCTURE PLOT
minEn 			= -50.0
maxEn 			= +50.0
show_Bext_box 	= False
show_Rash_box	= True

#show plots after saving theme
show_plots		= False


#NIU RESPONSE PLOT
nWfs 	= 6

plot_titles		= False
plot_k_labels	= True
plot_descriptor	= False

#sequenital colormaps
cmap 			= mpl.cm.viridis
#cmap 			= mpl.cm.plasma
#cmap 			= mpl.cm.inferno
#cmap 			= mpl.cm.magma

#diverging colormaps
#cmap			= mpl.cm.coolwarm

#well known
#cmap			= mpl.cm.jet



def plot_Dir(dirpath=".",show_Plots=False):
	print('search directory ',dirpath)
	info_path =dirpath+'/polOutput.txt'
	try:
		info_file = open(info_path)
	except IOError:
		print('IOError: did not find file ',info_path)
		axBands = []
	else:
		info_file.close()
		#
		#BANDSTRUCTURE
		axBands = plotBands(dirpath, dirpath, minEn, maxEn, show_Bext_box, show_Rash_box)
		if show_Plots:
			plt.show()
		#
		#FIRST ORDER
		for n in range(1,nWfs+1):
			axNiu = plotNiuColor(dirpath,n,cmap,plot_titles=plot_titles, plot_k_labels=plot_k_labels, plot_descriptor=plot_descriptor)
		if show_Plots:
			plt.show()
		print('successfully ploted ',dirpath)
	finally:
		print('done with',dirpath)
		return axBands

	
		



#fig = plt.figure()
#ax = plt.axes(xlim=(0,7), ylim=(minEn,maxEn))
#line, = ax.plot([], [], lw=2)
#
## initialization function: plot the background of each frame
#def init():
#	line.set_data([], [])
#	return line,
#
#
## animation function.  This is called sequentially
#def animate(i):
#	x = np.linspace(0, 2, 100)
#	y = np.sin(2 * np.pi * (x - 0.01 * i))
#	line.set_data(x, y)
#	return line,
#
#
### call the animator.  blit=True means only re-draw the parts that have changed.
#anim = animation.FuncAnimation(fig, animate, init_func=init,
#                               frames=200, interval=frame_delay_in_ms, blit=True)
##
### save the animation as an mp4.  This requires ffmpeg or mencoder to be
### installed.  The extra_args ensure that the x264 codec is used, so that
### the video can be embedded in html5.  You may need to adjust this for
### your system: for more information, see
### http://matplotlib.sourceforge.net/api/animation_api.html
###anim.save(out_file, fps=30, extra_args=['-vcodec', 'libx264'])
##
#plt.show()
##
#
#




axBands = []
for dirpath, dirnames, filenames in os.walk("."):
	#create the pdf plot
	if 'pycache' in dirpath:
		print('exclude diretory ',dirpath)
	else:
		axBands.append(	plot_Dir(dirpath,show_plots) )
		
	
print('len axBands=',len(axBands))		
		

		

