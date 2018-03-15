import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import griddata
from potwellInterface import read_Niu_shift
from potwellInterface import getData


#Perceptually Uniform Sequential colormaps






def plotNiuColor(searchDir, plotState, cmap=mpl.cm.viridis,
						plot_titles=False, title_size=12, 
						plot_k_labels=True, k_label_size=10,
						plot_descriptor=False, descriptor_size=14):
	f2_path	= searchDir+'/f2response.txt'
	f3_path = searchDir+'/f3response.txt'

	do_f2, f2_xi, f2_yi, f2_plot_x, f2_plot_y, f2_min, f2_max,f2_bz= read_Niu_shift(plotState, f2_path)
	do_f3, f3_xi, f3_yi, f3_plot_x, f3_plot_y, f3_min, f3_max,f3_bz= read_Niu_shift(plotState, f3_path)
	
	if abs(f2_bz-f3_bz) > 1e-8:
		print(searchDir+'/: WARNING different fields detected f2='+str(f2_bz)+' f3='+str(f3_bz))
	bz = f2_bz

	#get total max min of data
	zmin = min(f2_min,f3_min)
	zmax = max(f2_max,f3_max)

	# Plot each slice as an independent subplot
	fig, axes = plt.subplots(nrows=2, ncols=2, sharey='row')
	# Make an axis for the colorbar on the right side

	
	if do_f2:
		#
		# X SHIFT
		plt.subplot(221)
		CSx_f2 	= plt.contourf(f2_xi, f2_yi, f2_plot_x, 15, cmap=cmap,vmax=zmax, vmin=zmin)
		if plot_titles:
			plt.title('a_x (Ang)',fontsize=title_size)
		if plot_k_labels:
			plt.ylabel('ky (a.u.)',fontsize=k_label_size)
		#
		#
		# Y SHIFT
		plt.subplot(222)
		CSy_f2 	= plt.contourf(f2_xi,f2_yi, f2_plot_y, 15, cmap=cmap,vmax=zmax, vmin=zmin)
		if plot_titles:
			plt.title('a_y (Ang)',fontsize=title_size)
		
		plt.tick_params(axis='y',  which='both',direction='in',   left='off', right='off',     labelleft='off', labelright='off') 

	
	if do_f3:
		#
		# X SHIFT
		plt.subplot(223)
		CSx_f3 	= plt.contourf(f3_xi, f3_yi, f3_plot_x, 15, cmap=cmap, vmin=zmin,vmax=zmax)
		if plot_k_labels:
			plt.ylabel('ky (a.u.)',fontsize=k_label_size)
			plt.xlabel('kx (a.u.)',fontsize=k_label_size)
		#
		#
		# Y SHIFT
		plt.subplot(224)
		CSy_f3 	= plt.contourf(f3_xi,f3_yi, f3_plot_y, 15, cmap=cmap, vmin=zmin,vmax=zmax)
		if plot_k_labels:
			plt.xlabel('kx (a.u.)',fontsize=k_label_size)



	#set uniform colorbar
	print('colormap min_val=',zmin)
	print('colormap max_val=',zmax)
	
	
	
	#cb1.set_clim(zmin,zmax)
	a_Rashba = getData('alpha_rashba',searchDir+'/polOutput.txt')

	if plot_descriptor:
		plt.figtext(.25,.001,' n='+str(plotState)+' Bz='+str(bz)+'T'+' a_Rashba='+str(a_Rashba),fontsize=descriptor_size)

	#rescale the whole figure to add colorbar
	scale = .8
	plt.tight_layout(pad=1.25,rect=(0,0,scale,scale))

	cax = fig.add_axes([.81, 0.1, .03, scale*0.8])	# [left, bottom, width, height] 
	#cb1 = fig.colorbar(CSy_f2, cax=cax, cmap=cmap, norm=norm, label="1. order shift (ang)")
	norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
	cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
	



	plt.savefig(searchDir+'/a1_n'+str(plotState)+'_Bz'+str(bz)+'aRash'+str(a_Rashba)+'niuShift.pdf',bbox_inches='tight')
	#plt.show()

	return axes



#test
#searchDir 	= '.'
#fileName 	= 'f2response.txt'
#plotState	= 5
#
#ax = plotNiuColor(searchDir,  plotState)
#plt.show()
