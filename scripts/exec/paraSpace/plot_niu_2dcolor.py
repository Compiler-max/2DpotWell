import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import griddata
from util_2dPW_Interf import read_Niu_shift
from util_2dPW_Interf import getData


#Perceptually Uniform Sequential colormaps








def plotNiuColor(searchDir, plotState, 
						aX_ang=0,aY_ang=0, plot_percent=False,
						plot_essin=True,
						cmap=mpl.cm.viridis,
						plot_titles=False, title_size=12, 
						plot_k_labels=True, k_label_size=10,
						plot_descriptor=False, descriptor_size=14,
						save_dir="."):
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

	if plot_essin:
		#essin flips sign of response
		f2_plot_x = - f2_plot_x
		f2_plot_y = - f2_plot_y
		#
		f3_plot_x = - f3_plot_x
		f3_plot_y = - f3_plot_y

		old_max = zmax
		zmax = - zmin
		zmin = - old_max


	
	
	# Plot each slice as an independent subplot
	fig, axes = plt.subplots(nrows=2, ncols=2, sharey='row')
	# Make an axis for the colorbar on the right side

	#only plot percent if the lattice vectors are given
	do_plot_percent = False
	if plot_percent:
		if (aX_ang is not 0) and (aY_ang is not 0):
			do_plot_percent = True

	if do_plot_percent:
		#x
		f2_plot_x = f2_plot_x / aX_ang 
		f3_plot_x = f3_plot_x / aX_ang
		#y
		f2_plot_y = f2_plot_y / aY_ang
		f3_plot_y = f3_plot_y / aY_ang
		
		zmin = min(f2_min/max(aX_ang,aY_ang),f3_min/max(aX_ang,aY_ang))
		zmax = max(f2_max/min(aX_ang,aY_ang),f3_max/min(aX_ang,aY_ang))

		zmin= -(min(abs(zmin),abs(zmax)))
		zmax= -zmin



	
	
	if do_f2:

		#
		# X SHIFT
		plt.subplot(221)
		CSx_f2 	= plt.contourf(f2_xi, f2_yi, f2_plot_x, 15, cmap=cmap,vmax=zmax, vmin=zmin)
		if plot_titles:
			if do_plot_percent:
				plt.title(r'$a_x (p_{\mathrm{q}})$',fontsize=title_size)
			else:
				plt.title(r'$a_x$ (Ang)',fontsize=title_size)
		if plot_k_labels:
			plt.ylabel('ky (a.u.)',fontsize=k_label_size)
		#
		#
		# Y SHIFT
		plt.subplot(222)
		CSy_f2 	= plt.contourf(f2_xi,f2_yi, f2_plot_y, 15, cmap=cmap,vmax=zmax, vmin=zmin)
		if plot_titles:
			if do_plot_percent:
				plt.title(r'$a_y (p_{\mathrm{q}})$',fontsize=title_size)
			else:
				plt.title(r'$a_y$ (Ang)',fontsize=title_size)
		
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
	


	if plot_essin:
		plt.savefig(save_dir+'/a1_n'+str(plotState)+'_Bz'+str(bz)+'aRash'+str(a_Rashba)+'essinShift.pdf',bbox_inches='tight')
	else:
		plt.savefig(save_dir+'/a1_n'+str(plotState)+'_Bz'+str(bz)+'aRash'+str(a_Rashba)+'niuShift.pdf',bbox_inches='tight')
	#plt.show()
	plt.close()

	return axes


#



def plotDir():
	#test
	searchDir 	= '.'
	nwfs 		= 6



	aUtoAngstrm = 0.52917721092
	aX 			= 10.0 * aUtoAngstrm
	aY			= 5.0 * aUtoAngstrm

	print('unit cell=(',aX,', ',aY,'),(Ang)')

	#
	for plotState in range(1,nwfs+1):
		print('plot state=',plotState)
		ax = plotNiuColor(searchDir,  plotState, aX_ang=aX,aY_ang=aY, plot_percent=True, plot_essin=False, plot_titles=True)
		plt.show()

plotDir()