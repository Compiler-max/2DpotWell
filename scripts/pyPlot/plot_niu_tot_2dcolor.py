import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import griddata
from util_2dPW_Interf import read_Niu_shift
from util_2dPW_Interf import getData


#Perceptually Uniform Sequential colormaps








def plotNiu_orb(searchDir, plotState, 
						aX_ang=0,aY_ang=0, plot_percent=False,
						plot_essin=True,
						cmap=mpl.cm.viridis,
						plot_titles=False, title_size=12, 
						plot_k_labels=True, k_label_size=10,
						plot_descriptor=False, descriptor_size=14,
						save_dir="."):
	f2_path	= searchDir+'/f2response.txt'
	f3_path = searchDir+'/f3response.txt'

	do_f2, f2_xi, f2_yi, f2_plot_x, f2_plot_y, f2_min, f2_max,f2_bz, raw_f2_kx, raw_f2_ky, raw_f2_x, raw_f2_y	=	read_Niu_shift(plotState, f2_path)
	do_f3, f3_xi, f3_yi, f3_plot_x, f3_plot_y, f3_min, f3_max,f3_bz, raw_f3_kx, raw_f3_ky, raw_f3_x, raw_f3_y	= 	read_Niu_shift(plotState, f3_path)
	
	if abs(f2_bz-f3_bz) > 1e-8:
		print(searchDir+'/: WARNING different fields detected f2='+str(f2_bz)+' f3='+str(f3_bz))
	bz = f2_bz


	if do_f2 and do_f3:
		#get total max min of data
		zmin = min(f2_min,f3_min)
		zmax = max(f2_max,f3_max)

		#adjust prefactor
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

		#sum F2 & F3
		raw_tot_x = raw_f2_x + raw_f3_x
		raw_tot_y = raw_f2_y + raw_f3_y

		tot_plot_x = griddata((raw_f2_kx,raw_f2_ky), raw_tot_x, (f2_xi[None,:], f2_yi[:,None]), method='cubic')
		tot_plot_y = griddata((raw_f2_kx,raw_f2_ky), raw_tot_y, (f2_xi[None,:], f2_yi[:,None]), method='cubic')

	
	
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
			tot_plot_x = tot_plot_x  / aX_ang
			#y
			f2_plot_y = f2_plot_y / aY_ang
			f3_plot_y = f3_plot_y / aY_ang
			tot_plot_y = tot_plot_y  / aY_ang
			
			zmin = min(f2_min/max(aX_ang,aY_ang),f3_min/max(aX_ang,aY_ang))
			zmax = max(f2_max/min(aX_ang,aY_ang),f3_max/min(aX_ang,aY_ang))
	
			zmin= -(min(abs(zmin),abs(zmax)))
			zmax= -zmin



	#set uniform colorbar
	print('colormap min_val=',zmin)
	print('colormap max_val=',zmax)
	


	#cb1.set_clim(zmin,zmax)
	a_Rashba = getData('alpha_rashba',searchDir+'/polOutput.txt')



	#plot X shift
	plt.subplot(111)
	CSx_tot = plt.contourf(f2_xi, f2_yi, tot_plot_x, 15, cmap=cmap,vmax=zmax, vmin=zmin)
	if plot_titles:
		if do_plot_percent:
			plt.title(r'$a_x (p_{\mathrm{q}})$',fontsize=title_size)
		else:
			plt.title(r'$a_x$ (Ang)',fontsize=title_size)
	if plot_k_labels:
		plt.ylabel('ky (a.u.)',fontsize=k_label_size)

	plt.show()
	if plot_essin:
		plt.savefig(save_dir+'/a1_n'+str(plotState)+'_Bz'+str(bz)+'aRash'+str(a_Rashba)+'essinShiftX.pdf',bbox_inches='tight')
	else:
		plt.savefig(save_dir+'/a1_n'+str(plotState)+'_Bz'+str(bz)+'aRash'+str(a_Rashba)+'niuShiftX.pdf',bbox_inches='tight')


	#rescale the whole figure to add colorbar
	scale = .8
	plt.tight_layout(pad=1.25,rect=(0,0,scale,scale))

	cax = fig.add_axes([.81, 0.1, .03, scale*0.8])	# [left, bottom, width, height] 
	#cb1 = fig.colorbar(CSy_f2, cax=cax, cmap=cmap, norm=norm, label="1. order shift (ang)")
	norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
	cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
	
	




	#plot Y shift
	plt.subplot(111)
	CSy_tot = plt.contourf(f2_xi, f2_yi, tot_plot_y, 15, cmap=cmap,vmax=zmax, vmin=zmin)
	if plot_titles:
		if do_plot_percent:
			plt.title(r'$a_x (p_{\mathrm{q}})$',fontsize=title_size)
		else:
			plt.title(r'$a_x$ (Ang)',fontsize=title_size)
	if plot_k_labels:
		plt.ylabel('ky (a.u.)',fontsize=k_label_size)

	#rescale the whole figure to add colorbar
	scale = .8
	plt.tight_layout(pad=1.25,rect=(0,0,scale,scale))

	cax = fig.add_axes([.81, 0.1, .03, scale*0.8])	# [left, bottom, width, height] 
	#cb1 = fig.colorbar(CSy_f2, cax=cax, cmap=cmap, norm=norm, label="1. order shift (ang)")
	norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
	cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
	

	plt.show()
	if plot_essin:
		plt.savefig(save_dir+'/a1_n'+str(plotState)+'_Bz'+str(bz)+'aRash'+str(a_Rashba)+'essinShiftY.pdf',bbox_inches='tight')
	else:
		plt.savefig(save_dir+'/a1_n'+str(plotState)+'_Bz'+str(bz)+'aRash'+str(a_Rashba)+'niuShiftY.pdf',bbox_inches='tight')


	
	



	


	return axes


#



def plotDir():
	#test
	searchDir 	= './results'
	nwfs 		= 6



	aUtoAngstrm = 0.52917721092
	aX 			= 20.0 * aUtoAngstrm
	aY			= 10.0 * aUtoAngstrm

	print('unit cell=(',aX,', ',aY,'),(Ang)')

	#
	for plotState in range(1,nwfs+1):
		print('plot state=',plotState)
		ax = plotNiuColor(searchDir,  plotState, aX_ang=aX,aY_ang=aY, plot_percent=True, plot_essin=False, plot_titles=True)


plotDir()