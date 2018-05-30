import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import griddata
from util_2dPW_Interf import read_Niu_shift
from util_2dPW_Interf import getData


#Perceptually Uniform Sequential colormaps








def plotNiu_tot(searchDir, nWfs, 
						aX_ang=0,aY_ang=0, plot_percent=False,
						plot_essin=True,
						cmap=mpl.cm.viridis,
						plot_titles=False, title_size=12, 
						plot_k_labels=True, k_label_size=10,
						plot_descriptor=False, descriptor_size=14,
						axis_tick_size=12,
						show_rashba_box=False,
						save_dir="."):
	f2_path	= searchDir+'/f2response.txt'
	f3_path = searchDir+'/f3response.txt'






	for plotState in range(1, nWfs+1):
		do_f2, f2_xi, f2_yi, f2_plot_x, f2_plot_y, f2_min, f2_max,f2_bz, raw_f2_kx, raw_f2_ky, raw_f2_x, raw_f2_y	=	read_Niu_shift(plotState, f2_path)
		do_f3, f3_xi, f3_yi, f3_plot_x, f3_plot_y, f3_min, f3_max,f3_bz, raw_f3_kx, raw_f3_ky, raw_f3_x, raw_f3_y	= 	read_Niu_shift(plotState, f3_path)

		if plotState is 1:
			raw_x_shift = raw_f2_x + raw_f3_x
			raw_y_shift = raw_f2_y + raw_f3_y
		else:
			raw_x_shift = raw_x_shift + raw_f2_x + raw_f3_x
			raw_y_shift = raw_y_shift + raw_f2_y + raw_f3_y

		print('read data of state n=',plotState)



	
	if abs(f2_bz-f3_bz) > 1e-8:
		print(searchDir+'/: WARNING different fields detected f2='+str(f2_bz)+' f3='+str(f3_bz))
	bz = f2_bz




	tot_plot_x = griddata((raw_f2_kx,raw_f2_ky), raw_x_shift, (f2_xi[None,:], f2_yi[:,None]), method='cubic')
	tot_plot_y = griddata((raw_f2_kx,raw_f2_ky), raw_y_shift, (f2_xi[None,:], f2_yi[:,None]), method='cubic')
	print('finished setting up the 2d griddata')
	
	
	

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
		print('Will use the percentage plot')		
	
	
	#SET RASHBA INFO BOX
	a_Rashba = getData('alpha_rashba',searchDir+'/polOutput.txt')
	#
	box_string = r'$\alpha_{R}=$'+str(a_Rashba)
	box_font_size = 14
	box_x_rel = 0.80
	box_y_rel = 0.85
	box_opacity = 0.6
	


	zmin = -160
	zmax = 160





	##
	#plot X shift
	#plt.subplot(111)
	fig, axes = plt.subplots(1,1)
	CSx_tot = plt.contourf(f2_xi, f2_yi, tot_plot_x, 15, cmap=cmap)#,vmax=zmax, vmin=zmin)
	if plot_titles:
		if do_plot_percent:
			plt.title(r'$a_x (p_{\mathrm{q}})$',fontsize=title_size)
		else:
			plt.title(r'$a_x$ (Ang)',fontsize=title_size)
	if plot_k_labels:
		plt.ylabel('ky (a.u.)',fontsize=k_label_size)
		plt.xlabel('kx (a.u.)',fontsize=k_label_size)
	#
	#rescale the whole figure to add colorbar
	scale = .8
	plt.tight_layout(pad=1.25,rect=(0,0,scale,scale))
	#
	cax = fig.add_axes([.81, 0.1, .03, scale*0.8])	# [left, bottom, width, height] 
	cb1 = fig.colorbar(CSx_tot, cax=cax, cmap=cmap, label=r'polarizability $a_x$ ($\mu p_Q$)')
	cb1.ax.tick_params(labelsize=k_label_size)
	#
	if show_rashba_box:
		axes.text(box_x_rel, box_y_rel, box_string, 
			horizontalalignment='center',
			verticalalignment='center',
			transform = axes.transAxes,
			fontsize=box_font_size, bbox={'facecolor':'white', 'alpha':box_opacity, 'pad':10})
	
	axes.tick_params(axis = 'both', which = 'major', labelsize = axis_tick_size)
	#
	if plot_essin:
		plt.savefig(save_dir+'/a1_n'+str(plotState)+'_Bz'+str(bz)+'aRash'+str(a_Rashba)+'essinShiftX.pdf',bbox_inches='tight')
	else:
		plt.savefig(save_dir+'/a1_n'+str(plotState)+'_Bz'+str(bz)+'aRash'+str(a_Rashba)+'niuShiftX.pdf',bbox_inches='tight')
	plt.show()
	plt.close()
	




	#plot Y shift
	#plt.subplot(111)
	fig, axes = plt.subplots(1,1)
	CSy_tot = plt.contourf(f2_xi, f2_yi, tot_plot_y*1e6, 15, cmap=cmap,vmax=zmax, vmin=zmin)
	if plot_titles:
		if do_plot_percent:
			plt.title(r'$a_y $',fontsize=title_size)
		else:
			plt.title(r'$a_y$ (Ang)',fontsize=title_size)
	if plot_k_labels:
		plt.ylabel(r'$k_y$ (a.u.)',fontsize=k_label_size)
		plt.xlabel(r'$k_x$ (a.u.)',fontsize=k_label_size)
	#rescale the whole figure to add colorbar
	scale = .8
	plt.tight_layout(pad=1.25,rect=(0,0,scale,scale))
	cax = fig.add_axes([.81, 0.1, .03, scale*0.8])	# [left, bottom, width, height] 
	cb1 = fig.colorbar(CSy_tot, cax=cax, cmap=cmap, label=r'polarizability $a_y$ ($\mu p_Q$)')
	cb1.ax.tick_params(labelsize=k_label_size)
	#norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
	#cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
    #                           norm=norm,
    #                           orientation='vertical')

	
	if show_rashba_box:
		axes.text(box_x_rel, box_y_rel, box_string, 
			horizontalalignment='center',
			verticalalignment='center',
			transform = axes.transAxes,
			fontsize=box_font_size, bbox={'facecolor':'white', 'alpha':box_opacity, 'pad':10})


	axes.tick_params(axis = 'both', which = 'major', labelsize = axis_tick_size)

	if plot_essin:
		plt.savefig(save_dir+'/a1'+'_Bz'+str(bz)+'aRash'+str(a_Rashba)+'essinShiftY.pdf',bbox_inches='tight')
	else:
		plt.savefig(save_dir+'/a1'+'_Bz'+str(bz)+'aRash'+str(a_Rashba)+'niuShiftY.pdf',bbox_inches='tight')

	plt.show()
	plt.close()

	print('a_x=',sum(raw_x_shift)/aX_ang,r'\,$(p_Q)$')		
	print('a_y=',sum(raw_y_shift)/aY_ang,r'\,$(p_Q)$')



	
	



	




#



def plotDir():
	#test
	searchDir 	= '.'
	nWfs 		= 2



	aUtoAngstrm = 0.52917721092
	aX 			= 8.0 * aUtoAngstrm
	aY			= 8.0 * aUtoAngstrm


	print('unit cell=(',aX,', ',aY,'),(Ang)')

	#
	plotNiu_tot(searchDir,  nWfs, aX_ang=aX,aY_ang=aY, plot_percent=True, plot_essin=True, plot_titles=False,plot_k_labels=True,k_label_size=14, show_rashba_box=True)


plotDir()