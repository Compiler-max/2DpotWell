import numpy		as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from util_2dPW_Interf import read_UNKs




a0_to_Angstroem	= 0.52917721092







def plot_loop(nx, ny, nz, kpt_idx, nbnd, wvfct):

	#get unit cell
	aX = float(	input('please give aX in atomic units: ') )
	aY = float(	input('please give aY in atomic units: ') )
	#aX = aX * a0_to_Angstroem
	#aY = aY	* a0_to_Angstroem

	#print(aX)
	#print(aY)


	cnt = True
	while cnt:
		print('*')
		print('*')
		print('***********unk plotter*************************')
		#
		#what?
		bnd = int(		input('Which band do you want to plot? ')	)
		if len(wvfct)>1:
			kpt = int(		input('At which kpt do you want to plot?') 	)
			if kpt > kpt_idx[-1]:
				print('kpt out of bounds')
				kpt = 0
			if bnd < 0 or bnd > nbnd:
				print('band out of bound')
				bnd = 0
		else:
			kpt = 0

		print('will plot band #',bnd)
		print('at kpt#',kpt)
		#plot
		unk = wvfct[kpt]
		create_plot(kpt, nx, ny, nz, aX, aY, bnd, unk)
		#
		#continue?
		cnt_string = input('do you want to plot more (y/n)?')
		if cnt_string is not 'y':
			cnt = False
			print('by by')




def create_plot(kpt, nx, ny, nz, aX, aY, bnd, unk):

	#grid points
	x_lin	= np.linspace(0,aX,nx)
	y_lin	= np.linspace(0,aY,ny)
	X, Y = np.meshgrid(x_lin, y_lin)


	print('xlin:',x_lin)


	#project onto 2D plane
	unk_plane = unk[bnd,:,:,0]

	unk_plot  = np.empty([ny,nx])
	print(unk_plot)

	#create absoute value	
	for row_idx, row in enumerate(unk_plane):
		for elem_idx, elem in enumerate(row):
			unk_plot[elem_idx,row_idx,] = elem[0]**2 + elem[1]**2

	unk_plot.reshape([ny,nx])



	print('plot points=',len(unk_plane))
	print(' xgrid = ',len(X))
	print(' ygrid = ',len(Y))


	fig, ax = plt.subplots()
	CS = plt.contourf(X, Y, unk_plot, cmap=mpl.cm.magma)
	
	cbar = fig.colorbar(CS,orientation='vertical')


	plt.xlabel(r'x ($a_0$)')
	plt.ylabel(r'y ($a_0$)')


	plt.title(r'$u_{n='+str(bnd)+', k}$')

	plt.savefig('u_n'+str(bnd+1)+'_k'+str(kpt+1)+'.pdf')
	plt.show()

	#wvfct = griddata((x_lin,y_lin), plot_unk, (x_lin[None,:], y_lin[:,None]), method='cubic')


	#make plot
	#fig, axes 	= plt.subplots(nrows=1, ncols=1)
	#plt.subplot(111)
	#C_plot		= plt.contourf(x_lin, y_lin, wvfct, 15, cmap=cmap,vmax=zmax, vmin=zmin)





	#finalize
	print('created plot')
	return



def main():
	nx, ny, nz, nbnd, kpt_idx, wvfct	= read_UNKs(".")
	print('*')
	print('*')
	print('*')
	print('*')
	print('****data summary:****************')
	print('kpts=',kpt_idx)

	print('nx=',nx)
	print('ny=',ny)
	print('nz=',nz)
	print('nbnd=',nbnd)

	plot_loop(nx, ny, nz, kpt_idx, nbnd, wvfct)


	return



main()