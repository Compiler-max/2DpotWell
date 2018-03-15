import numpy as np
from scipy.interpolate import griddata
import re



def file_exists(fpath):
	try:
		open(fpath,'r')
	except IOError:
		print("IOError: ",fpath," could not be found")
		return False
	return True




def getData(descriptor, fpath="./polOutput.txt"):
	#read scalar, array between lines start and finish
	start	= 'begin '+descriptor
	finish	= 'end '+descriptor

	data = np.array([])
	try:
		open(fpath,'r')
	except IOError:
		print('IOError file ',fpath,' does not exist' )	
		return data

	parser	= False
	found	= False
	with open(fpath,'r') as f:
		for counter,line in enumerate(f):
			if finish in line:
				parser	= False
				found	= True
				break
			if parser:
				data = np.fromstring( line, dtype=np.float, sep=' ' )	
				
				#uncommend for debug output:
				#print('getData: found data assoc. to "',start,'" in line ',counter)
				#print(line)
			if start in line:
				parser = True

	if not found:
		print('WARNING: could not file descriptor: "'+descriptor+'" in file '+fpath)
	if(len(data) == 1):
		return data[0]
	return data



def print_Info(fpath="./polOutput.txt"):
	print('Gcut    =',getData('gCut'			,fpath)										)
	print('MP-Grid =',getData('mp_grid'			,fpath)										)
	print('atPot   =',getData('atPot'			,fpath)					,' (eV)'			)
	print('B_ext   =',getData('magnetic_field'	,fpath)					, ' (T)'			)
	print('a_rashba=',getData('alpha_rashba'	,fpath)					,' (eV Ang)'		)
	print('nSolve  =',getData('nSolve'			,fpath)										)
	print('p_0     =',getData('zero_order'		,fpath)					,r'($\mu C / cm$)'	)
	print('p_niu_F2=',getData('niu_f2'			,fpath)					,r'($\mu C / cm$)'	)
	print('p_niu_F3=',getData('niu_f3'			,fpath)					,r'($\mu C / cm$)'	)

	return


def read_AbIn_energies(fpath="./enABiN.txt"):
   	#
   	#read the abinit band structure
   	#
	#
	exists = False
	try:
		out_file     = open(fpath,'r')
	except IOError:
		print("IOError: ",fpath," could not be found")
		return 0, 0, np.empty([]), np.empty([]), np.array([0,0,0])
	#read data
	data        = np.genfromtxt(out_file, dtype=[('int'),('float64'),('float64'),('float64'),('float64')])
	out_file.close()
	exists = True
	#
	#read data lines
	qpts = []
	en_abi= []
	for counter, line in enumerate(data):
	    en_abi.append( line[4] )
	    if line[0] > len(qpts):
	        qpts.append(    ([line[1],line[2],line[3]])   )   
	#
	#convert to numpy arrays
	qpts    = np.array(qpts)
	nQ      = len(qpts)
	nSolve  = int(  len(en_abi)/ nQ )
	en_abi= np.reshape(en_abi,(nQ,nSolve))
	#
	return exists, nQ, nSolve, qpts, en_abi





def read_w90_energies(fpath="./wf1_geninterp.dat"):
	#read the genInterp.dat file
	#
	exists	= False
	nK		= 0
	nWfs 	= 0
	kpts 	= []
	en_w90	= []
	velo	= []
	try:
		w90_file  = open(fpath,'r')
	except IOError:
		print("IOError: ",fpath," could not be found")
		return exists, nK, nWfs, np.array(kpts), np.reshape(en_w90,(nK,nWfs)), velo
	#
	data        = np.genfromtxt(w90_file, dtype=[('int'),('float64'),('float64'),('float64'),('float64'),('float64'),('float64'),('float64')])
	w90_file.close()
	exits 		= True
	#read lines
	for counter, line in enumerate(data):
	    en_w90.append( line[4] )
	    velo.append( [line[5],line[6],line[7]])
	    if line[0] > len(kpts):
	        kpts.append(    ([line[1],line[2],line[3]])   )   
	#
	kpts    = np.array(kpts)
	nK      = len(kpts)
	nWfs    = int(  len(en_w90)/ nK )
	en_w90   = np.reshape(en_w90,(nK,nWfs))
	#
	return exists, nK, nWfs, kpts, en_w90, velo





def read_Niu_shift(nWf, fpath):
	#
	#	reads the niu first order k-resolved response files
	#
	exists = False
	try:
		open(fpath,'r')
	except IOError:
		print('IOError: could not find file',fpath)
		xi 			= np.empty([])
		yi 			= np.empty([])
		aNiuX_plot 	= np.empty([])
		aNiuY_plot	= np.emtpy([])
		zmin		= 0
		zmax 		= 0
		bz			= 0.0
		return exists, xi,yi, aNiuX_plot, aNiuY_plot, zmin, zmax, bz
	exists = True
	

	#read header
	with open(fpath,'r') as f:
		dummy		= f.readline()
		dummy		= f.readline()
		paraString	= f.readline()
		fieldString	= f.readline()
	
	nWfs, nKpts = np.fromstring(paraString,dtype=int,count=2,sep=" ")
	bz			= np.fromstring(fieldString,dtype=float,count=1,sep=" ")[0]
	
	
	#read body
	outFile     = open(fpath,'r')
	data        = np.genfromtxt(outFile, dtype=[('int'),('int'),('float64'),('float64'),('float64'),('float64'),('float64'),('float64')],skip_header=4)
	
	kptsX	= np.array([])
	kptsY	= np.array([])
	kptsZ	= np.array([])
	aNiu_x 	= np.array([])
	aNiu_y 	= np.array([])
	aNiu_z 	= np.array([])  
	
	for counter, line in enumerate(data):
		stat = line[0]
		kpt = line[1]
	
		if stat == nWf:
			kptsX	= np.append(kptsX,line[2])
			kptsY	= np.append(kptsY,line[3])
			kptsZ	= np.append(kptsZ,line[4])
	
			aNiu_x	= np.append(aNiu_x,line[5])
			aNiu_y	= np.append(aNiu_y,line[6])
			aNiu_z	= np.append(aNiu_z,line[7])
	

	#get min and max of data
	zmin 	= min(aNiu_x)
	tmp 	= min(aNiu_y)
	zmin	= min(zmin,tmp) 
	zmax 	= max(aNiu_x)
	tmp		= max(aNiu_y)
	zmax	= max(zmax,tmp)

	# create x-y points to be used in heatmap
	xi = np.linspace(kptsX.min(),kptsX.max(),len(kptsX))
	yi = np.linspace(kptsY.min(),kptsY.max(),len(kptsY))
	
	# Z is a matrix of x-y values
	aNiuX_plot = griddata((kptsX,kptsY), aNiu_x, (xi[None,:], yi[:,None]), method='cubic')
	aNiuY_plot = griddata((kptsX,kptsY), aNiu_y, (xi[None,:], yi[:,None]), method='cubic')


	return exists, xi,yi, aNiuX_plot, aNiuY_plot, zmin, zmax, bz


