import numpy as np
import os
import os.path
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
				print('getData: found data assoc. to "',start,'" in line ',counter)
				print(line)
				print('extracted data=',data)
			if start in line:
				parser = True

	if not found:
		print('WARNING: could not find descriptor: "'+descriptor+'" in file '+fpath)
	if(len(data) == 1):
		return data[0]
	return data









def read_w90_out(nAt, nWf, fpath="./wf1.wout"):
	with open(fpath,'r') as f:
		print('found file',fpath)
		foundSite	= False
		foundFinal	= False
		p = re.compile(r'\d+\.\d+')  # Compile a pattern to capture float values
		at_cent = []
		wf_cent = []
		wf_sprd = []


		for counter, line in enumerate(f):
			#read unit cell
			if 'a_1' in line:
				floats = [float(i) for i in p.findall(line)]  # Convert strings to float
				a_X = floats[0]	#get only first component (cubic cell)

			if 'a_2' in line:
				floats = [float(i) for i in p.findall(line)]  # Convert strings to float
				a_Y = floats[1]	#get only second component (cubic cell)

			if 'a_3' in line:
				floats = [float(i) for i in p.findall(line)]  # Convert strings to float
				a_Z = floats[2]	#get only third component (cubic cell)

			
			#read atom cites
			if 'Site' in line:
				foundSite = True
				siteLine = counter
			if foundSite:
				for at in range(nAt):
					if counter == siteLine+2+at:
						floats = [float(i) for i in p.findall(line)]  # Convert strings to float 
						at_cent.append([floats[3],floats[4],0.0])

			#read mp-grid
			if 'Grid size' in line:
				ints 	= [int(s) for s in line.split() if s.isdigit()]	# get all integers from line
				mpGrid 	= ints[0:3]

			#read final wf
			if 'Final State' in line:
				foundFinal = True
				finalState = counter
			if foundFinal:
				for wf in range(nWf):
					if counter == finalState+1+wf:
						#print(line)
						floats = [float(i) for i in p.findall(line)]  # Convert strings to float 
						wf_cent.append(floats[0:3])
						wf_sprd.append(floats[3])

		#convert to np array
		at_cent = np.array(		at_cent	)
		mpGrid	= np.array(		mpGrid	)
		wf_cent = np.array(		wf_cent )
		wf_sprd = np.array(		wf_sprd )

		#summary
		#print('aCell =(',a_X,', ',a_Y,', ',a_Z,').')
		#print('at_cent :',at_cent)
		#print('mp grid:',mpGrid)
		#print('wf_cent:',wf_cent)
		#print('wf_sprd:',wf_sprd)

		#sort wf to atom
		wf_shift = []
		for wf in wf_cent:
			norm = []
			for at in range(nAt):
				norm.append(np.linalg.norm(wf-at_cent[at]))
			atIndex = norm.index(min(norm))
			wf_shift.append(at_cent[atIndex]-wf)
		wf_shift = np.array( wf_shift )
		wf_shift_sum = sum(wf_shift)

		wf_pol = np.array([ wf_shift_sum[0]/a_X, wf_shift_sum[1]/a_Y, 0.0])
		
		wf_sprd_sum = sum(wf_sprd)

		#print('wf_shift sum:',wf_shift_sum)
		print('tot sprd: ',wf_sprd_sum,'; pol =',wf_pol," (pol Quant.)")



	return at_cent, mpGrid, wf_cent, wf_sprd, wf_shift, wf_pol


def read_gcut_probe_w90(nAt, nWfs, dirpath="."):
	gCut		= []
	at_cent		= []
	mpGrid		= []
	wf_cent		= []
	wf_sprd 	= []
	wf_shift	= []
	wf_pol		= []
	p = re.compile(r'\d+\.\d+')  # Compile a pattern to capture float values


	for dirpath, dirnames, filenames in os.walk(dirpath):
		if 'gCut' in dirpath and 'pycache' not in dirpath:
			print('search in dirpath: ',dirpath)
			floats = [float(i) for i in p.findall(dirpath)]  
			gCut.append( floats[0] )
	
			for filename in [f for f in filenames if f.endswith("wf1.wout")]:
				filepath = dirpath+'/'+filename
				print('found new file:'+filepath)
				at_cent_temp, mpGrid_temp, wf_cent_temp, wf_sprd_temp, wf_shift_temp, wf_pol_temp = read_w90_out(nAt,nWfs,filepath)

			at_cent.append(at_cent_temp)
			mpGrid.append(mpGrid_temp)			
			wf_cent.append(wf_cent_temp)
			wf_sprd.append(wf_sprd_temp)
			wf_shift.append(wf_shift_temp)
			wf_pol.append(wf_pol_temp)

	#sort with respect to gcut value
	try:
		s	= sorted(zip(gCut, at_cent, mpGrid, wf_cent, wf_sprd, wf_shift, wf_pol))
		gCut, at_cent, mpGrid, wf_cent, wf_sprd, wf_shift, wf_pol = map(list,zip(*s))
		print('read_gcut_probe_w90: sorted data successfully')
	except:
		print('read_gcut_probe_w90: could not sort data')
	finally:
		print('done reading data')


	return gCut, at_cent, mpGrid, wf_cent, wf_sprd, wf_shift, wf_pol




def get_pw_count(dirpath):
	gCut = []
	pw_count = []
	p = re.compile(r'\d+\.\d+')  # Compile a pattern to capture float values
	print('hello from get_pw_count')

	for dirpath, dirnames, filenames in os.walk(dirpath):
		if 'gCut' in dirpath and 'pycache' not in dirpath:
			print('search in dirpath: ',dirpath)
			floats = [float(i) for i in p.findall(dirpath)]  
			gCut.append( floats[0] )
	
			for filename in [f for f in filenames if f.endswith("polOutput.txt")]:
				pw_count.append(getData('Gmin',dirpath+'/'+filename))


			print('found Gcut=',floats[0],' and associated Gmin=',pw_count[-1])
	
	try:
		s	= sorted(zip(gCut, pw_count))
		gCut, pw_count= map(list,zip(*s))
		print('get_pw_count: sorted the data')
	except:
		print('get_pw_count: could not sort data in path='+dirpath)


	print('gCut=',gCut)
	print('#pw=',pw_count)

	for pw in pw_count:
		pw = int(pw)

	return pw_count




def get_All_subDirs(descriptor,path="."):
	data = []
	#sys
	uCell 	= []
	#numerics
	gCut 	= []
	mpGrid 	= [] 
	#results
	p0			= []
	pf2			= []
	pf3			= []
				
	#interpolation
	p0_interp	= []
	pf2_interp	= []
	pf3_interp 	= []

	#search for pol files
	for dirpath, dirnames, filenames in os.walk(path):
		if 'pycache' not in dirpath:
			print('search in dirpath: ',dirpath)
			#GET BERRY
			for filename in [f for f in filenames if f.endswith("polOutput.txt")]:
				filepath = dirpath+'/'+filename
				print('found new file:'+filepath)
				#
				uCell.append(		getData('unit_cell'			,filepath)			)
				#
				gCut.append(		getData('gCut'				,filepath)			)
				mpGrid.append(		getData('mp_grid'			,filepath)			)
				#
				data.append(		getData( descriptor			,filepath)			)
				p0.append(			getData('sum_zero_shifts'		,filepath)			)
				pf2.append(			getData('sum_niu_f2'			,filepath)			)
				pf3.append(			getData('sum_niu_f3'			,filepath)			)
			#GET INTERPOLATION
			for filename in [f for f in filenames if f.endswith("polInterp.txt")]:
				filepath = dirpath+'/'+filename
				print('found new interp file:'+filepath)
				#
				p0_interp.append(	getData('zero_order_sum',filepath)	)
				pf2_interp.append(	getData('f2_sum',filepath)			)
				pf3_interp.append(	getData('f3_sum',filepath)			)
			
	#sort only works for scalars
	if 'magnetic_field' in descriptor:
		print('magnetic field will be stripped to z-component')
		bz = []
		for bVec in data:
			bz.append(bVec[2])
		print('Bz=',bz)
		data = bz

	#print('data=',data)
	#print('p0=',len(p0))
	#print('pf2=',len(pf2))
	#print('pf3=',len(pf3))
	#print('p0_interp=',len(p0_interp))
	#print('pf2_interp=',len(pf2_interp))
	#print('pf3_interp=',len(pf3_interp))

	#sort the lists in ascending order of data list
	print('..now try to sort the data')
	try:
		s	= sorted(zip(data,p0,pf2,pf3,p0_interp, pf2_interp, pf3_interp))
		data, p0, pf2, pf3, p0_interp, pf2_interp, pf3_interp = map(list,zip(*s))
		print('data sorted successfully')
	except:
		print('could not sort the data')
	finally:
		print('finished sorting the data')


	


	#convert to numpy
	data 		= np.array( data		 	)
	p0			= np.array(	p0				)
	pf2			= np.array(	pf2				)
	pf3			= np.array(	pf3				)
	p0_interp	= np.array(	p0_interp		)
	pf2_interp	= np.array(	pf2_interp		)
	pf3_interp	= np.array(	pf3_interp		)


	return data, p0, pf2, pf3, p0_interp, pf2_interp, pf3_interp



def print_Info(fpath="./polOutput.txt"):
	print('Gcut    =',getData('gCut'			,fpath)										)
	print('MP-Grid =',getData('mp_grid'			,fpath)										)
	print('atPot   =',getData('atPot'			,fpath)					,' (eV)'			)
	print('B_ext   =',getData('magnetic_field'	,fpath)					, ' (T)'			)
	print('a_rashba=',getData('alpha_rashba'	,fpath)					,' (eV Ang)'		)
	print('nSolve  =',getData('nSolve'			,fpath)										)
	print('p_0     =',getData('sum_zero_shifts'		,fpath)				,r'(Ang)'			)
	print('p_niu_F2=',getData('sum_niu_f2'			,fpath)				,r'(Ang)'			)
	print('p_niu_F3=',getData('sum_niu_f3'			,fpath)				,r'(Ang)'			)

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
	exists 		= True
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
		aNiuY_plot	= np.empty([])
		aNiu_x		= np.empty([])
		aNiu_y		= np.empty([])
		kptsX		= np.array([])
		kptsY		= np.array([])
		zmin		= 0
		zmax 		= 0
		bz			= 0.0
		return exists, xi,yi, aNiuX_plot, aNiuY_plot, zmin, zmax, bz,kptsX, kptsY, aNiu_x, aNiu_y
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


	return exists, xi,yi, aNiuX_plot, aNiuY_plot, zmin, zmax, bz, kptsX, kptsY, aNiu_x, aNiu_y





def read_UNKs(directory=".",spin=1):

	kpt_idx = []
	wvfct	= []

	#search for pol files
	for dirpath, dirnames, filenames in os.walk(directory):
		if not ('pycache' in dirpath):
			print('search in dirpath: ',dirpath)
			for filename in [f for f in filenames if f.endswith("."+str(spin))]:
				#
				if 'UNK' in filename:
					print('found file', filename)
					fpath		= dirpath+'/'+filename
					#
					#read the header
					with open(fpath) as fp: 
						firstLine	= fp.readline()
					header	= np.fromstring(firstLine,dtype=int, count=5,sep=" ")
					nx	= header[0]
					ny 	= header[1]
					nz	= header[2]
					ik	= header[3]
					nbnd= header[4]
					#print('nx=',nx)
					#print('ny=',ny)
					#print('nz=',nz)
					#print('ik=',ik)
					#print('nbnd=',nbnd)
					#
					#
					#read the body
					try:
						outFile     = open(fpath,'r')
					except IOError:
						print('IOError: could not open ',fpath)
					else: 
						raw_data        = np.genfromtxt(outFile, dtype=[('float64'),('float64')],skip_header=1)
						outFile.close()
					data	= raw_data.reshape((nbnd,nx,ny,nz))



				kpt_idx.append(ik)
				wvfct.append(data)

	#sort ky k-pt index
	s	= sorted(zip(kpt_idx, wvfct))
	kpt_idx, wvfct = map(list,zip(*s))



	return nx, ny, nz, nbnd, kpt_idx, wvfct