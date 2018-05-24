import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches



def get_example(use_angstrom=False,use_eV=False):
	#get example data to play around
	aX	= 10
	aY	= 10
	#
	nAt		= 2
	relXpos = np.array(		[0.25,0.75]			) 
	relYpos = np.array(		[0.5,0.5]			) 
	atRx	= np.array(		[1.0,2.0]			)
	atRy	= np.array(		[0.5,3.0]			)
	atPot	= np.array(		[-1.0,-1.5]			)
	dVpot	= np.array(		[.0,.0]				)
	#
	#
	#HANDLE UNITS
	aUtoAngstrm = 0.52917721092
	aUtoEV		= 27.211385
	#
	if use_angstrom:
		aX 		= aUtoAngstrm * aX
		aY 		= aUtoAngstrm * aY
		atRx 	= aUtoAngstrm * atRx
		atRy 	= aUtoAngstrm * atRy
	#
	if use_eV:
		atPot = aUtoEV * atPot
		dVpot = aUtoEV * dVpot
	#
	return aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot



def read_line(line,dtype=float,count=1):
	#extracts all numbers from string
	#
	#extract value(s)
	string 	= line.split('=')[-1]
	raw		= np.fromstring(string,dtype=dtype,count=count,sep=" ")
	#
	#make scalars scalar
	if raw.size == 1:
		raw	= raw[-1]
	#
	return raw


def get_input_data(fpath,use_angstrom=False,use_eV=False):
	#traverse the input file and read the following vars
	nAt 	= 0
	aX		= 0
	aY		= 0
	relXpos = []
	relYpos = []
	atRx	= []
	atRy	= []
	atPot	= []
	dVpot	= []	
	#read file
	with open(fpath,'r') as f:
		for counter,line in enumerate(f):
			#lattice constants
			if 'aX' in line:
				aX	= read_line(line,dtype=float,count=1)
			if 'aY' in line:
				aY	= read_line(line,dtype=float,count=1)
			#
			#atom count
			if 'nAt' in line:
				nAt	= read_line(line,dtype=int,count=1)
			#
			#atom containers
			if nAt > 0:
				if 'relXpos' in line:
					relXpos 	= 			read_line(line,dtype=float,count=nAt)
				if 'relYpos' in line:
					relYpos 	= 			read_line(line,dtype=float,count=nAt)
				if 'atRx' in line:
					atRx 		=  			read_line(line,dtype=float,count=nAt)
				if 'atRy' in line:
					atRy 		=  			read_line(line,dtype=float,count=nAt)
				if 'atPot' in line:
					atPot		=			read_line(line,dtype=float,count=nAt)
				if 'dVpot' in line:
					dVpot		=			read_line(line,dtype=float,count=nAt)
	#
	#
	#HANDLE UNITS
	aUtoAngstrm = 0.52917721092
	aUtoEV		= 27.211385
	#
	if use_angstrom:
		aX 		= aUtoAngstrm * aX
		aY 		= aUtoAngstrm * aY
		atRx 	= aUtoAngstrm * atRx
		atRy 	= aUtoAngstrm * atRy
	#
	if use_eV:
		atPot = aUtoEV * atPot
		dVpot = aUtoEV * dVpot



	return aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot














#
#
#------------------------------------------------------------------------------------------------------------------------------------------------

#options
use_angstrom 	= False
use_eV			= False

fpath='./input.orig'




#read input file
aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot = get_input_data(fpath, use_angstrom, use_eV)

#uncomment for example data
#aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot = get_example(fpath, use_angstrom, use_eV)

#set unit strings 
if use_angstrom:
	length_unit = r' ($\AA$)'
else:
	length_unit = ' (a.u.)'

if use_eV:
	energy_unit = ' (eV)'
else:
	energy_unit = r' ($\mathrm{E}_\mathrm{H}$)'









fig1 = plt.figure()
ax1 = fig1.add_subplot(111, aspect='equal')
ax1.set_xlim([0,aX])
ax1.set_ylim([0,aY])

for at in range(nAt):
	#set rectangle
	xmin 	= relXpos[at]*aX - atRx[at]
	ymin 	= relYpos[at]*aY - atRy[at]
	dx		= atRx[at]*2.0
	dy		= atRy[at]*2.0
	Vpot	= atPot[at]
	#
	#plot rectangle
	ax1.add_patch(		patches.Rectangle( (xmin, ymin),  dx,  dy, fill=True, alpha=.4, 
												facecolor='blue', edgecolor='black',
												label='well #'+str(at),
												linestyle='solid', linewidth=2.0
										)
				)
	#create energy string
	pot_formatted = str(	"{:4.2f}".format(atPot[at]) )
	print( pot_formatted + energy_unit )

	#plot energy string
	ax1.text(xmin+dx*.3, ymin+dy*.4, pot_formatted+'\n'+energy_unit, fontsize=12 )

	#plot center indicator
	plt.plot([aX/2.0], [aY/2.0], marker='+', markersize=5, color="black")

ax1.set_xticks([0.0,aX/2.0,aX])
ax1.set_yticks([0.0,aY/2.0,aY])

plt.tick_params(axis='both', which='major', labelsize=11)
plt.xlabel(r'$x$' + length_unit,fontsize=14)
plt.ylabel(r'$y$' + length_unit,fontsize=14)

#plt.legend()
plt.show()