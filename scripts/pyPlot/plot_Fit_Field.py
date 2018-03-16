import numpy as np
import os
import os.path
import matplotlib.pyplot as plt
from scipy import optimize
from potwellInterface import getData


filename	= 'polOutput.txt'

BUnit		= 'T'
polUnit		= 'mu C/cm'
potUnit		= 'eV'


#plots
fitCol		= 'blue'
dataCol		= 'red'
interpCol	= 'orange'















#def getData(filename, descriptor):
#	#read array between lines start and finish
#	start	= 'begin '+descriptor
#	finish	= 'end '+descriptor
#
#	parser	= False
#	found	= False
#	data = []
#	with open(filename) as f:
#		for counter,line in enumerate(f):
#			if finish in line:
#				parser	= False
#				found	= True
#				break
#			if parser:
#				data = np.fromstring( line, dtype=np.float, sep=' ' )	
#				
#				#uncommend for debug output:
#				#print('getData: found data assoc. to "',start,'" in line ',counter)
#				#print(line)
#			if start in line:
#				parser = True
#
#	if not found:
#		print('WARNING: could not file descriptor: "'+descriptor+'" in file '+filename)
#	return data
	

#numerics
gCut 		= []
mpGrid		= []
nSolve		= []

#features
atPot 		= []
Bfield 		= []
aRashba		= []

#results
p0			= []
pf2			= []
pf3			= []
			
#interpolation
p0_interp	= []
pf2_interp	= []
pf3_interp 	= []

for dirpath, dirnames, filenames in os.walk("."):
	#GET BERRY
	for filename in [f for f in filenames if f.endswith("polOutput.txt")]:
		filepath = dirpath+'/'+filename
		print('found new file:'+filepath)
		
		gCut.append(		getData('gCut'				,filepath)			)
		mpGrid.append(		getData('mp_grid'			,filepath)			)
		nSolve.append(		getData('nSolve'			,filepath)			)

		atPot.append(		getData('atPot'				,filepath)			)
		Bfield.append(		getData('magnetic_field'	,filepath)			)
		p0.append(			getData('zero_order'		,filepath)			)
		pf2.append(			getData('niu_f2'			,filepath)			)
		pf3.append(			getData('niu_f3'			,filepath)			)
		aRashba.append(		getData('alpha_rashba'		,filepath)			)
		#
	#GET INTERPOLATION
	for filename in [f for f in filenames if f.endswith("polInterp.txt")]:
		filepath = dirpath+'/'+filename
		print('found new interp file:'+filepath)
		
		p0_interp.append(	getData('zero_order_sum',filepath)	)
		pf2_interp.append(	getData('f2_sum',filepath)			)
		pf3_interp.append(	getData('f3_sum',filepath)			)
		


print('found ', len(atPot)		,	' data file(s)')
print('found ', len(p0_interp)	,	' interpolation data file(s)')


print('gCut:'+str(gCut))
print('mpGrid:'+str(mpGrid))


foundInterp = False
if len(p0_interp) > 0:
	if len(p0_interp[0]) > 0:
		foundInterp = True
		print('found interpolation data')














#convert list of np.arrays into 2D np array
gCut		= np.array( gCut			)
Bfield 		= np.array(	Bfield			)
aRashba		= np.array(	aRashba			)
p0			= np.array(	p0				)
pf2			= np.array(	pf2				)
pf3			= np.array(	pf3				)

p0_interp 	= np.array(	p0_interp		)
pf2_interp	= np.array( pf2_interp		)
pf3_interp	= np.array( pf3_interp)




#Get the desired vectorial components
B_z			= Bfield[:,2]
p0_x 		= p0[:,0]
pf2_x		= pf2[:,0]
pf3_x		= pf3[:,0]

if foundInterp:
	p0_interp_x	= p0_interp[:,0]
	pf2_interp_x= pf2_interp[:,0]
	pf3_interp_x= pf3_interp[:,0]




#fit function
fitfunc = lambda p, x: p[0]*x**2 + p[1]*x + p[2] # Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
fitderiv= lambda p, x: 2.0*p[0]*x + p[1]





print('*****data **************************')
print('x   aRashba: ',aRashba)
print('x-axis  B_z: ',B_z)
print('y-axis p0_x: ',p0_x)
print('slope niuf2:',pf2_x)
print('slope niuf3:',pf3_x)
print(' ')

pGuess = [1.0, 0.1, min(p0_x)]
p0fit, success = optimize.leastsq(errfunc, pGuess[:], args=(B_z, p0_x))

if not success:
	print('WARNING: optimizer is stuggling')

print('*****parabula fit*******************')
print(str(p0fit[0])+'	B**2 + '+str(p0fit[1])+'	B +'+str(p0fit[2]))
# Initial guess for the parameters & fit





#plot parabula
fig, ax = plt.subplots(1,1)

#plot zero order
ax.plot(B_z, p0_x, '.', color=dataCol	)
if foundInterp:
	ax.plot(B_z, p0_interp_x, '+',	color=fitCol	)

#plot parabula interpolation
print('interpolate in range B_z=['+str(B_z.min())+':'+str(B_z.max())+'] (T)')
Blin = np.linspace(B_z.min(),B_z.max(),100)
ax.plot(Blin, fitfunc(p0fit,Blin))

#plot slops (niu)
dB = .25
for npt, Bpoint in enumerate(B_z):
	Bvecin = np.linspace(Bpoint-dB,Bpoint+dB,100)
	slope =  (pf2_x[npt]+pf3_x[npt])*(Bvecin-dB) + p0_x[npt] 

	ax.plot(Bvecin, slope, color=dataCol)


#plot limits
ax.set_xlim(B_z.min(),B_z.max())

xlabel = 'B (T)'

plt.ylabel('pol')
plt.xlabel(xlabel)

plt.show()





#slop plot:
pf_x 		= pf2_x + pf3_x
if foundInterp:
	pf_interp_x	= pf2_interp_x + pf3_interp_x


fig, ax 	= plt.subplots(1,1)

#plot fit slope
ax.plot(Blin, fitderiv(p0fit,Blin), color=fitCol	)

#plot niu slopes
ax.plot(B_z, pf_x ,'*',color=dataCol)

if foundInterp:
	ax.plot(B_z, pf_interp_x,'*',color=interpCol)

#plot limits
ax.set_xlim(B_z.min(),B_z.max())


#labels
plt.title('slopes')
plt.xlabel(xlabel)

plt.show()