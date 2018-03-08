import numpy as np
import os
import os.path
import matplotlib.pyplot as plt
from scipy import optimize


filename	= 'polOutput.txt'

BUnit		= 'T'
polUnit		= 'mu C/cm'
potUnit		= 'eV'


def getData(filename, descriptor):
	#read array between lines start and finish
	start	= 'begin '+descriptor
	finish	= 'end '+descriptor

	parser	= False
	with open(filename) as f:
		for counter,line in enumerate(f):
			if finish in line:
				parser = False
				break
			if parser:
				data = np.fromstring( line, dtype=np.float, sep=' ' )	
				
				#uncommend for debug output:
				#print('getData: found data assoc. to "',start,'" in line ',counter)
				#print(line)
			if start in line:
				parser = True

	return data
	


atPot 	= []
Bfield 	= []
p0		= []
pf2	= []
pf3	= []			

p0_interp	= []
pf2_interp	= []
pf3_interp 	= []

for dirpath, dirnames, filenames in os.walk("."):
	print('search directory:',dirpath)
	#GET BERRY
	for filename in [f for f in filenames if f.endswith("polOutput.txt")]:
		print('found new file:')
		print('')
		filepath = dirpath+'/'+filename
		
		atPot.append(		getData(filepath,'atPot')			)
		Bfield.append(		getData(filepath,'magnetic_field')	)
		p0.append(			getData(filepath,'zero_order')		)
		pf2.append(			getData(filepath,'niu_f2')			)
		pf3.append(			getData(filepath,'niu_f3')			)
		#
	#GET INTERPOLATION
	for filename in [f for f in filenames if f.endswith("polInterp.txt")]:
		print('found new interp file:')
		print('')
		filepath = dirpath+'/'+filename
		
		p0_interp.append(	getData(filepath,'zero_order_sum')	)
		pf2_interp.append(	getData(filepath,'f2_sum')			)
		pf3_interp.append(	getData(filepath,'f3_sum')			)
		


print('found ', len(atPot)		,	' data file(s)')
print('found ', len(p0_interp)	,	' interpolation data file(s)')
















#convert list of np.arrays into 2D np array
Bfield 	= np.array(Bfield)
p0		= np.array(p0)
pf2		= np.array(pf2)
pf3		= np.array(pf3)




#fit
fitfunc = lambda p, x: p[0]*x**2 + p[1]*x + p[2] # Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
fitderiv= lambda p, x: 2.0*p[0]*x + p[1]





B_z		= Bfield[:,2]
p0_x 	= p0[:,0]
pf2_x	= pf2[:,0]
pf3_x	= pf3[:,0]

print('*****data **************************')
print('x-axis  B_z: ',B_z)
print('y-axis p0_x: ',p0_x)
print('slope niuf2:',pf2_x)
print('slope niuf3:',pf3_x)
print(' ')

pGuess = [1.0, 0.1, min(p0_x)]
p0fit, success = optimize.leastsq(errfunc, pGuess[:], args=(B_z, p0_x))

print('*****parabula fit*******************')
print(str(p0fit[0])+'	B**2 + '+str(p0fit[1])+'	B +'+str(p0fit[2]))
# Initial guess for the parameters & fit





#plot parabula
fig, ax = plt.subplots(1,1)

#plot zero order
ax.plot(B_z, p0_x, '+', color='red'	)

#plot parabula interpolation
print('interpolate in range B_z=['+str(B_z.min())+':'+str(B_z.max())+'] (T)')
Blin = np.linspace(B_z.min(),B_z.max(),100)
ax.plot(Blin, fitfunc(p0fit,Blin))

#plot slops
dB = .25
for npt, Bpoint in enumerate(B_z):
	Bvecin = np.linspace(Bpoint-dB,Bpoint+dB,100)
	slope =  (pf2_x[npt]+pf3_x[npt])*(Bvecin-dB) + p0_x[npt] 

	ax.plot(Bvecin, slope, color='green')


#plot limits
ax.set_xlim(B_z.min(),B_z.max())

xlabel = 'B (T)'

plt.ylabel('pol')
plt.xlabel(xlabel)

plt.show()





#slop plot:
pf_x = pf2_x + pf3_x
fig, ax = plt.subplots(1,1)

#plot fit slope
ax.plot(Blin, fitderiv(p0fit,Blin), color='blue'	)

#plot niu slopes
ax.plot(B_z, pf_x ,'*',color='red')

#plot limits
ax.set_xlim(B_z.min(),B_z.max())


#labels
plt.title('slopes')
plt.xlabel(xlabel)

plt.show()