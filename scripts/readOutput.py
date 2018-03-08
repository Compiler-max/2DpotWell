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

for dirpath, dirnames, filenames in os.walk("."):
	print('search directory:',dirpath)
	for filename in [f for f in filenames if f.endswith("polOutput.txt")]:
		print('found new file:')
		print('')

		#print os.path.join(dirpath, filename)
		filepath = dirpath+'/'+filename
		atPot.append(	getData(filepath,'atPot')			)
		Bfield.append(	getData(filepath,'magnetic_field')	)
		p0.append(		getData(filepath,'zero_order')		)
		pf2.append(		getData(filepath,'niu_f2')			)
		pf3.append(		getData(filepath,'niu_f3')			)
		
		#print('atPot	=',	atPot[-1],		potUnit)
		#print('Bfield	=',	Bfield[-1],		BUnit)
		#print('p_0	=',		p0[-1],			polUnit)
		#print('p_f2	=',		pf2[-1],		polUnit)
		#print('p_f3	=',		pf3[-1],		polUnit)
		#print('')
		#print('+++')
		#print('+++')
		#print('+++')


print('found ', len(atPot), ' data file(s)')


#print('Bfield=',getData(filename,'magnetic_field'),BUnit)
#print('p_0 =',getData(filename,'zero_order'),		polUnit)
#print('p_f2=',getData(filename,'niu_f2'),			polUnit)
#print('p_f3=',getData(filename,'niu_f3'),			polUnit)
















#convert list of np.arrays into 2D np array
Bfield 	= np.array(Bfield)
p0		= np.array(p0)
pf2		= np.array(pf2)
pf3		= np.array(pf3)




#fit
fitfunc = lambda p, x: p[0]*x**2 + p[1]*x + p[2] # Target function

errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function








#plot
fig, ax = plt.subplots(1,1)

B_z		= Bfield[:,2]
p0_x 	= p0[:,0]
pf2_x	= pf2[:,0]
pf3_x	= pf3[:,0]

print('x-axis  B_z: ',B_z)
print('y-axis p0_x: ',p0_x)
print('niu f2     :',pf2_x)
print('niu f3     :',pf3_x)

# Initial guess for the parameters & fit
pGuess = [1.0, 0.1, min(p0_x)]
p0fit, success = optimize.leastsq(errfunc, pGuess[:], args=(B_z, p0_x))


#plot zero order
ax.plot(B_z, p0_x, '+', color='red'	)

#plot parabula interpolation
Blin = np.linspace(B_z.min(),B_z.max(),100)
ax.plot(Blin, fitfunc(p0fit,Blin))

#plot slops
dB = .25
for npt, Bpoint in enumerate(B_z):
	Bvecin = np.linspace(Bpoint-dB,Bpoint+dB,100)
	slope =  (pf2_x[npt]+pf3_x[npt])*(Bvecin-dB) + p0_x[npt] 

	ax.plot(Bvecin, slope, color='green')




plt.ylabel('pol')
plt.xlabel('B (T)')

plt.show()