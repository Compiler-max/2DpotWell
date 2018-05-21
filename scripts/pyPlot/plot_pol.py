import sys
import numpy as np
import matplotlib.pyplot as plt

from util_2dPW_Interf import get_All_subDirs

#GET DESCRITPOR
if len(sys.argv)==1:
	print('please give a descriptor to read from polOutput.txt files')
	sys.exit()


raw_descriptor 	= sys.argv[1]
descriptor		= raw_descriptor.strip()	#remove whitspace
print("will search for descriptor: '",str(descriptor[0]),"'.")


#read pol files
data, p0, pf2, pf3, p0_interp, pf2_interp, pf3_interp = get_All_subDirs(descriptor,".")

print(descriptor," :",data)
print('#pf2:',pf2)
print('#pf3:',pf3)





#p0_abs = np.sqrt(		p0[:,0]**2 + p0[:,1]**2 + p0[:,2]**2 )

#p0_max = []
#tmp_max = np.amax(p0_abs)
#for idx, p0 in enumerate(p0_abs):
#	p0_max.append(tmp_max)

#print('p0_max:',p0_max)
Bz = 1.0

pf2 	= Bz * pf2
pf3 	= Bz * pf3

#p0_x 	= p0[:,0]
#p0_table =  np.concatenate(data,p0_x)
#print(p0_table)

aUtoAngstrm = 0.52917721092
polQuantum = 1.144295347539959e-6
aX = 20.0 * aUtoAngstrm
aY = 10.0  * aUtoAngstrm

#PLOT
direction = ['x-pol', 'y-pol']
for ind, string in enumerate(direction):
	fig, ax  = plt.subplots(1,1) 

	a_ind = aX
	if ind is 0:
		a_ind = aX
	elif ind is 1:
		a_ind = aY



	#ax.set_xlim(min(data),max(data))
	#ax.set_xlim(0.0,2.0)

	plt.plot(data,  p0[:,ind],	marker='+',		color='black',	label='p0'		)
	#plt.plot(data,	p0_max[:],		color='black',	label='p0_abs')
	#plt.plot(data,	-1.0*p0_max[:],		color='black',	label='p0_abs')
	#plt.plot(data, pf2[:,ind],	marker='+',		color='red',	label='p_f2'	)
	#plt.plot(data, pf3[:,ind],	marker='+',		color='green',	label='p_f3'	)
	#plt.plot(data, (pf2[:,ind]+pf3[:,ind]) / a_ind, marker='+', color='blue', label='niu')
	#plt.plot(data,-(pf2[:,ind]+pf3[:,ind]) / a_ind, marker='+', color='orange', label='essin')


	#plt.plot(data,pf2[:,ind]-pf3[:,ind], marker='+', color='orange', label='f2-f3')

	#plt.plot(data, pf2_interp[:,ind],	 marker='+',	color='orange',	label='p_f2(I)')
	#plt.plot(data, pf3_interp[:,ind],	 marker='+',	color='lightgreen',	label='p_f2(I)')

	plt.title('OMP: '+string+' Bz='+str(Bz)+' T')
	plt.ylabel(r'electric pol. [$p_Q$]')
	plt.xlabel(descriptor)

	

	plt.legend()
	plt.show()

