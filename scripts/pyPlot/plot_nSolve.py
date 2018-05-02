import numpy as np
import matplotlib.pyplot as plt

from util_2dPW_Interf import get_All_subDirs

descriptor = 'nSolve'

#read pol files
data, p0, pf2, pf3, p0_interp, pf2_interp, pf3_interp = get_All_subDirs(descriptor,".")

print(descriptor," :",data)
print('#pf2:',pf2)
print('#pf3:',pf3)



aUtoAngstrm 	= 0.52917721092
polQuantum = 1.144295347539959e-6


aX = 10.0 *aUtoAngstrm
aY = 5.0 *aUtoAngstrm


p0_abs = np.sqrt(		p0[:,0]**2 + p0[:,1]**2 + p0[:,2]**2 )

p0_max = []
tmp_max = np.amax(p0_abs)
for idx, p0 in enumerate(p0_abs):
	p0_max.append(tmp_max)

print('p0_max:',p0_max)
Bz = 1.0

pf2 	= Bz * pf2
pf3 	= Bz * pf3

#p0_x 	= p0[:,0]
#p0_table =  np.concatenate(data,p0_x)
#print(p0_table)


#PLOT
fig, ax  = plt.subplots(1,1) 

plt.plot(data,(-pf2[:,1]-pf3[:,1]) * 1e6 / (polQuantum*aY), marker='+', color='orange', label='essin')


my_xticks = [6,12,24,48,96]

ax.set_xticks(my_xticks)
	
plt.title('convergence of interband mixing')
plt.ylabel(r'response $p^{(1)}_y  \: (\mu \: p_Q)$')
plt.xlabel('included bands')
#
#plt.legend()
plt.tight_layout()
plt.show()

