import os
import os.path
import matplotlib.pyplot as plt
from myBandPlot import plotBands
from niuPlot import plotNiuColor

#define energy window to be plotted
minEn = -50.0
maxEn = +50.0

nWfs 	= 6

for dirpath, dirnames, filenames in os.walk("."):
	print('search directory:',dirpath)
	#create the pdf plots
	if dirpath is not ".":
		ax = plotBands(dirpath,dirpath, minEn, maxEn)
		plt.show()
		for n in range(1,nWfs+1):
			plotNiuColor(dirpath,"f2response.txt",n)
			plotNiuColor(dirpath,"f3response.txt",n)

	#look for them
for dirpat, dirnames, filenames in os.walk("."):
	for filename in [f for f in filenames if f.endswith(".pdf")]:
		print('found file in '+dirpath)
		#plotBands()


