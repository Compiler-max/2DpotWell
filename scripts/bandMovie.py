import os
import os.path
import matplotlib.pyplot as plt
from myBandPlot import plotBands

#define energy window to be plotted
minEn = -50.0
maxEn = +50.0


for dirpath, dirnames, filenames in os.walk("."):
	print('search directory:',dirpath)
	#create the pdf plots
	if dirpath.endswith("results"):
		plotBands(dirpath,dirpath, minEn, maxEn)

	#look for them	
	for filename in [f for f in filenames if f.endswith(".pdf")]:
		print('found file in '+dirpath)
		#plotBands()