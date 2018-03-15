import os
import os.path
import matplotlib.pyplot as plt
import matplotlib as mpl
from myBandPlot 		import plotBands
from potwellInterface 	import print_Info
from niuPlot 			import plotNiuColor

#BANDSTRUCTURE PLOT
minEn 			= -50.0
maxEn 			= +50.0
show_Bext_box 	= False
show_Rash_box	= True

#show plots after saving theme
show_plots		= True


#NIU RESPONSE PLOT
nWfs 	= 6
cmap 			= mpl.cm.viridis
#cmap 			= mpl.cm.plasma
#cmap 			= mpl.cm.inferno
#cmap 			= mpl.cm.magma

plot_titles		= False
plot_k_labels	= True
plot_descriptor	= False
#diverging colormaps
#cmap			= mpl.cm.coolwarm

#well known
#cmap			= mpl.cm.jet



def make_Plots(dirpath,show_Plots=False):
	print('**********system info:********************')
	print_Info(dirpath+'/polOutput.txt')

	print('***********band plot**********************')
	ax = plotBands(dirpath, dirpath, minEn, maxEn, show_Bext_box, show_Rash_box)
	if show_Plots:
		plt.show()

	print('***********first Order plot*****************')
	for n in range(1,nWfs+1):
		ax = plotNiuColor(dirpath,n,cmap,plot_titles=plot_titles, plot_k_labels=plot_k_labels, plot_descriptor=plot_descriptor)
	if show_Plots:
		plt.show()












for dirpath, dirnames, filenames in os.walk("."):
	#create the pdf plot
	if 'pycache' in dirpath:
		print('exclude diretory ',dirpath)
	else:
		print('search directory ',dirpath)
		make_Plots(dirpath,show_plots)

		
		

		
		
		

		#for n in range(1,nWfs+1):

	#look for them
#for dirpat, dirnames, filenames in os.walk("."):
#	for filename in [f for f in filenames if f.endswith(".pdf")]:
#		print('found file in '+dirpath)
		#plotBands()



	
