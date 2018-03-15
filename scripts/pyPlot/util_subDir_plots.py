import os
import os.path
import matplotlib.pyplot as plt
import matplotlib as mpl
#my stuff
from plot_bandStruct 	import plotBands
from plot_niu_2dcolor 	import plotNiuColor


#BANDSTRUCTURE PLOT
minEn 			= -50.0
maxEn 			= +50.0
show_Bext_box 	= False
show_Rash_box	= True

#show plots after saving theme
show_plots		= False


#NIU RESPONSE PLOT
nWfs 	= 6

plot_titles		= False
plot_k_labels	= True
plot_descriptor	= False

#sequenital colormaps
cmap 			= mpl.cm.viridis
#cmap 			= mpl.cm.plasma
#cmap 			= mpl.cm.inferno
#cmap 			= mpl.cm.magma

#diverging colormaps
#cmap			= mpl.cm.coolwarm

#well known
#cmap			= mpl.cm.jet



def plot_Dir(dirpath=".",show_Plots=False):
	print('search directory ',dirpath)
	info_path =dirpath+'/polOutput.txt'
	try:
		info_file = open(info_path)
	except IOError:
		print('IOError: did not find file ',info_path)
	else:
		info_file.close()
		#
		#BANDSTRUCTURE
		axBands = plotBands(dirpath, dirpath, minEn, maxEn, show_Bext_box, show_Rash_box)
		if show_Plots:
			plt.show()
		#
		#FIRST ORDER
		for n in range(1,nWfs+1):
			axNiu = plotNiuColor(dirpath,n,cmap,plot_titles=plot_titles, plot_k_labels=plot_k_labels, plot_descriptor=plot_descriptor)
		if show_Plots:
			plt.show()
		print('successfully ploted ',dirpath)
	finally:
		print('done with',dirpath)
		












for dirpath, dirnames, filenames in os.walk("."):
	#create the pdf plot
	if 'pycache' in dirpath:
		print('exclude diretory ',dirpath)
	else:
		plot_Dir(dirpath,show_plots)
	
		
		

		

