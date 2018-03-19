import os
import os.path
import shutil
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import animation
#my stuff
from plot_bandStruct 	import plotBands
from plot_niu_2dcolor 	import plotNiuColor



#target directories
band_dir 	= "./bands"
niu_dir 	= "./1stShift"


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

movie_path		= band_dir+"/anim_bands.mp4"
movie_fps		= 1












def save_makedir(dirpath):
	if not os.path.exists(dirpath):
		try:
			os.makedirs(dirpath)
		except OSError:
			print('OSError: could not make ',dirpath)
	else:
		print('dir ',dirpath,' already exists!')
		raw_user = input("do you want to wipe dir? (y/n)")
		if raw_user is 'y':
			shutil.rmtree(dirpath)
			print('removed ',dirpath)
			try:
				os.makedirs(dirpath)
			except OSError:
				print('OSError: could not make ',dirpath)
		else:
			print('WARNING different results may be collected in ',dirpath)





def plot_Dir(dirpath=".",show_Plots=False, save_bands_dir=".",save_niuColor_dir="."):
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
		ax = plotBands(dirpath, dirpath, minEn, maxEn, show_Bext_box, show_Rash_box, save_bands_dir)
		if show_Plots:
			plt.show()
		#
		###FIRST ORDER
		#for n in range(1,nWfs+1):
		#	axNiu = plotNiuColor(dirpath,n,cmap,plot_titles=plot_titles, plot_k_labels=plot_k_labels, plot_descriptor=plot_descriptor,save_dir=save_niuColor_dir)
		#if show_Plots:
		#	plt.show()
		print('successfully ploted ',dirpath)
	finally:
		print('done with',dirpath)



def find_between_r( s, first, last ):
	#find substring in string s, between the first and last string
    try:
        start = s.rindex( first ) + len( first )
        end = s.rindex( last, start )
        return s[start:end]
    except ValueError:
        return ""



def convert_bands_to_png(band_dir,png_quality=500):
	#FIRST SORT THE PDFs
	aRashba		= [] 
	bandNames 	= []
	for dirpath, dirnames, filenames in os.walk(band_dir):
		for filename in [f for f in filenames if f.endswith("bands.pdf")]:
			filepath = dirpath+'/'+filename
			raw_a	 = float(		find_between_r(str(filename),'aRashb','bands.pdf')	)

			bandNames.append(filename)
			aRashba.append(raw_a)	
	s	= sorted(zip(aRashba,bandNames))	
	aRashba, bandNames = map(list,zip(*s))
	
	print('detected aRashba:',aRashba)
	print('found the files: ',bandNames)

	#converted to png
	for idx, filename in enumerate(bandNames):
		inPath	= band_dir+'/'+filename
		outPath = band_dir+'/enPlot'+"%04d" % idx+'.png'

		if os.path.exists(outPath):
			try:
				os.remove(outPath)
				print('removed old ',outPath)
			except IOError:
				print(' could not delete ',outPath)

		command = " convert -density "+str(png_quality)+" "+inPath+" "+outPath
		try:
			os.system(command)
		except OSError:
			print('could not execute ',command)
		else:
			print('created ',outPath)




def animate_pngs(png_folder, movie_path, movie_fps=1):
	if os.path.exists(movie_path):
		try:
			os.remove(movie_path)
			print('removed old ',movie_path)
		except IOError:
			print(' could not delete ',movie_path)

	command =  "ffmpeg -framerate "+str(movie_fps)+" -i "+str(band_dir)+"/enPlot%04d.png -c:v libx264 -r 30 "+str(movie_path)
	try:
		os.system(command)
	except OSError:
		print('could not create ',movie_path)




def search_dirs():
	for dirpath, dirnames, filenames in os.walk("."):
		#create the pdf plot
		if ('pycache' in dirpath) or (band_dir in dirpath) or (niu_dir in dirpath):
			print('exclude diretory ',dirpath)
	
		else:
			plot_Dir(dirpath,show_plots, save_bands_dir=band_dir, save_niuColor_dir=niu_dir)	


	
save_makedir(band_dir)
save_makedir(niu_dir)

do_search = input("do you want to search for new data? (y/n)")
if do_search is 'y':
	print('start search dirs')
	search_dirs()


print('*')
print('*')
print('*')
print('*')
print('*')
print('*')
print('now try to convert to pdfs to png')
convert_bands_to_png(band_dir)
print('*')
print('*')
print('*')
print('*')
print('*')
print('*')
print('now try to animate pngs')
animate_pngs(band_dir, movie_path)

