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
		for n in range(1,nWfs+1):
			axNiu = plotNiuColor(dirpath,n,cmap,plot_titles=plot_titles, plot_k_labels=plot_k_labels, plot_descriptor=plot_descriptor,save_dir=save_niuColor_dir)
		if show_Plots:
			plt.show()
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

def sort_files(dir, begin_identifier, end_identifier, file_end_identifier):
	#sorts with respect to parameter between begin/end identifier in filename
	parameter	= []
	fileList	= []

	for dirpath, dirnames, filenames in os.walk(dir):
		for filename in [f for f in filenames if f.endswith(file_end_identifier)]:
			filepath	 = dirpath+'/'+filename
			raw_para	 = float(		find_between_r(filename,begin_identifier,end_identifier)	)

			fileList.append(filename)
			parameter.append(raw_para)	
	s	= sorted(zip(parameter,fileList))	
	parameter, fileList = map(list,zip(*s))
	
	return parameter, fileList

def animate_pngs(png_folder, identifier, movie_path, movie_fps=1):
	if os.path.exists(movie_path):
		try:
			os.remove(movie_path)
			print('removed old ',movie_path)
		except IOError:
			print(' could not delete ',movie_path)

	command =  "ffmpeg -framerate "+str(movie_fps)+" -i "+str(png_folder)+identifier+"%04d.png -c:v libx264 -r 30 "+str(movie_path)
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


def pdf_to_png(dir, fileList, out_file_name, png_quality=500):
	for idx, filename in enumerate(fileList):
		inPath	= dir+'/'+filename
		outPath = dir+'/'+out_file_name+"%04d" % idx+'.png'

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


def animate_bands(band_dir, identifier, png_quality=500):
	#FIRST SORT THE PDFs
	aRashba, bandNames = sort_files(band_dir, 'aRashb','bands','bands.pdf')
	print('detected aRashba:',aRashba)
	print('found the files: ',bandNames)

	#converted to png
	pdf_to_png(band_dir, bandNames, identifier)
	
	#create movie
	animate_pngs(band_dir,identifier,band_dir+'/band_anim.mp4')

	#remove pngs
	for fname in bandNames:
		try:
			os.remove(filename)
		except OSError:
			print('OSError: could not remove file ',filename)

def animate_niu(niu_dir,identifier, png_quality=500):
		#FIRST SORT THE PDFs

		#convert to png
		#pdf_to_png(niu_dir,niuNames,'niuPlot')

		for n in range(1,nWfs+1):
			#make sub-directory
			curr_dir	= niu_dir+'/n'+str(n)
			save_makedir(curr_dir)
			print('state n='+str(n))

			#move to sub-dir
			for dirpath, dirnames, filenames in os.walk(niu_dir):
				for filename in [f for f in filenames if f.endswith('niuShift.pdf')]:
					filepath	 = dirpath+'/'+filename
					if '_n'+str(n)+'_' in filename:				
						os.rename(filepath, curr_dir+'/'+filename)

			#sort current sub-dir and make png
			aRashba, niuNames = sort_files(curr_dir,'aRash','niuShift','niuShift.pdf')
			pdf_to_png(curr_dir,niuNames,'n'+str(n)+identifier)
			
			#animate sub-dir
			movie_path = niu_dir+'/niuShift_n'+str(n)+'.mp4'
			animate_pngs(curr_dir+'/n'+str(n),identifier,niu_dir+'/n'+str(n)+'_anim.mp4')

			#remove png files
			for dirpath, dirnames, filenames in os.walk(curr_dir):
				for filename in [f for f in filenames if f.endswith('.png')]:
					try:
						os.remove(filename)
					except OSError:
						print('OSError: could not remove file ',filename)












def main():	
	print('*')
	print('*')
	print('********************PREPARATIONS***********************************')
	#create target directories
	save_makedir(band_dir)
	save_makedir(niu_dir)
	
	#search for new data
	print('*')
	print('*')
	print('*')
	print('*')
	print('*')
	print('********************GREP DATA***********************************')
	do_search = input("do you want to search for new data? (y/n)")
	if do_search is 'y':
		print('start search dirs')
		search_dirs()

	#png filenames
	png_ident_bands	= 'enPlot'
	png_ident_niu	= 'niuShift'

	#make the animiations
	print('*')
	print('*')
	print('*')
	print('*')
	print('*')
	print('********************ANIMATIONS***********************************')
	print('*')
	
	do_search = input("do you want animate the bands? (y/n)")
	if do_search is 'y':
		animate_bands(				band_dir,png_ident_bands)
	
	print('*')
	print('*')
	print('*')
	print('*')

	do_search = input("do you want animate the first order shifts? (y/n)")
	if do_search is 'y':
		animate_niu(		niu_dir, png_ident_niu)
	print('*')
	print('*')
	print('*')
	print('*')
	print('*')
	print('*')
	



main()
