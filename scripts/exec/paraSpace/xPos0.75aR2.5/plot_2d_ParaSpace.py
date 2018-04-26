import os
import os.path
import re











def getData(dirpath="."):
	p = re.compile(r'\d+\.\d+')  # Compile a pattern to capture float values


	#search for pol files
	for dirpath, dirnames, filenames in os.walk(dirpath):
		print('search in dirpath: ',dirpath)
		


		#GET BERRY
		for filename in [f for f in filenames if f.endswith("polOutput.txt")]:
			filepath = dirpath+'/'+filename
			floats = [float(i) for i in p.findall(dirpath)]  # Convert strings to float
			print('found new file:'+filepath+' x2pos='+str(floats[0])+' aR='+str(floats[1])+'(eV AA)')







#test 

getData()