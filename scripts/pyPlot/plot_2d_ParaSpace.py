import os
import os.path
import re











def getData(dirpath="."):
	p = re.compile(r'\d+\.\d+')  # Compile a pattern to capture float values

	descriptor_f2 = 'niu_f2'
	descriptor_f3 = 'niu_f3' 

	niu_f2	= []
	niu_f3 	= []
	x2pos 	= []
	aR		= []



	#search for pol files
	for dirpath, dirnames, filenames in os.walk(dirpath):
		


		#GET BERRY
		for filename in [f for f in filenames if f.endswith("polOutput.txt")]:
			filepath = dirpath+'/'+filename
			print('dirpath=',dirpath)
			floats = [float(i) for i in p.findall(dirpath)]  # Convert strings to float
			print('raw floats',floats)
			if len(floats)>1:
				print('found new file:'+filepath+' x2pos='+str(floats[0])+' aR='+str(floats[1])+'(eV AA)')

				niu_f2.append(		getData(descriptor_f2, fpath=filepath)	)
				niu_f3.append(		getData(descriptor_f3, fpath=filepath)	)


	return x2pos, aR, niu_f2, niu_f3



#test 

x2pos, aR, niu_f2, niu_f3 = getData()