import numpy as np



def getData(filename, descriptor):
	#read array between lines start and finish
	start	= 'begin '+descriptor
	finish	= 'end '+descriptor

	parser	= False
	found	= False
	data = []
	with open(filename) as f:
		for counter,line in enumerate(f):
			if finish in line:
				parser	= False
				found	= True
				break
			if parser:
				data = np.fromstring( line, dtype=np.float, sep=' ' )	
				
				#uncommend for debug output:
				#print('getData: found data assoc. to "',start,'" in line ',counter)
				#print(line)
			if start in line:
				parser = True

	if not found:
		print('WARNING: could not file descriptor: "'+descriptor+'" in file '+filename)
	return data
	