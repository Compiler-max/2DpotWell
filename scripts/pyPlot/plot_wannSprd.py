#scripts that plots the total spread over wannierization cycle
#
# with ability to plot init guess vs random phases
#
#
import matplotlib.pyplot as plt
import re















def get_w90_spreads(fpath="./wf1.wout"):
	start_id 	= 'Initial State'
	cycle_id 	= 'Cycle:'
	sprd_id 	= 'Sum of centres and spreads'
	final_id 	= 'Final State'

	p = re.compile(r'\d+\.\d+')  # Compile a pattern to capture float values
	trigger_start = False

	cycles	= []
	spreads = []

	#read from .wout file
	with open(fpath,'r') as f:
		#traverse lines in file
		for counter, line in enumerate(f):
			if start_id in line:
				trigger_start = True
				start_line	= counter
				cycles.append(0)

			if cycle_id in line:
				ints 	= [int(s) for s in line.split() if s.isdigit()]	# get all integers from line
				cycles.append(ints[0])

			if sprd_id in line:
				floats = [float(i) for i in p.findall(line)]  # Convert strings to float
				spreads.append(floats[3])

			if final_id in line:
				break

	if len(cycles)>len(spreads):
		print('warning did not catch all spreads')
	if len(cycles)<len(spreads):
		print('too many spreads found') 

	#print data
	print('cycles :',cycles)
	print('spreads:',spreads)

	return cycles, spreads



def plot_w90_spreads(fpath_proj="./wf1.wout", fpath_rand="./wf1.wout"):
	#
	# data with initial projection
	cycles_proj, spreads_proj = get_w90_spreads(fpath_proj)
	#
	#data with random phases
	cycles_rand, spreads_rand = get_w90_spreads(fpath_rand)

	#plot figure
	fig, ax  = plt.subplots(1,1) 
	plt.semilogy(cycles_proj,spreads_proj,'.-', label='projection')
	plt.semilogy(cycles_rand,spreads_rand,'.-',	label='random ')

	#figure labels
	plt.title('convergence of Wannierization')
	plt.xlabel('iterations')
	plt.ylabel(r'total spread $\Omega$')

	#plt.ylim([0.1,5.0])

	#finalize
	plt.legend(title='init guess')
	plt.tight_layout()
	plt.show()



plot_w90_spreads("./wf1.wout","./wf1.wout")