#!/bin/bash
root=$PWD



#set mpi processes
nProcs=16


#data to calculate
Bz=( -1.0 0.0 1.0 ) #10.0 100.0 1000.0 )

Vpot=-2.8				#atom 1 = Vpot
Vpot2=( -2.8 -2.9 ) #)-3.0 -3.2 -3.4 )	#atom 2 = Vpot - dV


#set numerical paras
gCut=10.0	#set Gcut
nSolve=96
nQx=4
nShell=1

#
nSCx=$((nQx+1))
nKx=$((nQx+1))


#set numerics (for i1 & i2)
infoString='# changed by probeB.sh'
sed -i  "/Gcut/c\   Gcut        = $gCut $infoString"	 		./input.txt
sed -i  "/nSolve/c\   nSolve      = $nSolve $infoString"	./input.txt
sed -i 	"/nQx/c\   nQx         = $nQx $infoString" 			./input.txt
sed -i 	"/nQy/c\   nQy         = $nQx $infoString" 			./input.txt
sed -i 	"/nSCx/c\   nSCx        = $nSCx $infoString" 		./input.txt
sed -i 	"/nSCy/c\   nSCy        = $nSCx $infoString" 		./input.txt
sed -i 	"/nKx/c\   nKx         = $nKx $infoString" 			./input.txt
sed -i 	"/nKy/c\   nKy         = $nKx $infoString" 			./input.txt
sed -i 	"/shells/c\   shells       = $nShell $infoString" 	./input.txt




#folder to save data
label='Bprobe'
mkdir -p $label
now="$(date +"%T")"
echo '['$now']: created directory '$label
printf "*\n*\n"







#body
for bfield in ${Bz[*]}; do
	for v2 in ${Vpot2[*]}; do
		#
		echo '['$(date +"%T")']: start Bz='$bfield'T; V1='$Vpot'; V2='$v2'; nShell='$nShell
		#	
		#creates directories
		dir=$label'/b'$bfield'dV'$v2
		if [[ ! -e $dir ]]; then
   		 	mkdir $dir
   		 	echo '['$(date +"%T")']: created dir '$dir
		else
    		rm -r -f $dir/*
    		echo '['$(date +"%T")']: '$dir' was wiped'
		fi
		#
		#cpy relevant files to subdirectory
		cp main.exe $dir
		cp input.txt $dir
		cp subStream.sh $dir
		wait
		#
		#prepare input file
		cd $dir
		sed -i "/B0/c\   B0		= $Bz $infoString"	./input.txt
		sed -i "/Bext/c\   Bext	= 0.0 0.0 1.0 $infoString" ./input.txt
		sed -i "/atPot/c\   atPot	= $Vpot $v2 $infoString" ./input.txt
		wait
		#	
		#execute calculation
		./subStream.sh $nProcs
		wait
		echo '['$(date +"%T")']: finished subStream'
		#
		#jump back to root
		cd $root
		echo '['$(date +"%T")']: done current calc'
		printf '*\n*\n'
	done
done

echo '['$(date +"%T")']: all done, by'