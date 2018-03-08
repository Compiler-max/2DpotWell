#!/bin/bash
root=$PWD



#set mpi processes
nProcs=16


#data to calculate
Bz=( -1.0 0.0 1.0 ) #10.0 100.0 1000.0 )




#set numerical paras
gCut=28.0	#set Gcut
nSolve=96
nQx=16
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
		#
		echo '['$(date +"%T")']: start Bz='$bfield'T; V1='$Vpot'; V2='$v2'; nShell='$nShell
		#	
		#creates directories
		dir=$label'/b'$bfield
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
		sed -i "/B0/c\   B0		= $bfield $infoString"	./input.txt	
		wait
		#	
		#execute calculation
		./subStream.sh $nProcs
		wait
		echo '['$(date +"%T")']: finished subStream'
		#
		rm ./main.exe
		rm ./subStream.sh
		rm -r w90files
		rm -r output
		rm -r oldInput
		rm -r rawData


		#jump back to root
		cd $root
		echo '['$(date +"%T")']: done current calc'
		printf '*\n*\n'
done

echo '['$(date +"%T")']: all done, by'