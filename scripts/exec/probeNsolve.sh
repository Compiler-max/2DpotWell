#!/bin/bash
root=$PWD



#set mpi processes
nProcs=16


#data to calculate
nSolve=( 6 12 24 48 96 192 384 786	)


#set numerical paras
gCut=28.0	#set Gcut
nSolve=96
bfield=1.0
nQx=8
nQy=16
nShell=1
aRashba=2.0

#
nSCx=$((nQx+1))
nSCy=$((nQy+1))
nKx=$((nQx+1))
nKy=$((nQy+1))

#set numerics
infoString='# changed by probeB.sh'
sed -i  "/Gcut/c\   Gcut        = $gCut 	$infoString"	 	./input.txt
sed -i 	"/nQx/c\   nQx         = $nQx 	$infoString" 			./input.txt
sed -i 	"/nQy/c\   nQy         = $nQy 	$infoString" 			./input.txt
sed -i 	"/nSCx/c\   nSCx        = $nSCx 	$infoString" 		./input.txt
sed -i 	"/nSCy/c\   nSCy        = $nSCy 	$infoString" 		./input.txt
sed -i 	"/nKx/c\   nKx         = $nKx 	$infoString" 			./input.txt
sed -i 	"/nKy/c\   nKy         = $nKy 	$infoString" 			./input.txt
sed -i 	"/shells/c\   shells       = $nShell 	$infoString" 	./input.txt

#set Hamiltionian
sed -i 	"/doVdesc/c\	doVdesc		= f 	$infoString"		./input.txt
sed -i 	"/doZeeman/c\	doZeeman	= f 	$infoString"		./input.txt
sed -i 	"/doMagHam/c\	doMagHam	= f 	$infoString"		./input.txt
sed -i "/doRashba/c\    doRashba    = t     $infoString"        ./input.txt


#folder to save data
label='nSolve'
mkdir -p $label
now="$(date +"%T")"
echo '['$now']: created directory '$label
printf "*\n*\n"







#body
for nState in ${nSolve[*]}; do
	#
	echo '['$(date +"%T")']: start nSolve='$nState
	#	
	#creates directories
	dir=$label'/n'$nState
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
	#sed parameter
	cd $dir
	sed -i "/nSolve/c\    nSolve  =  $nState " ./input.txt
	wait
	#	
	#execute calculation
	./subStream.sh $nProcs
	wait
	echo '['$(date +"%T")']: finished subStream'
	#
	#clean up the folder
	rm ./main.exe
	rm ./subStream.sh
	rm -r w90files
	rm -r output
	rm -r oldInput
	rm -r rawData
	rm ./input.txt
	wait
	cp results/* ./
	wait
	rm -r results


	#jump back to root
	cd $root
	echo '['$(date +"%T")']: done current calc'
	printf '*\n*\n'
done

echo '['$(date +"%T")']: all done, by'