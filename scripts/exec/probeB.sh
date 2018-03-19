#!/bin/bash
root=$PWD



#set mpi processes
nProcs=16


#data to calculate
Bz=( -10.0 -5.0 -1.0 0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 25.0 ) #10.0 100.0 1000.0 )




#set numerical paras
gCut=28.0	#set Gcut
nSolve=96
nQx=8
nQy=16
nKx=33
nKy=65
nShell=1

#
nSCx=$((nQx+1))
nSCy=$((nQy+1))


#set numerics
infoString='# changed by probeB.sh'
sed -i  "/Gcut/c\   Gcut        = $gCut 	$infoString"	 	./input.txt
sed -i  "/nSolve/c\   nSolve      = $nSolve 	$infoString"	./input.txt
sed -i 	"/nQx/c\   nQx         = $nQx 	$infoString" 			./input.txt
sed -i 	"/nQy/c\   nQy         = $nQy 	$infoString" 			./input.txt
sed -i 	"/nSCx/c\   nSCx        = $nSCx 	$infoString" 		./input.txt
sed -i 	"/nSCy/c\   nSCy        = $nSCy 	$infoString" 		./input.txt
sed -i 	"/nKx/c\   nKx         = $nKx 	$infoString" 			./input.txt
sed -i 	"/nKy/c\   nKy         = $nKy 	$infoString" 			./input.txt
sed -i 	"/shells/c\   shells       = $nShell 	$infoString" 	./input.txt

#set Hamiltionian
sed -i 	"/doVdesc/c\	doVdesc		= f 	$infoString"		./input.txt
sed -i 	"/doZeeman/c\	doZeeman	= t 	$infoString"		./input.txt
sed -i 	"/doMagHam/c\	doMagHam	= f 	$infoString"		./input.txt
sed -i "/doRashba/c\	doRashba	= f		$infoString"		./input.txt

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
	#set magnetic field amplitude
	cd $dir
	sed -i "/B0/c\    B0  =  $bfield " ./input.txt
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