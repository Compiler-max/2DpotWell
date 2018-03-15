#!/bin/bash
root=$PWD



#set mpi processes
nProcs=16


#data to calculate
Gvalues=( -2.0 0.0 1.0 1.5 2.0 2.5 3.0) #10.0 100.0 1000.0 )




#set numerical paras
nSolve=96
#nQx=16
nShell=1

#
nSCx=$((nQx+1))
nKx=$((nQx+1))


#set numerics
infoString='# changed by probeB.sh'
sed -i  "/nSolve/c\   nSolve      = $nSolve 	$infoString"	./input.txt
#sed -i 	"/nQx/c\   nQx         = $nQx 	$infoString" 			./input.txt
#sed -i 	"/nQy/c\   nQy         = $nQx 	$infoString" 			./input.txt
#sed -i 	"/nSCx/c\   nSCx        = $nSCx 	$infoString" 		./input.txt
#sed -i 	"/nSCy/c\   nSCy        = $nSCx 	$infoString" 		./input.txt
#sed -i 	"/nKx/c\   nKx         = $nKx 	$infoString" 			./input.txt
#sed -i 	"/nKy/c\   nKy         = $nKx 	$infoString" 			./input.txt
sed -i 	"/shells/c\   shells       = $nShell 	$infoString" 	./input.txt

#set Hamiltionian
sed -i 	"/doVdesc/c\	doVdesc		= f 	$infoString"		./input.txt
sed -i 	"/doZeeman/c\	doZeeman	= f 	$infoString"		./input.txt
sed -i 	"/doMagHam/c\	doMagHam	= f 	$infoString"		./input.txt
sed -i "/doRashba/c\    doRashba    = f     $infoString"        ./input.txt

#folder to save data
label='GcutTest'
mkdir -p $label
now="$(date +"%T")"
echo '['$now']: created directory '$label
printf "*\n*\n"







#body
for gCut in ${Gvalues[*]}; do
	#
	echo '['$(date +"%T")']: start Gcut='$gCut' a0^-1; '
	#
	#creates directories
	dir=$label'/gCut'$gCut
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
	sed -i  "/Gcut/c\   Gcut        = $gCut 	$infoString"	 	./input.txt
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
