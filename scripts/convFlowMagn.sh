#!/bin/bash

#vars
infoString='# edited by convFlow'
root=$PWD


#numerics
label='convTest_Magn'
gVal=( 5.0 )		#( 20.0 25.0 30.0 )
qVal=( 4 8 16 )
nSolve=24
cells=3





#initial prepartions


printf "*\n***** script for performing several Gcut iterations************\n*\n"

mkdir -p $label
now="$(date +"%T")"
echo '['$now']: created directory '$label
printf "*\n*\n"



#loop data
for g in ${gVal[*]}; do
	for q in ${qVal[*]}; do
		#
		#get q values
		nQx=$q
		nQy=$(($nQx * $cells))
		nSCx=$(($nQx + 1))
		nSCy=$(($nQy + 1))
		#
		#determine nShell depending of q grid spacing
		if (( $q <= 4 )); then
			nShell=1
		elif (( $q <= 8 )); then
			nShell=1
		else
			nShell=1
		fi
		#
		#get directory name	
		dir=$label'/g'$g'q'$q
		echo '['$(date +"%T")']: start Gcut='$g' nQx='$nQx' nQy='$nQy' nShell='$nShell
		#	
		#creates directories
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
		cp magnStream.sh $dir
		wait
		#
		#prepare input file
		cd $dir
		#set numerics (for i1 & i2)
		sed -i  "/Gcut/c\    Gcut        = $g $infoString"	 		./input.txt
		sed -i  "/nSolve/c\    nSolve      = $nSolve $infoString"	./input.txt
		#
		sed -i 	"/nQx/c\    nQx         = $nQx $infoString" 			./input.txt
		sed -i 	"/nQy/c\    nQy         = $nQy $infoString" 			./input.txt
		sed -i 	"/nSCx/c\    nSCx        = $nSCx $infoString" 		./input.txt
		sed -i 	"/nSCy/c\    nSCy        = $nSCy $infoString" 		./input.txt
		sed -i 	"/nKx/c\    nKx         = $nSCx $infoString" 			./input.txt
		sed -i 	"/nKy/c\    nKy         = $nSCy $infoString" 			./input.txt
		#
		sed -i 	"/shells/c\    shells       = $nShell $infoString" 	./input.txt
		wait
		#	
		#execute calculation
		./magnStream.sh
		wait
		echo '['$(date +"%T")']: finished magnStream'
			
	
		#
		#jump back to root
		cd $root
		echo '['$(date +"%T")']: done with '$label'= '$g
		printf '*\n*\n'
	done
done



echo '['$(date +"%T")']: finished all calculations, by'	
echo '*******************************************************'
