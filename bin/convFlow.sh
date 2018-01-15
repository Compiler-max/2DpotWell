#!/bin/bash

#vars
infoString='# edited by convFlow'
root=$PWD


#numerics
label='convTest'
gVal=( 20.0 25.0 30.0 )
qVal=( 8 16 32 64 )
nSolve=60






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
		nQ=$q
		nSC=$((nQ+1))
		nK=$nSC
		#
		#determine nShell depending of q grid spacing
		if (( $q <= 4 )); then
			nShell=3
		elif (( $q <= 8 )); then
			nShell=2
		else
			nShell=1
		fi
		#
		#get directory name	
		dir=$label'/g'$g'q'$q
		echo '['$(date +"%T")']: start Gcut='$g' nQ='$q' nShell='$nShell
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
		cp subStream.sh $dir
		wait
		#
		#prepare input file
		cd $dir
		#set numerics (for i1 & i2)
		sed -i  "/Gcut/c\    Gcut        = $g $infoString"	 		./input.txt
		sed -i  "/nSolve/c\    nSolve      = $nSolve $infoString"	./input.txt
		#
		sed -i 	"/nQx/c\    nQx         = $nQ $infoString" 			./input.txt
		sed -i 	"/nQy/c\    nQy         = $nQ $infoString" 			./input.txt
		sed -i 	"/nSCx/c\    nSCx        = $nSC $infoString" 		./input.txt
		sed -i 	"/nSCy/c\    nSCy        = $nSC $infoString" 		./input.txt
		sed -i 	"/nKx/c\    nKx         = $nK $infoString" 			./input.txt
		sed -i 	"/nKy/c\    nKy         = $nK $infoString" 			./input.txt
		#
		sed -i 	"/shell/c\    shell       = $nShell $infoString" 	./input.txt
		wait
		#	
		#execute calculation
		./subStream.sh
		wait
		echo '['$(date +"%T")']: finished subStream'
			
	
		#
		#jump back to root
		cd $root
		echo '['$(date +"%T")']: done with '$label'= '$g
		printf '*\n*\n'
	done
done



echo '['$(date +"%T")']: finished all calculations, by'	
echo '*******************************************************'
