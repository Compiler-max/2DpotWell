#!/bin/bash

#vars
infoString='# edited by GcutFlow'
root=$PWD


#numerics
label='convTest'
gVal=( 5.0 10.0) #15.0 20.0 25.0 30.0 )
qVal=( 4 8) # 16 32 64 128)
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
		case $nQ in 
			[1-4]*)
				nShell=3
				;;
			[5-8]*)
				nShell=2
				;;
			*)
				nShell=1
				;;
		esac
		echo 'set nShell as '$nShell
		#
		#get directory name	
		dir=$label'/g'$g'q'$q
		echo '['$(date +"%T")']: start Gcut='$g' nQ='$q
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
		echo '['$(date +"%T")']: start subStream'
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