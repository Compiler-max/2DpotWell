#!/bin/bash

#vars
infoString='# edited by statFlow'
root=$PWD



#numerics
gCut=25.0
nQ=16
nSC=$((nQ+1))
nK=$nSC
states=( 6 12  24 48 96 192 384 768 1536	3072 )


#derive shell to use from nQ
#determine nShell depending of q grid spacing
if (( $nQ <= 4 )); then
	nShell=3
elif (( $nQ <= 8 )); then
	nShell=2
else
	nShell=1
fi



printf "*\n***** script for performing several nState iterations************\n*\n"
label='nStates'
mkdir -p $label
now="$(date +"%T")"
echo '['$now']: created directory '$label
printf "*\n*\n"

echo '['$(date +"%T")']: numerics info Gcut='$gCut' nQ='$nQ' nShell='$nShell
printf "*\n*\n"


for s in ${states[*]}; do
	#get directory name	
	dir=$label'/nStat'$s
	echo '['$(date +"%T")']: start nStates='$s
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
	sed -i  "/Gcut/c\    Gcut        = $gCut $infoString"	 		./input.txt
	sed -i  "/nSolve/c\    nSolve      = $s $infoString"		./input.txt
	#
	sed -i 	"/nQx/c\    nQx         = $nQ $infoString" 			./input.txt
	sed -i 	"/nQy/c\    nQy         = $nQ $infoString" 			./input.txt
	sed -i 	"/nSCx/c\    nSCx        = $nSC $infoString" 		./input.txt
	sed -i 	"/nSCy/c\    nSCy        = $nSC $infoString" 		./input.txt
	sed -i 	"/nKx/c\    nKx         = $nK $infoString" 			./input.txt
	sed -i 	"/nKy/c\    nKy         = $nK $infoString" 			./input.txt
	#
	sed -i 	"/shells/c\    shells       = $nShell $infoString" 	./input.txt
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




echo '['$(date +"%T")']: finished all calculations, by'	
echo '*******************************************************'