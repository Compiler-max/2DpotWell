#!/bin/bash


#VARS
wannDIR=$work/thirdparty/wannier90
infoString='# edited by subStream'
nprocs=$1


#FUNCTIONS
function searchWann {
	#search for wannier function subdirectory
	if [ -d "$wannDIR" ]; then
		if [[ -f "$wannDIR/wannier90.x" && -f "$wannDIR/postw90.x" ]]; then
			success=true
		else
			echo 'please compile wannier90 and postw90 in '$wannDIR
		fi
	else
		echo 'could not find wannier90 in '$wannDIR
	fi
}

function prepareInput {
	#setup workDir & tasks for the input files
	#
	#create working directory
	mkdir rawData
	rm -r -f rawData/*

	mkdir -p oldInput
	rm -r -f oldInput/*

	#save initial input file
	cp input.txt input.orig

	#e-Structure
	sed -i "/doSolveHam/c\    doSolveHam  = t $infoString" 	./input.txt
	sed -i "/doPw90/c\    doPw90      = f $infoString" 		./input.txt
	sed -i "/doBerry/c\    doBerry     = f $infoString" 	./input.txt
	
	#post w90 - berry method
	cp input.txt inputBerry.txt
	sed -i "/doSolveHam/c\    doSolveHam  = f $infoString" 	./inputBerry.txt
	sed -i "/doPw90/c\    doPw90      = f $infoString" 		./inputBerry.txt
	sed -i "/doBerry/c\    doBerry     = t $infoString" 	./inputBerry.txt
	sed -i "/doNiu/c\        doNiu = t $infoString" 		./inputBerry.txt
}

function runCalc {
	#calls the exectutables
	#echo 'calculations disabled'
	

	#electronic structure
	rm w90files/*
	wait
	mpirun -np $nprocs ./main.exe > outABin.txt
	wait
	mv input.txt oldInput/inputABiN.txt
	rm wf1.win #remove input for wannier setup
	mv wf1* w90files
	echo '['$(date +"%T")']: finished electronic structure'
	
	
	#run wannier90
	cd w90files
	$wannDIR/wannier90.x wf1
	wait
	$wannDIR/postw90.x wf1
	wait
	cd ..
	echo '['$(date +"%T")']: finished wannier90 calculations'
	
	
	#Berry
	mv inputBerry.txt input.txt
	wait
	./main.exe > outBerry.txt
	wait
	mv input.txt oldInput/inputBerry.txt
	echo '['$(date +"%T")']: finished berry method'	

	#finalize
	mv out*.txt output
	mv input.orig input.txt
	rm rawData/Mmn*
	rm rawData/Amn*
	rm rawData/ck*
	rm rawData/gVec*
	wait
	echo '['$(date +"%T")']: all done'	
}





#BODY
#
#search foor w90
echo '['$(date +"%T")']: start subStream'
success=false
searchWann
if [ ! success ]; then
	echo 'problems detecting wannier90, please set work=/path/2DpotWell and make sure wannier90.x and postw90.x are compiled'
else
	prepareInput #this is a function call
	wait
	runCalc	#this is a function call
	wait
fi














