#!/bin/bash


#VARS
wannDIR=$work/thirdparty/wannier90
infoString='# edited by subStream'


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
	mkdir -p rawData
	rm -r -f rawData/*
	#i1
	sed -i "/doSolveHam/c\    doSolveHam  = t $infoString" ./input.txt
	sed -i "/doPw90/c\    doPw90      = f $infoString" ./input.txt
	sed -i "/doBerry/c\    doBerry     = f $infoString" ./input.txt
	sed -i "/useRot/c\    useRot     = f $infoString" ./input.txt
	#i2
	cp input.txt input2.txt
	sed -i "/doSolveHam/c\    doSolveHam  = f $infoString" ./input2.txt
	sed -i "/doPw90/c\    doPw90      = t $infoString" ./input2.txt
	sed -i "/doBerry/c\    doBerry     = t $infoString" ./input2.txt
	#i3
	cp input2.txt input3.txt
	sed -i "/useRot/c\    useRot      = t $infoString" ./input3.txt

}

function runCalc {
	#calls the exectutables
	#echo 'calculations disabled'
	#electronic structure
	./main.exe > outABin.txt
	wait
	echo '['$(date +"%T")']: finished electronic structure'
	
	#wannier90
	$wannDIR/wannier90.x wf1
	wait
	$wannDIR/postw90.x wf1
	wait
	echo '['$(date +"%T")']: finished wannier90 calculations'
	
	#post-w90
	mv input.txt inputABinit.txt
	mv input2.txt input.txt
	wait
	./main.exe > out.txt
	wait
	mv input.txt inputNoROT.txt
	mv input3.txt input.txt
	mv polOutput.txt polOutputNoROT.txt
	mkdir rawDataNoROT
	cp ./rawData/* rawDataNoROT
	wait
	./main.exe > outROT.txt
	wait
	mv rawData rawDataROT
	mv input.txt inputROT.txt
	mv polOutput.txt polOutputROT.txt
	echo '['$(date +"%T")']: finished post wannier90 run'	
}





#BODY
#
#search foor w90
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














