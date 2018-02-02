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

	mkdir -p oldInput
	rm -r -f oldInput/*
	#i1
	sed -i "/doSolveHam/c\    doSolveHam  = t $infoString" ./input.txt
	sed -i "/doPrepW90/c\    doPrepW90 = f $infoString" ./input.txt
	sed -i "/doPw90/c\    doPw90      = f $infoString" ./input.txt
	sed -i "/doBerry/c\    doBerry     = f $infoString" ./input.txt
	sed -i "/useRot/c\    useRot     = f $infoString" ./input.txt
	#i2
	cp input.txt input2.txt
	sed -i "/doSolveHam/c\    doSolveHam  = f $infoString" ./input2.txt
	sed -i "/doPrepW90/c\     doPrepW90 = t $infoString"	./input2.txt
	sed -i "/doPw90/c\    doPw90      = f $infoString" ./input2.txt
	sed -i "/doBerry/c\    doBerry     = f $infoString" ./input2.txt
	#i3
	cp input2.txt input3.txt
	sed -i "/doPrepW90/c\     doPrepW90 = f $infoString"	./input3.txt
	sed -i "/doBerry/c\    doBerry     = t $infoString" ./input3.txt
	
	cp input3.txt input4.txt
	sed -i "/useRot/c\    useRot      = t $infoString" ./input4.txt

}

function runCalc {
	#calls the exectutables
	#echo 'calculations disabled'
	

	#electronic structure
	mpirun -n 4 ./main.exe > outABin.txt
	wait
	echo '['$(date +"%T")']: finished electronic structure'
	
	#prepW90
	mv input.txt oldInput/inputABinit.txt
	mv input2.txt input.txt
	wait
	./main.exe > out.txt

	#run wannier90
	cd w90files
	$wannDIR/wannier90.x wf1
	wait
	$wannDIR/postw90.x wf1
	wait
	cd ..
	echo '['$(date +"%T")']: finished wannier90 calculations'
	
	
	#Berry (no rotation = random gauge)
	mv input.txt oldInput/inputW90prep.txt
	mv input3.txt input.txt
	wait
	./main.exe > outNoROT.txt
	wait
	echo '['$(date +"%T")']: finished berry random gauge'	

	#Berry2 (rotation = wannier gauge)
	mv input.txt oldInput/inputNoROT.txt
	mv input4.txt input.txt
	mv output/polOutput.txt output/polOutputNoROT.txt
	wait
	./main.exe > outROT.txt
	wait
	echo '['$(date +"%T")']: finished berry wannier gauge'	
	
	mv input.txt oldInput/inputROT.txt
	mv output/polOutput.txt output/polOutputROT.txt
	mv out*.txt output
	wait
	rm -r rawData 
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














