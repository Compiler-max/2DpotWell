program main
	!TWO dimensional potential well code
	use mathematics, 	only: 		dp, PI_dp, eigSolver
	!
	use sysPara    , 	only: 		readInp, &
									aX, aY,vol, nAt, relXpos, relYpos, atRx, atRy, atPot,&
									nG, Gcut, nK, nWfs, nSC, nR, dx, dy, dkx, dky, &
									Gvec, atPos, atR, kpts, rpts
	!
	use potWellModel, 	only: 		solveHam 
	!
	use wannier	   ,	only: 		isNormal			!, calcWcent, calcWsprd,calc0ElPol, &
	!						genWUnk ,bandInterpol,gaugeUnk, calcConnViaK, gaugeConnToHam, calcPolViaA
	use semiClassics,	only:		dummy
	!
	use output	   , only:		writeMeshInfo, writeMeshBin, writeWannFiles, printTiming	!printMat, printInp, printWannInfo,& 
	!						writeSysInfo  
	implicit none

	
    real(dp)   	, allocatable, dimension(:,:)		:: wCent, wSprd
    complex(dp)	, allocatable, dimension(:,:,:)	:: wnF, unkW, Uh , Aconn!, ukn basCoeff,
    real(dp) 									:: pEl, pIon, pTot, pElViaA 
    real										:: mastT0,mastT1,mastT,aT0,aT1,aT,kT0,kT1,kT,wT0,wT1,wT, oT0, oT1, oT
    integer 									:: xi,ki
    call cpu_time(mastT0)







    !INPUT & ALLOCATION SECTION
    call cpu_time(aT0)
	call readInp()
	!
	allocate(		wnF( 	nR, nSC, nWfs	)		)
	allocate(		unkW(	nR, nK, nWfs	)		) 
	allocate(		Uh(		nWfs, nWfs,	nK	)		)
	allocate(		Aconn(	nWfs, nWfs, nK	)		)
	allocate(		wCent(	2, nWfs			)		)
	allocate(		wSprd(	2, nWfs			)		)

	!
	call cpu_time(aT1)
	aT = aT1 - aT0

	
	!SOLVE ELECTRONIC STRUCTURE & GENERATE THE WANNIER FUNCTIONS ON THE FLY
	write(*,*)"[main]: start solving Schroedinger eq."
	call cpu_time(kT0)
	!
	call solveHam(Uh, wnF, unkW)
	!
	call cpu_time(kT1)
	write(*,*)"[main]: done solving Schroedinger eq."
	kT = kT1-kT0


	if( .not. isNormal(wnF) ) then
		write(*,*)"[main]: the used Wannier functions are not properly normalized"
	else
		write(*,*)"[main]: Wannier functions in home unit cell are properly normalized"
	end if

	




	
	!!WANNIER CENTERS & POLARIZATION
	call cpu_time(wT0)
	!call calcWcent(wnF,	wCent)
	!call calcWsprd(wnF, wCent, wSprd)
	!call calc0ElPol(wCent, pEl, pIon, pTot)
	!write(*,*)"[main]: done with center polarization calc"

	
	!call calcConnViaK(unkW, Aconn)   
	!write(*,*)"[main]: connection calculated"
	!!call gaugeConnToHam(unkW, Uh, Aconn)

	!call calcPolViaA(Aconn, pElViaA)
	!write(*,*)"[main]: pol via connection calculated"
	call cpu_time(wT1)
	wT = wT1 - wT0

	!OUTPUTING RESULTS SECTION
	call cpu_time(oT0)
	!call writeSysInfo() 
	call writeMeshInfo() 
	call writeMeshBin()
	call writeWannFiles(wnF)
	!call writeWannFiles(gnr, wnF, wCent, wSprd, Aconn, unkW)
	call cpu_time(oT1)
	oT = oT1 - oT0
	
	!TIMING INFO SECTION
	!call printInp()
	!call printWannInfo(wCent, wSprd, pEl, pIon, pTot, pElViaA)
	call cpu_time(mastT1)
	mastT= mastT1-mastT0
	write(*,*) "[main]: timing info:"
	call printTiming(aT, kT, wT, oT, mastT)

	stop
end 

