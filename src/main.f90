program main
	!TWO dimensional potential well code
	use mathematics, only: dp, PI_dp, eigSolver
	!
	use sysPara    , only: 		readInp, &
									aX, aY,vol, nAt, relXpos, relYpos, atRx, atRy, atPot,&
									nG, Gcut, nK, nWfs, nSC, nR, dx, dy, dkx, dky, &
									Gvec, atPos, atR, kpts, rpts
	!
	use potWellModel, only: 		solveHam 
	!
	use wannier	   , only: isNormal,calcWcent, calcWsprd,calc0ElPol , &
							genUnkW, calcConn, calcPolViaA  !,bandInterpol,gaugeUnk, calcConnViaK, gaugeConnToHam
	!
	use output	   , only:		writeMeshInfo, writeMeshBin, writeWaveFunc, writeWannFiles,writePolFile,& 
								printTiming	!printMat, printInp, printWannInfo,writeSysInfo  


						
	implicit none

	
    real(dp), 		allocatable,	dimension(:,:)		:: wCent, wSprd
    complex(dp),	allocatable,	dimension(:,:,:)	:: wnF, unk, Uh, Aconn !, ukn basCoeff,
    real(dp) 											:: pEl(2), pIon(2), pTot(2), pElViaA(2) 
    real												:: mastT0,mastT1,mastT,aT0,aT1,aT,kT0,kT1,kT,wT0,wT1,wT, oT0, oT1, oT
    integer 											:: xi,ki
    call cpu_time(mastT0)







    !INPUT & ALLOCATION SECTION
    call cpu_time(aT0)
	call readInp()
	!
	allocate(			wnF( 		nR		, 	nSC		, nWfs	)				)
	allocate(			unk(		nR		, 	nK		, nWfs	)				) 
	allocate(			Uh(			nWfs	, 	nWfs	,	nK	)				)
	allocate(			Aconn(		2		,	nK		, nWfs	)				)
	allocate(			wCent(		2		, 	nWfs			)				)
	allocate(			wSprd(		2		, 	nWfs			)				)

	!
	call cpu_time(aT1)
	aT = aT1 - aT0

	
	!SOLVE ELECTRONIC STRUCTURE & GENERATE THE WANNIER FUNCTIONS ON THE FLY
	write(*,*)"[main]:**************************ELECTRONIC STRUCTURE PART*************************"
	call cpu_time(kT0)
	!
	call solveHam(Uh, wnF, unk)
	!
	call cpu_time(kT1)
	write(*,*)"[main]: done solving Schroedinger eq."
	kT = kT1-kT0

	




	
	!!WANNIER CENTERS & POLARIZATION
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************WANNIER CENTERS & POLARIZATION*************************"
	if( .not. isNormal(wnF) ) then
		write(*,*)"[main]: the used Wannier functions are not properly normalized"
	else
		write(*,*)"[main]: Wannier functions in home unit cell are properly normalized"
	end if

	call cpu_time(wT0)
	call calcWcent(wnF,	wCent)
	call calcWsprd(wnF, wCent, wSprd)
	call calc0ElPol(wCent, pEl, pIon, pTot)
	write(*,*)"[main]: done with center polarization calc"

	
	call calcConn(unk, Aconn)   
	call calcPolViaA(Aconn, pElViaA)
	call cpu_time(wT1)
	wT = wT1 - wT0




	!OUTPUTING RESULTS SECTION
	write(*,*)"[main]: everything done, start writing results"

	call cpu_time(oT0)
	
	!call writeSysInfo() 
	call writeMeshInfo() 
	call writeMeshBin()
	call writeWaveFunc(unk, Aconn)
	call writeWannFiles(wnF, wCent, wSprd)			!call writeWannFiles(gnr, wnF, wCent, wSprd, Aconn, unkW)
	call writePolFile(pEl, pIon, pTot, pElViaA )
	
	call cpu_time(oT1)
	oT = oT1 - oT0
	




	!TIMING INFO SECTION
	!call printInp()
	!call printWannInfo(wCent, wSprd, pEl, pIon, pTot, pElViaA)
	call cpu_time(mastT1)
	mastT= mastT1-mastT0
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*) '**************TIMING INFORMATION************************'
	call printTiming(aT, kT, wT, oT, mastT)

	stop
end 

