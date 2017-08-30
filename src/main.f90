program main
	!TWO dimensional potential well code
	use mathematics, 	only: 		dp, PI_dp, eigSolver
	!
	use sysPara,	 	only: 		readInp, &
									aX, aY,vol, nAt, relXpos, relYpos, atRx, atRy, atPot,&
									nG, Gcut, nK, nKx, nKy, nKw, nKxW, nKyW, nWfs, nSC, nR, dx, dy, dkx, dky, &
									Gvec, atPos, atR, kpts, rpts
	!
	use potWellModel, 	only: 		solveHam, calcVeloMat, calcConn 
	!
	use wannier,	 	only: 		isNormal,calcWcent, calcWsprd,calc0ElPol , &
									genUnkW, calcPolViaA, interpConnCurv  !,bandInterpol,gaugeUnk, calcConnViaK, gaugeConnToHam
	!
	use semiclassics,	only:		calcFirstOrdP
	!
	use output,		 	only:		writeMeshInfo, writeMeshBin, writeWaveFunc, writeWannFiles,writePolFile,& 
									printTiming	!printMat, printInp, printWannInfo,writeSysInfo  



	implicit none

	
    real(dp), 		allocatable,	dimension(:,:)		:: 	wCent, wSprd
    complex(dp),	allocatable,	dimension(:,:,:)	:: 	wnF, unk, unkW, Uh, Aconn, Fcurv, Aint, veloBwf, bWf !, ukn basCoeff,
    complex(dp),	allocatable,	dimension(:,:,:,:)	::	Velo
    real(dp),		allocatable,	dimension(:,:)		:: 	En
    real(dp) 											:: 	pEl(2), pIon(2), pTot(2), pElViaA(2), pInt(2), p1(3)
    real												:: 	mastT0,mastT1,mastT,aT0,aT1,aT,kT0,kT1,kT,wT0,wT1,wT, oT0, oT1, oT,&
    														wI0, wI1, wI, scT0, scT1, scT
    integer 											::	xi,ki
    call cpu_time(mastT0)




   write(*,*)nWfs


    !INPUT & ALLOCATION SECTION
    call cpu_time(aT0)
	call readInp()
	!
	allocate(			En(			nWfs		,	nK)							)
	allocate(			wnF( 		nR		, 	nSC		, nWfs	)				)
	allocate(			unk(		nR		, 	nK		, nWfs	)				) 
	allocate(			Uh(			nWfs	, 	nWfs	,	nK	)				)
	allocate(			Aconn(		2		,	nK		, nWfs	)				)
	allocate(			Fcurv(		3		,	nK		, nWfs	)				)
	allocate(			wCent(		2		, 	nWfs			)				)
	allocate(			wSprd(		2		, 	nWfs			)				)
	!
	!allocate(			Aint(		2		,	nKw		, nWfs	)				)
	allocate(			VeloBwf(	nR		,	nK		, 2*nWfs)				)
	allocate(			Velo(		3		,	nK		,nWfs, nwFs)			)
	
	
	
	!
	call cpu_time(aT1)
	aT = aT1 - aT0

	
	!SOLVE ELECTRONIC STRUCTURE & GENERATE THE WANNIER FUNCTIONS ON THE FLY
	write(*,*)"[main]:**************************ELECTRONIC STRUCTURE PART*************************"
	call cpu_time(kT0)
	!
	call solveHam(Uh, wnF, unk, En, veloBwf)
	!Get Velocity operator matrix
	call calcVeloMat(unk, veloBwf, Velo)

	!
	call calcConn(unk, nKx, nKy, Aconn)   
	call calcPolViaA(Aconn, pElViaA)
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

	
	
	call cpu_time(wT1)
	wT 	= wT1 - wT0



	write(*,*)"[main]:**************************WANNIER INTERPOLATION*************************"
	call cpu_time(wI0)
	!call interpConnCurv(wnF, Aint, Fcurv)
	!call calcPolViaA(Aint, pInt)
	call cpu_time(wI1)
	wI	= wI1 - wI0


	write(*,*)"[main]:**************************SEMICLASSICS*************************"
	call cpu_time(scT0)
	call calcFirstOrdP(Fcurv, Aconn, Velo, En, p1)
	call cpu_time(scT1)
	scT	= scT1 - scT0

	!OUTPUTING RESULTS SECTION
	write(*,*)"[main]:**************************WRITE OUTPUT*************************"
	call cpu_time(oT0)
	!call writeSysInfo() 
	call writeMeshInfo() 
	call writeMeshBin()
	write(*,*)"[main]: ...wrote mesh info"
	call writeWaveFunc(unk, Aconn)
	write(*,*)"[main]: ...wrote unks and connection"
	call writeWannFiles(wnF, wCent, wSprd)			!call writeWannFiles(gnr, wnF, wCent, wSprd, Aconn, unkW)
	write(*,*)"[main]: ...wrote wannier functions"
	call writePolFile(pEl, pIon, pTot, pElViaA, pInt )
	write(*,*)"[main]: ...wrote polarization txt file"
	
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
	call printTiming(aT, kT, wT, oT, wI, scT, mastT)

	stop
end 

