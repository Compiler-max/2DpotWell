program main
	!TWO dimensional potential well code
	use mathematics, 	only: 		dp, PI_dp, acc, eigSolver

	use sysPara
	use potWellModel, 	only: 		solveHam
	use berry,			only:		isKperiodic!, calcCurv, calcVeloMat 
	use wannier,	 	only: 		isNormal, calcWcent, calcWsprd, genUnkW 	
	use gaugeTrafo,		only:		calcConnCurv,  testIfReal
	use polarization,	only:		calcPolWannCent, calcPolViaA
	use semiclassics,	only:		calcFirstOrdP
	use peierls,		only:		peierlsSub
	use output,		 	only:		writeMeshInfo, writeMeshBin, writeWaveFunc, writeWannFiles,writePolFile,& 
									printTiming	!printMat, printInp, printWannInfo,writeSysInfo  



	implicit none

	
    real(dp), 		allocatable,	dimension(:,:)		:: 	wCent, wSprd
    real(dp),		allocatable,	dimension(:,:,:)	::	Fcurv
    complex(dp),	allocatable,	dimension(:,:,:)	:: 	wnF, unk, Uh, Aint, veloBwf, bWf !, ukn basCoeff,
    complex(dp),	allocatable,	dimension(:,:,:,:)	::	Velo, Ah, Fh, Vh
    real(dp),		allocatable,	dimension(:,:)		:: 	En, EnH
    real(dp),		allocatable,	dimension(:,:,:,:)	::	Aconn
    real(dp) 											:: 	pEl(2), pIon(2), pTot(2), pElViaA(2), pInt(2), pNiu(3), pPei(3)
    real												:: 	mastT0,mastT1,mastT,aT0,aT1,aT,kT0,kT1,kT,wT0,wT1,wT, oT0, oT1, oT,&
    														wI0, wI1, wI, scT0, scT1, scT, peiT, peiT0, peiT1
    integer 											::	xi,ki
    call cpu_time(mastT0)






    !INPUT & ALLOCATION SECTION
    call cpu_time(aT0)
	call readInp()
	!electronic structure arrays
	allocate(			unk(	nR		, 	nQ		, nWfs		)				) 
	!wannier functions
	allocate(			wnF( 		nR	, 	nSC		, nWfs		)				)
	allocate(			wCent(		2	, 	nWfs				)				)
	allocate(			wSprd(		2	, 	nWfs				)				)
	!wannier interpolation arrays
	allocate(			Ah(		3		,	nK		, nWfs, nWfs	)			)
	allocate(			Fh(		3		,	nK		, nWfs, nWfs	)			)
	allocate(			Vh(		3		,	nK		, nWfs, nWfs	)			)
	allocate(			EnH(				nK		,	 nWfs		)			)
	
	!
	!allocate(			En(					nQ		,	nWfs		)			)
	!allocate(			Uh(		nWfs	, 	nWfs	,	nQ		)			)
	!allocate(			Aint(		2		,	nK		, nWfs	)				)
	!allocate(			VeloBwf(	nR		,	nQ		, 2*nWfs)				)
	!allocate(			Velo(		3		,	nQ		,nWfs, nwFs)			)
	
	write(*,*)"[main]: nK=",nK
	write(*,*)"[main]: nQ=",nQ
	!
	call cpu_time(aT1)
	aT = aT1 - aT0

	
	!SOLVE ELECTRONIC STRUCTURE & GENERATE THE WANNIER FUNCTIONS ON THE FLY
	write(*,*)"[main]:**************************ELECTRONIC STRUCTURE PART*************************"
	call cpu_time(kT0)
	!
	!
	call solveHam(wnF, unk)
	!check boundary condition on unk files
	if( isKperiodic(unk)	 /= 0 ) then
		write(*,*)	"[main]: problem with unk gauge, wave functions are NOT periodic in k space"
	end if
	
	call cpu_time(kT1)
	write(*,*)"[main]: done solving Schroedinger eq."
	kT = kT1-kT0

	




	
	!!WANNIER CENTERS & POLARIZATION
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************WANNIER CENTERS & POLARIZATION*************************"
	call cpu_time(wT0)
	!
	!
	!if(  isNormal(wnF) /= 0 ) then
	!	write(*,*)"[main]: the used Wannier functions are not properly normalized"
	!else
	!	write(*,*)"[main]: Wannier functions in home unit cell are properly normalized"
	!end if
	!
	write(*,'(a)')"[main]: note that wannier function normalization is not tested currently!"
	!
	call calcWcent(wnF,	wCent)
	call calcWsprd(wnF, wCent, wSprd)
	call calcPolWannCent(wCent, pEl, pIon, pTot)
	!
	!
	call cpu_time(wT1)
	write(*,*)"[main]: done with center polarization calc"
	wT 	= wT1 - wT0






	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************WANNIER INTERPOLATION*************************"
	call cpu_time(wI0)
	call calcConnCurv(unk, wnF, Ah, Fh, Vh, EnH)
	write(*,*)"[main]: interpolation of connection & curvature done, test if they are real now and calc pol"
	call testIfReal(Ah, Fh)
	call calcPolViaA(dreal(Ah), pInt)		
	call cpu_time(wI1)
	write(*,*)"[main]: done with wannier interpolation"
	wI	= wI1 - wI0







	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"[main]:**************************SEMICLASSICS*************************"
	!call cpu_time(scT0)
	!call calcFirstOrdP(dreal(Fh), dreal(Ah), Vh, EnH, pNiu)
	!call cpu_time(scT1)
	!write(*,*)"[main]: done with first order polarization calculation"
	!scT	= scT1 - scT0








	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"[main]:**************************PEIERLS SUB*************************"
	!call cpu_time(peiT0)
	!!
	!call peierlsSub(wnF, 	unk, Ah, Fh, Vh, EnH, pPei )
	!!
	!call cpu_time(peiT1)
	!write(*,*)"[main]: done with peierls substitution"
	!peiT = peiT1 - peiT0







	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************WRITE OUTPUT*************************"
	call cpu_time(oT0)
	!call writeSysInfo() 
	call writeMeshInfo() 
	call writeMeshBin()
	write(*,*)"[main]: ...wrote mesh info"
	call writeWaveFunc(unk, Aconn, Fcurv)
	write(*,*)"[main]: ...wrote unks and connection"
	call writeWannFiles(wnF, wCent, wSprd)			!call writeWannFiles(gnr, wnF, wCent, wSprd, Aconn, unkW)
	write(*,*)"[main]: ...wrote wannier functions"
	call writePolFile(pEl, pIon, pTot, pElViaA, pInt, pNiu, pPei )
	write(*,*)"[main]: ...wrote polarization txt file"
	!
	write(*,'(a,E10.3)')"[main]: overlap warnings etc. where given for diviations from expected above the threshold: ",acc
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
	call printTiming(aT, kT, wI, scT, peiT, wT, oT, mastT)

	stop
end program

