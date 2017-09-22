program main
	!TWO dimensional potential well code
	use omp_lib
	use mathematics, 	only: 		dp, PI_dp, acc, eigSolver

	use sysPara
	use potWellModel, 	only: 		solveHam
	use wannier,	 	only: 		wannMethod	
	use berry,			only:		berryMethod

	use polarization,	only:		calcPolWannCent, calcPolViaA, calcIonicPol
	use semiclassics,	only:		calcFirstOrdP
	use peierls,		only:		peierlsSub
	use output,		 	only:		writeMeshInfo, writeMeshBin, writeWaveFunc, writeWannFiles,writePolFile,& 
									printTiming	!printMat, printInp, printWannInfo,writeSysInfo  



	implicit none

	
    real(dp), 		allocatable,	dimension(:,:)		:: 	wCent, wSprd
    complex(dp),	allocatable,	dimension(:,:,:)	:: 	wnF, unk	!, ukn basCoeff,
    complex(dp),	allocatable,	dimension(:,:,:,:)	::	vW,Velo	, veloBwf			!, Ah, Fh, Vh
    real(dp),		allocatable,	dimension(:,:)		:: 	En					!, EnH
    real(dp),		allocatable,	dimension(:,:,:)	::	Fw
    real(dp),		allocatable,	dimension(:,:,:,:)	::	Aw
    real(dp) 											:: 	pWann(2), pIon(2), pTot(2), pBerry(2), pInt(2), pNiu(3), pPei(3)
    real												:: 	mastT0,mastT1,mastT,aT0,aT1,aT,kT0,kT1,kT,wT0,wT1,wT, oT0, oT1, oT,&
    														wI0, wI1, wI, scT0, scT1, scT, peiT, peiT0, peiT1
    integer 											::	xi,ki
    call cpu_time(mastT0)






    !INPUT & ALLOCATION SECTION
    call cpu_time(aT0)
	call readInp()
	!electronic structure arrays
	allocate(			unk(		nR 		,	nWfs	, nQ		)				)
	allocate(			En(						nWfs	, nQ		)				)
	allocate(			VeloBwf(2,	nR		, 	nG		, nQ		)				) 
	!allocate(			wnF( 		nR		, 	nSC		, nWfs		)				)
	!allocate(			wCent(		2		, 	nWfs				)				)
	!allocate(			wSprd(		2		, 	nWfs				)				)
	

	!coarse mesh quanties
	
	!allocate(			Aw(			3		,	nQ		, nWfs,nWfs	)				)
	!allocate(			vW(			3		,	nQ		, nWfs,nWfs	)				)
	!allocate(			FW(			3		,	nQ		, nWfs		)				)

	!wannier interpolation arrays
	!allocate(			Ah(		3		,	nK		, nWfs, nWfs	)			)
	!allocate(			Fh(		3		,	nK		, nWfs, nWfs	)			)
	!allocate(			Vh(		3		,	nK		, nWfs, nWfs	)			)
	!allocate(			EnH(				nK		,	 nWfs		)			)
	
	!
	write(*,*)"[main]:**************************Infos about this run*************************"
	write(*,*)"[main]: nK=",nK
	write(*,*)"[main]: nQ=",nQ
	
	!
	call cpu_time(aT1)
	aT = aT1 - aT0



	
	!SOLVE ELECTRONIC STRUCTURE & GENERATE THE WANNIER FUNCTIONS ON THE FLY
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************ELECTRONIC STRUCTURE PART*************************"
	call cpu_time(kT0)
	!
	!
	call solveHam(unk, En, VeloBwf)
	
	!check boundary condition on unk files
	!if( isKperiodic(unk)	 /= 0 ) then			!test is questionable since k mesh goes from kmin till -kmin-dk 
	!	write(*,*)	"[main]: problem with unk gauge, wave functions are NOT periodic in k space"
	!end if
	!
	call cpu_time(kT1)
	write(*,*)"[main]: done solving Schroedinger eq."
	kT = kT1-kT0
	



	
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	call cpu_time(wT0)
	if( doWanni ) then
		write(*,*)	"[main]:**************************WANNIER FUNCTION METHOD*************************"
		!
		call wannMethod(unk, pWann)
		!
		write(*,*)	"[main]: done with center polarization calc"
	else
		write(*,*)	"[main]: wannier method disabled"
	end if	
	
	call cpu_time(wT1)
	
	wT 	= wT1 - wT0




	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	call cpu_time(wI0)
	if ( doBerry ) then
		write(*,*)"[main]:**************************WAVEFUNCTION METHOD*************************"
		call berryMethod(unk, En, veloBwf, pBerry, pNiu)
		write(*,*)"[main]: done with wavefunction method "
	else
		write(*,*)"[main]: berry method disabled"
	end if
	call cpu_time(wI1)
	wI	= wI1 - wI0
	











	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"[main]:**************************WANNIER INTERPOLATION*************************"
	!call cpu_time(wI0)
	!call calcConnCurv(unk, wnF, Ah, Fh, Vh, EnH)
	!write(*,*)"[main]: interpolation of connection & curvature done, test if they are real now and calc pol"
	!call testIfReal(Ah, Fh)
	!call calcPolViaA(dreal(Ah), pInt)		
	!call cpu_time(wI1)
	!write(*,*)"[main]: done with wannier interpolation"
	!wI	= wI1 - wI0
!
!











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
	
	
	call calcIonicPol(pIon)
	pTot	= pIon + pWann
	call writePolFile(pWann, pIon, pTot, pBerry, pInt, pNiu, pPei )
	write(*,*)"[main]: ...wrote polarization txt file"
	!
	
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

