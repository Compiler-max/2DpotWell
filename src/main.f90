program main
	!TWO dimensional potential well code
	use omp_lib
	use mathematics, 	only: 		dp, PI_dp, acc, eigSolver

	use sysPara
	use potWellModel, 	only: 		solveHam
	use projection,		only:		projectUnk
	use wannier,	 	only: 		wannMethod	
	use berry,			only:		berryMethod

	use polarization,	only:		calcIonicPol
	use semiclassics,	only:		calcFirstOrdP
	use peierls,		only:		peierlsMethod
	use output,		 	only:		writeMeshInfo, writeMeshBin, writeUNKs, writeInterpBands, writeWannFiles,writePolFile,& 
									printTiming	!printMat, printInp, printWannInfo,writeSysInfo  



	implicit none

	
    real(dp), 		allocatable,	dimension(:,:)		:: 	wCent, wSprd
    complex(dp),	allocatable,	dimension(:,:,:)	:: 	wnF, unk, unkP, tHopp	!, ukn basCoeff,
    complex(dp),	allocatable,	dimension(:,:,:,:)	::	veloBwf, veloP			!, Ah, Fh, Vh
    real(dp),		allocatable,	dimension(:,:)		:: 	En, EnP					!, EnH
    real(dp),		allocatable,	dimension(:,:,:)	::	Fw
    real(dp),		allocatable,	dimension(:,:,:,:)	::	Aw
    real(dp) 											:: 	pWann(2), pIon(2), pTot(2), pBerry(2), pInt(2), pNiu(3), pPei(3)
    real												:: 	mastT0,mastT1,mastT,aT0,aT1,aT,kT0,kT1,kT,wT0,wT1,wT, oT0, oT1, oT,&
    														bT0, bT1, bT, scT0, scT1, scT, peiT, peiT0, peiT1, pT0, pT1, pT
    integer 											::	xi,ki
    call cpu_time(mastT0)




    !INPUT & ALLOCATION SECTION
    call cpu_time(aT0)
	call readInp()
	!electronic structure arrays
	allocate(			unk(		nR 		,	nG		, 	nQ		)				)
	allocate(			unkP(		nR 		,	nWfs	, 	nQ		)				)
	allocate(			En(						nBands	, 	nQ		)				)
	allocate(			EnP(					nWfs	,	nQ		)				)
	allocate(			tHopp(		nWfs	, 	nWfs	,	nSc		)				)
	allocate(			veloP(	3,	nWfs	,	nWfs	,	nQ		)				)

	allocate(			veloBwf(2,	nR		, 	nBands	, 	nQ		)				) 

	!wannier interpolation arrays
	!allocate(			Ah(		3		,	nK		, nWfs, nWfs	)			)
	!allocate(			Fh(		3		,	nK		, nWfs, nWfs	)			)
	!allocate(			Vh(		3		,	nK		, nWfs, nWfs	)			)
	!allocate(			EnH(				nK		,	 nWfs		)			)
	
	!
	write(*,*)"[main]:**************************Infos about this run*************************"
	write(*,*)"[main]: nK=",nK
	write(*,*)"[main]: nQ=",nQ

    write(*,*)"[main]: nBands=", nBands
	write(*,*)"[main]: nWfs  =", nWfs	
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
	!
	call cpu_time(kT1)
	write(*,*)"[main]: done solving Schroedinger eq."
	kT = kT1-kT0
	



	!PROJECT THE STATES
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************PROJECT STATES *************************"
	call cpu_time(pT0)
	!
	!
	call projectUnk(En, unk, EnP, unkP,tHopp)
	!
	call cpu_time(pT1)
	write(*,*)"[main]: done with projections."
	pT = pT1-pT0


	
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	call cpu_time(wT0)
	if( doWanni ) then
		write(*,*)	"[main]:**************************WANNIER FUNCTION METHOD*************************"
		!
		call wannMethod(unkP, pWann)
		!
		write(*,*)	"[main]: done with center polarization calc"
	else
		write(*,*)	"[main]: wannier method disabled"
	end if	
	
	call cpu_time(wT1)
	
	wT 	= wT1 - wT0



	!TODO REWRITE WITH EnP, veloP
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	call cpu_time(bT0)
	if ( doBerry ) then
		write(*,*)"[main]:**************************WAVEFUNCTION METHOD*************************"
		call berryMethod(unkP, EnP, pBerry, pNiu)
		write(*,*)"[main]: done with wavefunction method "
	else
		write(*,*)"[main]: berry method disabled"
	end if
	call cpu_time(bT1)
	bT	= bT1 - bT0
	















	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************PEIERLS SUB*************************"
	call cpu_time(peiT0)
	!
	if(doPei)  then
		call peierlsMethod(tHopp, pPei)
	end if
	!
	call cpu_time(peiT1)
	write(*,*)"[main]: done with peierls substitution"
	peiT = peiT1 - peiT0







	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************WRITE OUTPUT*************************"
	call cpu_time(oT0)
	!call writeSysInfo() 
	call writeMeshInfo() 
	call writeMeshBin()
	call writeUNKs(unkP)
	call writeInterpBands(EnP)
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
	call printTiming(aT, kT, pT, wT, bT, oT, mastT)

	deallocate(			unk			)
	deallocate(			unkP		)
	deallocate(			En			)
	deallocate(			veloBwf		) 
	deallocate(			EnP			)
	deallocate(			tHopp		)

	stop
end program

