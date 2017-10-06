program main
	!TWO dimensional potential well code
	use omp_lib
	use mathematics, 	only: 		dp

	use sysPara
	use potWellModel, 	only: 		solveHam
	use projection,		only:		projectUnk
	use wannier,	 	only: 		wannMethod	
	use berry,			only:		berryMethod

	use polarization,	only:		calcIonicPol
	use peierls,		only:		peierlsMethod
	use output,		 	only:		writeMeshInfo, writeMeshBin, writeUNKs, writeInterpBands, writeWannFiles,writePolFile,& 
									printTiming	!printMat, printInp, printWannInfo,writeSysInfo  



	implicit none

	
    complex(dp),	allocatable,	dimension(:,:,:)	:: 	unk, unkW, wnf, Uq!, ukn basCoeff,
    real(dp),		allocatable,	dimension(:,:)		:: 	En
    real(dp) 											:: 	pWann(2), pBerry(2), pNiu(3), pPei(3)
    real												:: 	mastT0, mastT1, mastT, T0, T1, &
    															aT,kT,wT, oT, bT, peiT, pT
    call cpu_time(mastT0)




    !INPUT & ALLOCATION SECTION
    call cpu_time(T0)
	call readInp()
	!electronic structure arrays
	allocate(			unk(		nR 		,	nG		, 	nQ		)				)
	allocate(			unkW(		nR 		,	nWfs	, 	nQ		)				)
	allocate(			Uq(			nBands	,	nWfs	, 	nQ		)				)
	allocate(			wnF( 		nR		, 	nSC		, nWfs		)				)
	allocate(			En(						nBands	, 	nQ		)				)
	
	


	!wannier interpolation arrays
	!allocate(			Ah(		3		,	nK		, nWfs, nWfs	)			)
	!allocate(			Fh(		3		,	nK		, nWfs, nWfs	)			)
	!allocate(			Vh(		3		,	nK		, nWfs, nWfs	)			)
	!allocate(			EnH(				nK		,	 nWfs		)			)
	
	!
	write(*,*)"[main]:**************************Infos about this run*************************"
	write(*,*)"[main]: nK=",nK
	write(*,*)"[main]: nQ=",nQ
	write(*,*)"[main]: nG=",nG

    write(*,*)"[main]: nBands=", nBands
	write(*,*)"[main]: nWfs  =", nWfs	
	!
	call cpu_time(T1)
	aT = T1 - T0



	
	!SOLVE ELECTRONIC STRUCTURE & GENERATE THE WANNIER FUNCTIONS ON THE FLY
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************ELECTRONIC STRUCTURE PART*************************"
	call cpu_time(T0)
	!
	!
	call solveHam(unk, En)
	!
	call cpu_time(T1)
	write(*,*)"[main]: done solving Schroedinger eq."
	kT = T1-T0
	



	!PROJECT THE STATES
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************PROJECT STATES *************************"
	call cpu_time(T0)
	!
	!
	call projectUnk(unk, unkW, Uq)
	!
	call cpu_time(T1)
	write(*,*)"[main]: done with projections."
	pT = T1-T0


	
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	call cpu_time(T0)
	if( doWanni ) then
		write(*,*)	"[main]:**************************WANNIER FUNCTION METHOD*************************"
		!
		call wannMethod(unkW, wnf,pWann)
		!
		write(*,*)	"[main]: done with center polarization calc"
	else
		write(*,*)	"[main]: wannier method disabled"
	end if	
	
	call cpu_time(T1)
	
	wT 	= T1 - T0



	!TODO REWRITE WITH EnP, veloP
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	call cpu_time(T0)
	if ( doBerry ) then
		write(*,*)"[main]:**************************WAVEFUNCTION METHOD*************************"
		call berryMethod(unkW, En, Uq, pBerry, pNiu)
		write(*,*)"[main]: done with wavefunction method "
	else
		write(*,*)"[main]: berry method disabled"
	end if
	call cpu_time(T1)
	bT	= T1 - T0
	















	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"[main]:**************************PEIERLS SUB*************************"
	!call cpu_time(T0)
	!!
	!if(doPei)  then
	!	call peierlsMethod(wnf, pPei)
	!end if
	!!
	!call cpu_time(T1)
	!write(*,*)"[main]: done with peierls substitution"
	!peiT = T1 - T0







	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************WRITE OUTPUT*************************"
	call cpu_time(T0)
	!call writeSysInfo() 
	call writeMeshInfo() 
	call writeMeshBin()
	call writeUNKs(unkW)
	write(*,*)"[main]: ...wrote mesh info"
	!call calcIonicPol(pIon)
	pTot	= pIon + pWann
	call writePolFile(pWann, pBerry, pNiu, pPei )
	write(*,*)"[main]: ...wrote polarization txt file"
	!
	call cpu_time(T1)
	oT = T1 - T0
	




	!TIMING INFO SECTION
	call cpu_time(mastT1)
	mastT= mastT1-mastT0
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*) '**************TIMING INFORMATION************************'
	call printTiming(aT, kT, pT, wT, bT,peiT, oT, mastT)
	write(*,*)	"[main]: all done, deallocate and exit"
	!
	!
	stop
end program

