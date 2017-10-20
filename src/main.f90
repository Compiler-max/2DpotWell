program main
	!TWO dimensional potential well code
	use omp_lib
	use mathematics, 	only: 		dp, PI_dp

	
	use sysPara
	use input,			only:		filesExist, readHam 
	use potWellModel, 	only: 		solveHam
	use projection,		only:		projectUnk
	use wannier,	 	only: 		wannMethod	
	use berry,			only:		berryMethod

	use polarization,	only:		calcIonicPol
	use peierls,		only:		peierlsMethod
	use output,		 	only:		writeMeshInfo, writeMeshBin, writeUNKs, writePolFile,& 
									printTiming	!printMat, printInp, printWannInfo,writeSysInfo  



	implicit none

	
    complex(dp),	allocatable,	dimension(:,:,:)	:: 	ck, ckW, Uq
    real(dp),		allocatable,	dimension(:,:)		:: 	En
    real(dp) 											:: 	pWann(2), pBerry(2), pNiu(3), pPei(3)
    real												:: 	mastT0, mastT1, mastT, T0, T1, &
    															aT,kT,wT, oT, bT, peiT, pT
    call cpu_time(mastT0)




    !INPUT & ALLOCATION SECTION
    call cpu_time(T0)
	call readInp()
	!electronic structure arrays
	allocate(			En(						nBands	, 	nQ		)				)
	allocate(			ck(			nG		,	nBands  	,	nQ	)				)
	allocate(			ckW(		nG		,	nWfs		,	nQ	)				)
	allocate(			Uq(			nBands	,	nWfs	, 	nQ		)				)
	
	
	!
	write(*,*)"[main]:**************************Infos about this run*************************"
	write(*,*)"[main]: electronic structure mesh nQ=",nQ
	write(*,*)"[main]: interpolation mesh        nK=",nK
	write(*,*)"[main]: basis function per dim nGdim=",nGdim
	write(*,*)"[main]: basis functions           nG=",nG
	write(*,*)"[main]: real space points per cell  =",nR/nSC
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
	if( .not. doSolveHam .and. filesExist() ) then
		write(*,*)	"[main]: electronic structure disabled. Read in unks and energies"
		call readHam( ck, En)
	else
		write(*,*)	"[main]: start electronic structure calculation now"
		call solveHam(ck, En)
	end if
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
	call projectUnk(ck, ckW, Uq)
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
		call wannMethod(ckW, pWann)
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
		call berryMethod(ckW, En, Uq, pBerry, pNiu, pPei)
		write(*,*)"[main]: done with wavefunction method "
	else
		write(*,*)"[main]: berry method disabled"
	end if
	call cpu_time(T1)
	bT	= T1 - T0
	



	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************WRITE OUTPUT*************************"
	call cpu_time(T0)
	call writeMeshInfo() 
	
	if( writeBin )	then
		call writeMeshBin()
		!call writeUNKs(unkW)
	end if
	write(*,*)"[main]: ...wrote mesh info"
	call writePolFile(pWann, pBerry, pNiu, pPei )
	write(*,*)"[main]: ...wrote polarization txt file"
	!
	call cpu_time(T1)
	oT = T1 - T0
	
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************BASIS SET INFO*************************"
	if(			nG		< 		vol * Gcut * dsqrt(Gcut) /	(2.0_dp*PI_dp**2)		) then
		write(*,*)	"[main]: increase nG or decrease Gcut"
	end if
	if(		nR	<		vol * dsqrt(Gcut) 	/ PI_dp)	write(*,*)"[main]: need more real space points or smaller Gcut"

	if(		nRx / nSCx <  nQx) write(*,*)"[main]: need more reals space points or less k points"

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

