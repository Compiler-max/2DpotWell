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

	use output,		 	only:		writeMeshInfo, writeMeshBin, writeCkASunk, writePolFile,& 
									printTiming, printBasisInfo	!printMat, printInp, printWannInfo,writeSysInfo  



	implicit none

    complex(dp),	allocatable,	dimension(:,:,:)	:: 	ck, ckW, Uq
    real(dp),		allocatable,	dimension(:,:)		:: 	En
    real(dp) 											:: 	pWann(2), pBerry(2), pNiu(3), pPei(3)
    real												:: 	mastT0, mastT1, mastT, T0, T1, &
    															aT,kT,wT, oT, bT, peiT, pT
    
    !INPUT & ALLOCATION SECTION
    call cpu_time(mastT0)
    call cpu_time(T0)
	call readInp()
	!
	
	allocate(			En(						nBands	, 	nQ		)				)
	allocate(			ck(			nG		,	nBands  	,	nQ	)				)
	allocate(			ckW(		nG		,	nWfs		,	nQ	)				)
	allocate(			Uq(			nBands	,	nWfs	, 	nQ		)				)
	!
	write(*,*)"[main]:**************************Infos about this run*************************"
	write(*,*)"[main]: electronic structure mesh nQ=",nQ
	write(*,*)"[main]: interpolation mesh        nK=",nK
	write(*,*)"[main]: basis cutoff parameter  Gcut=",Gcut
	write(*,*)"[main]: basis function   maximum  nG=",nG
	write(*,*)"[main]: only solve for        nSolve=",nSolve
	write(*,*)"[main]: real space points per cell  =",nR/nSC
    write(*,*)"[main]: nBands=", nBands
	write(*,*)"[main]: nWfs  =", nWfs	
	!
	call cpu_time(T1)
	aT = T1 - T0



	
	!ELECTRONIC STRUCTURE
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
	



	!PROJECTIONS
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************PROJECT STATES *************************"
	call cpu_time(T0)
	!
	call projectUnk(ck, ckW, Uq)
	!
	call cpu_time(T1)
	write(*,*)"[main]: done with projections."
	pT = T1-T0


	!REAL SPACE METHOD
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



	!K SPACE METHOD
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!call cpu_time(T0)
	!if ( doBerry ) then
	!	write(*,*)"[main]:**************************WAVEFUNCTION METHOD*************************"
	!	call berryMethod(ckW, En, Uq, pBerry, pNiu, pPei)
	!	write(*,*)"[main]: done with wavefunction method "
	!	write(*,'(a,f12.8,a,f12.8,a)')	"[main]: calculated zero order pBerry=(",pBerry(1),", ",pBerry(2),")."
	!else
	!	write(*,*)"[main]: berry method disabled"
	!end if
	!call cpu_time(T1)
	!bT	= T1 - T0
	


	!OUTPUT
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************WRITE OUTPUT*************************"
	call cpu_time(T0)
	!
	call writePolFile(pWann, pBerry, pNiu, pPei )
	write(*,*)"[main]: ...wrote polarization txt file"
	call writeMeshInfo() 
	write(*,*)"[main]: ...wrote mesh info"
	if( writeBin )	then
		call writeMeshBin()
		write(*,*)"[main]: ...wrote mesh bin"
		!call writeUNKs(unkW)
		call writeCkASunk(ck, ckW)
		write(*,*)"[main]: ...wrote binary files for meshes and unks"
	end if
	!
	call cpu_time(T1)
	oT = T1 - T0
	
	
	!WARNINGS IF GCUT IS TO HIGH
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"[main]:**************************BASIS SET DEBUG*************************"
	!call printBasisInfo()
	!write(*,*)"[main]: ...wrote basis set debug info"


	
	!TIMING INFO SECTION
	call cpu_time(mastT1)
	mastT= mastT1-mastT0
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*) '**************TIMING INFORMATION************************'
	call printTiming(aT, kT, pT, wT, bT,peiT, oT, mastT)
	write(*,*)	"[main]: all done, exit"
	!
	!
	stop
end program

