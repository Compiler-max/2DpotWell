program main
	!TWO dimensional potential well code
	use omp_lib
	use mathematics, 	only: 		dp, PI_dp

	use sysPara
	use input,			only:		filesExist, readHam 
	use potWellModel, 	only: 		solveHam
	use w90Interface,	only:		w90Interf
	use postW90,		only:		effTBmodel
	use berry,			only:		berryMethod

	use output,		 	only:		writeMeshInfo, writeMeshBin, writeCkASunk, writePolFile,& 
									printTiming, printBasisInfo	!printMat, printInp, printWannInfo,writeSysInfo  



	implicit none

    complex(dp),	allocatable,	dimension(:,:,:)	:: 	ck, ckW, Uq
    real(dp),		allocatable,	dimension(:,:)		:: 	En
    real(dp) 											:: 	pWann(2), pBerry(2), pNiuF2(3), pNiuF3(3), pPei(3), &
    														pTBwann(3), pTBconn(3), pTBniuF2(3), pTBniuF3(3), pTBpei(3)
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
	write(*,*)"[main]: w90 seed_name= ", seedName	
	!
	call cpu_time(T1)
	aT = T1 - T0



	
	!ELECTRONIC STRUCTURE
	call cpu_time(T0)
	!
	if( doSolveHam ) then
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"[main]:**************************ELECTRONIC STRUCTURE PART*************************"
		write(*,*)	"[main]: start electronic structure calculation now"
		call solveHam(ck, En)
		write(*,*)"[main]: done solving Schroedinger eq."
		!W90
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"[main]:**************************WANNIER90 SETUP*************************"
		call w90Interf(ck,En)
		write(*,*)"[main]: done setting up wannier. please execute wannier90 now"
	end if
	!
	call cpu_time(T1)
	kT = T1-T0
	
	
	!EFF TB - post w90
	if(	doPw90 ) then
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"[main]:**************************POST WANNIER90 *************************"
		

		write(*,*)	"[main]: start with eff TB model calculations"
		call effTBmodel(pTBwann, pTBconn, pTBniuF2, pTBniuF3, pTBpei)
		write(*,*)"[main]: done with effective tight binding calculations"
	end if


	!PROJECTIONS
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"[main]:**************************PROJECT STATES *************************"
	!call cpu_time(T0)
	!!
	!call projectUnk(ck, ckW, Uq)
	!!
	!call cpu_time(T1)
	!write(*,*)"[main]: done with projections."
	!pT = T1-T0


	!K SPACE METHOD
	call cpu_time(T0)
	if ( doBerry ) then
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"[main]:**************************WAVEFUNCTION METHOD*************************"
		call berryMethod(pBerry, pNiuF2, pNiuF3, pPei)
		write(*,*)"[main]: done with wavefunction method "
		write(*,'(a,f12.8,a,f12.8,a)')	"[main]: calculated zero order pBerry=(",pBerry(1),", ",pBerry(2),")."
	else
		write(*,*)"[main]: berry method disabled"
	end if
	call cpu_time(T1)
	bT	= T1 - T0



	!REAL SPACE METHOD
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!write(*,*)"*"
	!call cpu_time(T0)
	!if( doWanni ) then
	!	write(*,*)	"[main]:**************************WANNIER FUNCTION METHOD*************************"
	!	!
	!	call wannMethod(ckW, pWann)
	!	!
	!	write(*,*)	"[main]: done with center polarization calc"
	!else
	!	write(*,*)	"[main]: wannier method disabled"
	!end if	
	!call cpu_time(T1)
	!wT 	= T1 - T0



	
	


	!OUTPUT
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************WRITE OUTPUT*************************"
	call cpu_time(T0)
	!
	call writeMeshInfo() 
	write(*,*)"[main]: ...wrote mesh info"
	if( writeBin )	then
		call writeMeshBin()
		write(*,*)"[main]: ...wrote mesh bin"
		write(*,*)"[main]: ...wrote binary files for meshes and unks"
	end if
	!
	call cpu_time(T1)
	oT = T1 - T0
	
	
	!WARNINGS IF GCUT IS TO HIGH
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"*"
	write(*,*)"[main]:**************************BASIS SET DEBUG*************************"
	call printBasisInfo()
	write(*,*)"[main]: ...wrote basis set debug info"


	
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

