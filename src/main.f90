program main
	!TWO dimensional potential well code
	use omp_lib
	use mathematics, 	only: 		dp, PI_dp

	use sysPara
	use potWellModel, 	only: 		solveHam
	use w90Interface,	only:		w90Interf
	use postW90,		only:		effTBmodel
	use berry,			only:		berryMethod

	use output,		 	only:		writeMeshInfo, writeMeshBin, writeCkASunk, writePolFile,& 
									printTiming, printBasisInfo	!printMat, printInp, printWannInfo,writeSysInfo  



	implicit none

    complex(dp),	allocatable,	dimension(:,:,:)	:: 	ck
    real(dp),		allocatable,	dimension(:,:)		:: 	En    														
    real												:: 	mastT0, mastT1, mastT, T0, T1, &
    															aT,kT,wT,pwT, oT, bT
    
    !timing zero init
    aT		= 0.0
    kT		= 0.0
    wT		= 0.0
    pwT		= 0.0
    bT		= 0.0
    oT 		= 0.0
    mastT	= 0.0

    !INPUT & ALLOCATION SECTION
    call cpu_time(mastT0)
    call cpu_time(T0)
	call readInp()
	!
	!
	write(*,*)"[main]:**************************Infos about this run*************************"
	write(*,*)"[main]: electronic structure mesh nQ=",nQ
	write(*,*)"[main]: interpolation mesh        nK=",nK
	write(*,*)"[main]: basis cutoff parameter  Gcut=",Gcut
	write(*,*)"[main]: basis function   maximum  nG=",nG
	write(*,*)"[main]: only solve for        nSolve=",nSolve
    write(*,*)"[main]: nBands=", nBands
	write(*,*)"[main]: nWfs  =", nWfs
	write(*,*)"[main]: w90 seed_name= ", seedName	
	!
	call cpu_time(T1)
	aT = T1 - T0



	
	!ELECTRONIC STRUCTURE
	if( doSolveHam ) then
		call cpu_time(T0)	
		allocate(	En(						nBands	, 	nQ		)	)
		allocate(	ck(			nG		,	nBands  	,	nQ	)	)
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"[main]:**************************ELECTRONIC STRUCTURE PART*************************"
		write(*,*)	"[main]: start electronic structure calculation now"
		call solveHam(ck, En)
		write(*,*)"[main]: done solving Schroedinger eq."
		call cpu_time(T1)
		kT = T1-T0
		!W90
		call cpu_time(T0)	
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"[main]:**************************WANNIER90 SETUP*************************"
		call w90Interf(ck,En)
		write(*,*)"[main]: done setting up wannier. please execute wannier90 now"
		call cpu_time(T1)
		wT = T1-T0
	end if
	!
	
	
	
	!EFF TB - post w90
	call cpu_time(T0)
	if(	doPw90 ) then
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"[main]:**************************POST WANNIER90 *************************"
		

		write(*,*)	"[main]: start with eff TB model calculations"
		call effTBmodel()
		write(*,*)"[main]: done with effective tight binding calculations"
	end if
	call cpu_time(T1)
	pwT	= T1-T0


	!K SPACE METHOD
	call cpu_time(T0)
	if ( doBerry ) then
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"[main]:**************************WAVEFUNCTION METHOD*************************"
		call berryMethod()
		write(*,*)"[main]: done with wavefunction method "
	else
		write(*,*)"[main]: berry method disabled"
	end if
	call cpu_time(T1)
	bT	= T1 - T0

	
	


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
	call printTiming(aT,kT,wT,pwT,bT,oT,mastT)
	write(*,*)	"[main]: all done, exit"
	!
	!
	stop
end program

