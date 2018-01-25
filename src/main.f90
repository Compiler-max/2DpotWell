program main
	!TWO dimensional potential well code
	use mpi
	use mathematics, 	only: 		dp, PI_dp

	use sysPara
	use potWellModel, 	only: 		solveHam
	use w90Interface,	only:		w90Interf
	use postW90,		only:		effTBmodel
	use berry,			only:		berryMethod

	use output,		 	only:		writeMeshInfo, writeMeshBin, writePolFile,& 
									printTiming, printBasisInfo	!printMat, printInp, printWannInfo,writeSysInfo  


							

	implicit none

	


    complex(dp),	allocatable,	dimension(:,:,:)	:: 	ck
    real(dp),		allocatable,	dimension(:,:)		:: 	En    														
    real												:: 	mastT0, mastT1, mastT, T0, T1, &
    															aT,kT,wT,pwT, oT, bT	
    logical												::	mpiSuccess		

    !MPI INIT
	call MPI_INIT( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, myID, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, nProcs, ierr)
    root = 0
    mpiSuccess = .true.
    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    !
    !
    if( myID == root) then
    	aT		= 0.0
    	kT		= 0.0
    	wT		= 0.0
    	pwT		= 0.0
    	bT		= 0.0
    	oT 		= 0.0
    	mastT	= 0.0
    	!
   		write(*,*)"[main]:**************************setup Grids*************************"
   		call cpu_time(mastT0)
   		call cpu_time(T0)
    end if

    !READ INPUT FILE & DISTRIBUTE
  	call readInp()
	!
	!CHECK IF QPTS CAN BE EQUALLY DISTRIBUTED -> if not break
	if( mod(nQ,nProcs)/=0)  then
		if(myID == root) write(*,*)"[main]: CRITICAL WARNING: mpi threads have to be integer fraction of nQ"
		mpiSuccess = .false.
	end if
	!
	!PRINT SOME INFO
	if( myID == root) then
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"[main]:**************************Infos about this run*************************"
		write(*,*)"[main]: electronic structure mesh nQ=",nQ
		write(*,*)"[main]: interpolation mesh        nK=",nK
		write(*,*)"[main]: basis cutoff parameter  Gcut=",Gcut
		write(*,*)"[main]: basis function   maximum  nG=",GmaxGLOBAL," of ",nG," trial basis functions"
		write(*,*)"[main]: only solve for        nSolve=",nSolve
    	write(*,*)"[main]: nBands=", nBands
		write(*,*)"[main]: nWfs  =", nWfs
		write(*,*)"[main]: w90 seed_name= ", seedName	
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		!
		call cpu_time(T1)
		aT = T1 - T0
	end if
	
		
	
	
	!ELECTRONIC STRUCTURE
	if( mpiSuccess .and. doSolveHam ) then
		!call cpu_time(T0)	
		allocate(	En(						nSolve	, 	qChunk	)	)
		allocate(	ck(			nG		,	nSolve 	,	qChunk	)	)
		if( myID == root ) 	write(*,*)"[main]:**************************ELECTRONIC STRUCTURE PART*************************"
		call MPI_BARRIER( MPI_COMM_WORLD, ierr )	
		!
		!
		call solveHam(ck, En)
		if( myID == root ) then
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
		end if
		call w90Interf( ck,En)
		if( myID == root ) then
			write(*,*)"[main]: done setting up wannier. please execute wannier90 now"
			write(*,*)"*"
			write(*,*)"*"
			write(*,*)"*"
			write(*,*)"*"
			call cpu_time(T1)
			wT = T1-T0
		end if
	end if
	call MPI_BARRIER( MPI_COMM_WORLD, ierr )




	if( myID == root ) then		
		!EFF TB - post w90
		call cpu_time(T0)
		write(*,*)"[main]:**************************POST WANNIER90 *************************"
		if(	doPw90 ) then
			
			write(*,*)	"[main]: start with eff TB model calculations"
			call effTBmodel()
			write(*,*)	"[main]: done with effective tight binding calculations"
		else
			write(*,*)	"[main]: effective TB model disabled"
		end if
		!
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		call cpu_time(T1)
		pwT	= T1-T0
	
	
		!K SPACE METHOD
		call cpu_time(T0)
		write(*,*)"[main]:**************************WAVEFUNCTION METHOD*************************"
		if ( doBerry ) then
			call berryMethod()
			write(*,*)"[main]: done with wavefunction method "
		else
			write(*,*)"[main]: berry method disabled"
		end if
		!
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		call cpu_time(T1)
		bT	= T1 - T0

		!OUTPUT
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
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		call cpu_time(T1)
		oT = T1 - T0
		
		
		!WARNINGS IF GCUT IS TO HIGH
		write(*,*)"[main]:**************************BASIS SET DEBUG*************************"
		call printBasisInfo()
		write(*,*)"[main]: ...wrote basis set debug info"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
	
	
	
		!TIMING INFO SECTION
		call cpu_time(mastT1)
		mastT= mastT1-mastT0
		write(*,*) '**************TIMING INFORMATION************************'
		call printTiming(aT,kT,wT,pwT,bT,oT,mastT)
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
	end if

	call MPI_Barrier( MPI_COMM_WORLD, ierr )
	write(*,*)	"[",myID,";main]: all done, exit"


	call MPI_FINALIZE ( ierr )
	!
	!
	stop
end program

