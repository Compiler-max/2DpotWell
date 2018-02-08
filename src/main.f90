program main
	!TWO dimensional potential well code
	use mpi
	use mathematics, 	only: 		dp, PI_dp

	use sysPara
	use potWellModel, 	only: 		solveHam
	use basisIO,		only:		readAbIn, readBasis
	use w90Interface,	only:		setup_w90
	use postW90,		only:		effTBmodel
	use berry,			only:		berryMethod

	use output,		 	only:		writeMeshInfo, writeMeshBin, writePolFile, write_K_lattices, & 
									printTiming, printBasisInfo	!printMat, printInp, printWannInfo,writeSysInfo  


							

	implicit none
	!#include "mpif.h"
	


    complex(dp),	allocatable,	dimension(:,:,:)	:: 	ck
    real(dp),		allocatable,	dimension(:,:)		:: 	En    														
    real												:: 	mastT0, mastT1, mastT, T0, T1, &
    															alloT,hamT,wannT,postWT, outT, berryT	

    !MPI INIT
	call MPI_INIT( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, myID, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, nProcs, ierr)
    !
    !SETUP
    if( myID == root) then
    	alloT	= 0.0
    	hamT	= 0.0
    	wannT	= 0.0
    	postWT	= 0.0
    	berryT	= 0.0
    	outT 	= 0.0
    	mastT	= 0.0
    	!
   		write(*,*)"[main]:**************************setup Grids*************************"
   		call cpu_time(mastT0)
   		call cpu_time(T0)
    end if

    !read & distribute input
  	call readInp()
	!
	!check if equal kpt distribution among mpi procs is possible -> if not break
	if( mod(nQ,nProcs)/=0)  stop '[main]: ERROR mpi threads have to be integer fraction of nQ'
	!
	!print info
	if( myID == root) then
		write(*,*)					"*"
		write(*,*)					"*"
		write(*,*)					"[main]:**************************Infos about this run*************************"
		write(*,'(a,i3,a,i3,a)')	"[main]: q mesh nQx=",nQx," nQy=",nQy
		write(*,'(a,f8.3,a,f8.3)')	"[main]: dqx=",dqx, " dqy=",dqy
		write(*,*)					"[main]: interpolation mesh        nK=",nK
		write(*,*)					"[main]: basis cutoff parameter  Gcut=",Gcut
		write(*,'(a,i7,a,i7,a)')	"[main]: basis function   maximum  nG=",GmaxGLOBAL," of ",nG," trial basis functions"
		write(*,*)					"[main]: only solve for        nSolve=",nSolve
    	write(*,*)					"[main]: nBands=", nBands
		write(*,*)					"[main]: nWfs  =", nWfs
		write(*,*)					"[main]: w90 seed_name= ", seedName	
		write(*,*)					"*"
		write(*,*)					"*"
		write(*,*)					"*"
		write(*,*)					"*"
		!
		call cpu_time(T1)
		alloT = T1 - T0


		!try to print some WARNINGs for to small Gcut
		write(*,*)"[main]:**************************BASIS SET DEBUG*************************"
		call printBasisInfo()
		write(*,*)"[main]: ...wrote basis set debug info"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		call cpu_time(T1)
		outT = T1 - T0
	end if

	
	
	!HAM SOLVER
	call MPI_BARRIER( MPI_COMM_WORLD, ierr )	
	if( doSolveHam ) then
		!
		if( myID == root )	call cpu_time(T0)
		if( myID == root ) 	write(*,*)"[main]:**************************ELECTRONIC STRUCTURE RUN*************************"
		!
		call solveHam()
		call MPI_BARRIER( MPI_COMM_WORLD, ierr )
		if( myID == root ) then
			write(*,*)"[main]: done solving Schroedinger eq."
			call cpu_time(T1)
			hamT = T1-T0
		end if
		!
	end if
	




	!POST HAM SOLVER
	if( .not. doSolveHam .and. myID == root ) then	
		write(*,*)"[main]:**************************READ E-STRUCTURE*************************"
		allocate(	En(						nSolve	, 	nQ	)	)
		allocate(	ck(			GmaxGLOBAL,	nSolve 	,	nQ	)	)

		!READ IN ELECTRONIC STRUCTURE
		call readAbIn(ck, En)
		!call readBasis() !reads in Gvec, nGq (optiónal)


	

		call cpu_time(T0)
		if( doPrepW90 )	 then
			write(*,*)	"[main]:**************************WANNIER90 INTERFACE*************************"
			call setup_w90(ck,En)
			call write_K_lattices()
			write(*,*)	"[main]: please run w90 now"

		end if
		call cpu_time(T1)
		wannT = T1 - T0




		!EFF TB - post w90
		call cpu_time(T0)
		if(	doPw90 ) then
			write(*,*)"[main]:**************************POST WANNIER90 *************************"
			write(*,*)	"[main]: start with eff TB model calculations"
			call effTBmodel()
			write(*,*)	"[main]: done with effective tight binding calculations"
			write(*,*)"*"
			write(*,*)"*"
			write(*,*)"*"
			write(*,*)"*"
		end if
		call cpu_time(T1)
		postWT	= T1-T0
	
	
		!K SPACE METHOD
		call cpu_time(T0)
		
		if ( doBerry ) then
			write(*,*)"[main]:**************************BERRY METHOD*************************"
			call berryMethod(ck, En)
			write(*,*)"[main]: done with wavefunction method "
			write(*,*)"*"
			write(*,*)"*"
			write(*,*)"*"
			write(*,*)"*"
		end if
		call cpu_time(T1)
		berryT	= T1 - T0

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
		!
		!TIMING INFO SECTION
		call cpu_time(mastT1)
		mastT= mastT1-mastT0
		write(*,*) '**************TIMING INFORMATION************************'
		call printTiming(alloT,hamT,wannT,postWT,berryT,outT,mastT)
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
	end if

	!MPI FINALIZE
	call MPI_Barrier( MPI_COMM_WORLD, ierr )
	write(*,'(a,i3,a)')	"[#",myID,";main]: all done, exit"
	call MPI_FINALIZE ( ierr )
	!
	!
	stop
end program

