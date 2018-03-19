program main
	!TWO dimensional potential well code
	use mpi
	use util_math, 	only: 		dp, PI_dp, aUtoTesla, aUtoEv, aUtoAngstrm

	use util_sysPara
	use ham_Solver, 	only: 		solveHam
	use pol_Berry,		only:		berryMethod
	use tb_interp,		only:		tb_method

	use util_basisIO,	only:

	use util_output,	only:		writeMeshInfo, writePolFile, write_K_lattices, & 
									printTiming, printBasisInfo	!printMat, printInp, printWannInfo,writeSysInfo  


							

	implicit none
	!#include "mpif.h"
	
														
    real	:: 	mastT0, mastT1, mastT, T0, T1, &
    		    	alloT,hamT,postWT, berryT	

    !MPI INIT
	call MPI_INIT( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, myID, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, nProcs, ierr)
    !
    !SETUP
    if( myID == root ) then
    	alloT	= 0.0
    	hamT	= 0.0
    	postWT	= 0.0
    	berryT	= 0.0
    	mastT	= 0.0
    	!
   		write(*,*)								"[main]:**************************setup Grids*************************"
   		call cpu_time(mastT0)
   		call cpu_time(T0)
    end if

    !read & distribute input
  	call readInp()
	if(myID == root ) then
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"*"
	end if
	!
	!check if equal kpt distribution among mpi procs is possible -> if not break
	if( mod(nQ,nProcs)/=0)  stop 				'[main]: ERROR mpi threads have to be integer fraction of nQ'
	!
	!print info
	if( myID == root  ) then
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"[main]:**************************Infos about this run*************************"
		write(*,'(a,i3,a,i3,a)')				"[main]: q mesh nQx=",nQx," nQy=",nQy
		write(*,'(a,f8.3,a,f8.3)')				"[main]: dqx=",dqx, " dqy=",dqy
		write(*,*)								"[main]: interpolation mesh        nK=",nK
		write(*,*)								"[main]: basis cutoff parameter  Gcut=",Gcut
		write(*,'(a,i7,a,i7,a)')				"[main]: basis function   maximum  nG=",GmaxGLOBAL," of ",nG," trial basis functions"
		write(*,*)								"[main]: only solve for        nSolve=",nSolve
    	write(*,*)								"[main]: nBands=", nBands
		write(*,*)								"[main]: nWfs  =", nWfs
		write(*,*)								"[main]: nAt   =", nAt
		write(*,'(a,i3,a)')						"[main]: project ",nWfs/nAt," states onto each atom"
		write(*,'(a,f6.3,a,f6.3,a,f6.3,a)')		"[main]: Bext = (",Bext(1)*aUtoTesla,", ",Bext(2)*aUtoTesla, ", ",Bext(3)*aUtoTesla,") (T)"
		write(*,*)								"[main]: w90 seed_name= ", seedName	
		
		write(*,*)								"*"
		write(*,*)								"[main]:**************************Hamiltonian info*************************"
		write(*,'(a,i3,a,e12.4,a)')				"[main]:	featuring ",nAt," well potentials, deepest well	",minval(atPot(:)),			" Hartree"
		if( doRashba)	write(*,'(a,e12.4,a)')	"[main]:	featuring a Rashba term with prefact 		",aRashba,							" Hartree a0"
		if( doZeeman)	write(*,'(a,e12.4,a)')	"[main]:	featuring a Zeeman term with prefact 		",0.5_dp*Bext(3),					" Hartree"
		if( doMagHam)	write(*,'(a,e12.4,a)')	"[main]:	featuring a osc. mag. field  prefact 		",0.5_dp*Bext(3)*aX/(2.0_dp*PI_dp),	" Hartree a0"
		write(*,*)								"*"
		if( debugHam)	write(*,*)				"[main]:	debug on: will test Ham for hermiticity"
		write(*,*)								"*"
		write(*,*)								"*"
		!
		call cpu_time(T1)
		alloT = T1 - T0


		!try to print some WARNINGs for to small Gcut
		write(*,*)								"[main]:**************************BASIS SET DEBUG*************************"
		call printBasisInfo()
		write(*,*)								"[main]: ...wrote basis set debug info"
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"*"
	end if

	!HAM SOLVER
	call MPI_BARRIER( MPI_COMM_WORLD, ierr )	
	if( doSolveHam ) then
		!
		if( myID == root )	call cpu_time(T0)
		if( myID == root ) 	write(*,*)			"[main]:**************************ELECTRONIC STRUCTURE RUN*************************"
		!
		!
		call solveHam()
		call MPI_BARRIER( MPI_COMM_WORLD, ierr )
		!
		!finalize:
		if( myID == root ) then
			call write_K_lattices()
			write(*,*)							"[main]:	...wrote k lattices to file"
			call writeMeshInfo() 
			write(*,*)							"[main]: ...wrote meshInfo.txt"
			write(*,*)							"*"
			write(*,*)							"*"
			write(*,*)							"*"
			write(*,*)							"*"
			write(*,*)							"[main]: done solving Schroedinger eq., please execute wannier90 now"
			call cpu_time(T1)
			hamT = T1-T0
		end if
		!
	end if
	




	!POST HAM SOLVER
	if( .not. doSolveHam .and. myID == root ) then	

	
	
		!K SPACE METHOD
		call cpu_time(T0)
		
		if ( doBerry ) then
			write(*,*)							"[main]:**************************BERRY METHOD*************************"
			call berryMethod()
			write(*,*)							"[main]: done with wavefunction method "
			write(*,*)							"*"
			write(*,*)							"*"
			write(*,*)							"*"
			write(*,*)							"*"
		end if
		call cpu_time(T1)
		berryT	= T1 - T0


		!EFF TB - post w90
		call cpu_time(T0)
		if(	doPw90 ) then
			write(*,*)							"[main]:**************************POST WANNIER90 (TB interpolation) *************************"
			call tb_method()
			write(*,*)							"[main]: done with interpolation"
			write(*,*)							"*"
			write(*,*)							"*"
			write(*,*)							"*"
			write(*,*)							"*"
		end if
		call cpu_time(T1)
		postWT	= T1-T0
	
		!
		!TIMING INFO SECTION
		call cpu_time(mastT1)
		mastT= mastT1-mastT0
		write(*,*) 								"[main]:**************TIMING INFORMATION************************"
		call printTiming(alloT,hamT,berryT, postWT,mastT)
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"*"
	end if

	!MPI FINALIZE
	call MPI_Barrier( MPI_COMM_WORLD, ierr )
	call MPI_FINALIZE ( ierr )
	!
	!
	stop
end program

