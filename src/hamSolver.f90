module ham_Solver
	!this modules sets up the Hamiltonian matrix and solves it for each k point
	!	in the process the wannier functions are generated aswell with routines from wannGen module
	use mpi
	use omp_lib
	use util_math,		only:	dp, PI_dp,i_dp, machineP, myExp, &
								eigSolverPART, isUnit, isHermitian, aUtoEv
	use util_sysPara				
	
	use	ham_PotWell,	only:	add_potWell
	use	ham_Zeeman,		only:	add_Zeeman
	use ham_Magnetic,	only:	add_magHam



	use ham_PWbasis,	only:	calcVeloGrad, calcAmatANA, calcMmat
	use util_basisIO,	only:	writeABiN_basVect, writeABiN_energy, writeABiN_basCoeff, writeABiN_velo, writeABiN_Amn, writeABiN_Mmn, &
								read_coeff, read_gVec
	use util_w90Interf,	only:	setup_w90, write_w90_matrices, printNNinfo
	implicit none	
	!#include "mpif.h"

	private
	public ::			solveHam			



		

	


	contains
!public:
	subroutine solveHam()   !call solveHam(wnF, unk, EnW, VeloBwf)
		!solves Hamiltonian at each k point
		!also generates the Wannier functions on the fly (sum over all k)
		!																

		integer							::	nntotMax, nntot
		integer,		allocatable		::	nnlist(:,:), nncell(:,:,:)
		!	
		!
		nntotMax = 12
		allocate(	nnlist(		nQ, nntotMax)		)
		allocate(	nncell(3, 	nQ, nntotMax)		)
		!
		

		!SOLVE HAM at each k-point
		call worker()
		!
		!get FD scheme from w90
		if( myID == root )  then
			call setup_w90(nntot, nnlist, nncell)
			call printNNinfo(nntot, nnlist)
		end if
		!
		!Bcast FD scheme
		call MPI_Bcast(nntot,			1		, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
		call MPI_Bcast(nnlist,		nntotMax*nQ	, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
		call MPI_Bcast(nncell,	3*	nntotMax*nQ	, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
		if( myID == root) write(*,'(a,i3,a)')	"[#",myID,";potWellMethod]: broadcasted fd scheme"
		!
		!calc Mmat
		call calc_Mmat(nntot, nnlist, nncell)
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)	!without root may try to read Mmat files which are not written


		!call w90 interface to write input files & .mmn etc. files

		if( myID == root ) then
			write(*,*)	"*"
			write(*,*)	"*"
			write(*,*)	"*"
			write(*,'(a,i3,a)')		"[#",myID,";solveHam]: wrote Mmn files, now collect files to write wannier90 input files"
			call write_w90_matrices()
			write(*,'(a,i3,a)')		"[#",myID,";solveHam]: wrote w90 matrix input files (.amn, .mmn, .eig)"
		end if


		!
		return
	end subroutine




!private:
	subroutine worker()
		!			solve Ham, write results and derived quantites														
		complex(dp),	allocatable		::	Hmat(:,:) , ck_temp(:,:), Amn_temp(:,:), velo_temp(:,:,:)
		real(dp),		allocatable		::	En_temp(:)
		integer							:: 	qi, qLoc, found, Gsize, boundStates
		!	
		!
		allocate(	Hmat(				Gmax,	Gmax				)	)
		allocate(	ck_temp(		GmaxGLOBAL, nSolve				)	)
		allocate(	En_temp(				Gmax					)	)
		allocate(	velo_temp(	3, 	nSolve,		nSolve				)	)	
		allocate(	Amn_temp(		nBands,		nWfs				)	)
		!
		qLoc = 1
		call MPI_BARRIER( MPI_COMM_WORLD, ierr)
		do qi = myID*qChunk +1, myID*qChunk + qChunk
			!
			!!SETUP HAM
			call populateH(qLoc, Hmat) 
			!
			!SOLVE HAM
			ck_temp	= dcmplx(0.0_dp)
			En_temp	= 0.0_dp
			Gsize 	= nGq(qLoc)
			call eigSolverPART(Hmat(1:Gsize,1:Gsize),En_temp(1:Gsize), ck_temp(1:Gsize,:), found)
			!
			!get derived quantities
			call calcAmatANA(qLoc, ck_temp, Amn_temp)
			call calcVeloGrad(qLoc, ck_temp, velo_temp)

			!
			!DEBUG TESTS
			if( found /= nSolve )	write(*,'(a,i3,a,i5,a,i5)'	)	"[#",myID,";solveHam]:WARNING only found ",found," bands of required ",nSolve
			if( nBands > found	)	stop	"[solveHam]: ERROR did not find required amount of bands"
			if( Gsize < nSolve	) 	stop	"[solveHam]: cutoff to small to get bands! only get basis functions"
			if( Gsize > Gmax	)	stop	"[solveHam]: critical error in solveHam please contact developer. (nobody but dev will ever read this^^)"
			!
			!WRITE COEFF TO FILE
			call writeABiN_basVect(qi, Gvec(:,:,qLoc))
			call writeABiN_basCoeff(qi, ck_temp)
			call writeABiN_energy(qi, En_temp(1:nSolve))
			call writeABiN_Amn(qi, Amn_temp)
			call writeABiN_velo(qi, velo_temp)

			!FINALIZE
			boundStates = 0
			do while (	En_temp(boundStates+1) < 0.0_dp )
				boundStates = boundStates + 1
			end do
			write(*,'(a,i3,a,i5,a,f6.2,a,a,i5,a,i5,a,i5,a)')"[#",myID,", solveHam]: qi=",qi," lowest energy=",En_temp(1)*aUtoEv,"[eV];",&
														" found #",boundStates," bound states. done tasks=(",qLoc,"/",qChunk,")"
			qLoc = qLoc + 1		
		end do
		!
		!
		return
	end subroutine


	subroutine calc_Mmat(nntot, nnlist, nncell)
		integer,		intent(in)		::	nntot, nnlist(:,:), nncell(:,:,:)
		integer,		allocatable		::	nGq_glob(:)
		complex(dp),	allocatable		::	ck_qi(:,:), cK_nn(:,:), Mmn(:,:,:)
		real(dp),		allocatable		::	Gvec_qi(:,:), Gvec_nn(:,:)
		real(dp)						::	gShift(2)
		integer							::	qi, nn, q_nn, nG_qi, nG_nn, qLoc
		!
		allocate(	nGq_glob(nQ)					)
		allocate(	Gvec_qi(dim, nG)				)
		allocate(	Gvec_nn(dim, nG)				)
		allocate(	ck_nn(GmaxGLOBAL, nSolve)		)
		allocate(	ck_qi(GmaxGLOBAL, nSolve)		)
		allocate(	Mmn(nBands, nBands, nntot)		)
		!
		call MPI_ALLGATHER( nGq	, qChunk, MPI_INTEGER, nGq_glob		, qChunk, MPI_INTEGER, MPI_COMM_WORLD, ierr)
		!
		!
		qLoc = 1
		do qi = myID*qChunk +1, myID*qChunk + qChunk
			!
			do nn = 1, nntot
				!calc overlap of unks
				if( nncell(3,qi,nn)/= 0 ) stop '[w90prepMmat]: out of plane nearest neighbour found. '
				!
				q_nn		=	nnlist(qi,nn)
				nG_qi		= 	nGq_glob(qi)
				nG_nn		= 	nGq_glob(q_nn)
				gShift(1)	= 	real(nncell(1,qi,nn),dp) * 2.0_dp * PI_dp / aX
				gShift(2)	= 	real(nncell(2,qi,nn),dp) * 2.0_dp * PI_dp / aY
				!
				!read basis coefficients
				call read_coeff(qi,	ck_qi)
				call read_coeff(q_nn, ck_nn)
				!
				!read Gvec
				call read_gVec(qi, 		Gvec_qi)
				call read_gVec(q_nn,	Gvec_nn)
				!
				!
				call calcMmat(qi, q_nn, gShift, nG_qi, nG_nn, Gvec_qi, Gvec_nn, ck_qi, ck_nn, Mmn(:,:,nn)	)
			end do
			!
			!write result to file, alternative: gather on root and write .mmn directly
			call writeABiN_Mmn(qi, Mmn)
			write(*,'(a,i3,a,i5,a,i5,a,i5,a)')		"[#",myID,",calc_Mmat]: wrote M_matrix for qi=",qi,"; task (",qLoc,"/",qChunk,") done"
			qLoc = qLoc + 1
		end do
		!		
		!
		return
	end subroutine
















!POPULATION OF H MATRIX
	subroutine populateH(qLoc, Hmat)
		!populates the Hamiltonian matrix by adding 
		!	1. kinetic energy terms (onSite)
		!	2. potential terms V(qi,i,j)
		!and checks if the resulting matrix is hermitian( for debugging)
		integer,		intent(in)	:: 	qLoc
		complex(dp), intent(inout) 	:: 	Hmat(:,:)
		integer						::	i
		!init to zero
		Hmat = dcmplx(0.0_dp) 
		!
		!
		!KINETIC ENERGY
		do i = 1, size(Hmat,1)
			Hmat(i,i) = 0.5_dp * dot_product(	Gvec(:,i,qLoc),	Gvec(:,i,qLoc)	)
		end do
		!
		!POTENTIAL WELLS
		call add_potWell(qLoc, Hmat)
		!
		!ZEEMAN LIKE
		call add_Zeeman(qLoc, Hmat)

		!ADD PEIERLS
		if( doMagHam )	 call add_magHam( qLoc, Hmat)
		!
		!DEBUG
		if(debugHam) then
			if ( .not.	isHermitian(Hmat)	) 	write(*,'(a,i3,a,i3)')	"[#",myID,";populateH]: WARNING Hamiltonian matrix is not Hermitian at qLoc=",qLoc
		end if
		!
		!		
		return
	end subroutine












end module ham_Solver
