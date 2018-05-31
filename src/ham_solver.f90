module ham_Solver
	!this modules sets up the Hamiltonian matrix and solves it for each k point
	!	in the process the wannier functions are generated aswell with routines from wannGen module
	use mpi
	use util_math,		only:	dp, PI_dp,i_dp, machineP, myExp, 							&
								eigSolverPART, isUnit, isHermitian, aUtoEv
	use util_sysPara,	only:	myID, root, ierr, 											&
								dim, aX, aY, 												&
								nQ, qChunk, nGq, GmaxGLOBAL, Gmax, nSolve, nBands, nWfs,	&
								glob_min_gap, glob_min_gap_indirect,						&
								do_w90plot,	doZeeman, doMagHam, doRashba, debugHam
	
	use	ham_PotWell,	only:	add_potWell
	use	ham_Zeeman,		only:	add_Zeeman
	use ham_Magnetic,	only:	add_magHam
	use ham_Rashba,		only:	add_rashba
	use ham_PWbasis,	only:	writeUNKs



	use ham_PWbasis,	only:	calcVeloGrad, calcAmatANA, calcMmat
	use util_basisIO,	only:	writeABiN_energy, writeABiN_basCoeff, writeABiN_Mmn,		 &
								read_coeff, read_gVec, read_energies
	use util_w90Interf,	only:	setup_w90, write_w90_matrices, printNNinfo
	use util_output,	only:	writeEnTXT
	implicit none	
	!#include "mpif.h"

	private
	public ::			solveHam			



		

	


	contains
!public:
	subroutine solveHam()   !call solveHam(wnF, unk, EnW, VeloBwf)
		!solves Hamiltonian at each k point
		!	and prepares the wannier90 input files (.mmn, .amn, .eig, .win, ...)
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
		!
		!
		!call w90 interface to write input files & .mmn etc. files
		if( myID == root ) then
			call prepare_w90_run()
		end if
		!
		!
		return
	end subroutine




!private:
	subroutine worker()
		!			solve Ham, write results and derived quantites														
		complex(dp),	allocatable		::	Hmat(:,:) , ck_temp(:,:)
		real(dp),		allocatable		::	En_temp(:), Gvec_qi(:,:)
		integer							:: 	qi_glob, qi_loc, found, Gsize, boundStates, loc_minBound, loc_maxBound
		real(dp)						::	loc_min_gap, max_valence, min_conduct
		!	
		!
		allocate(	Gvec_qi(	dim,	GmaxGLOBAL					)	)
		allocate(	Hmat(				Gmax,	Gmax				)	)
		allocate(	ck_temp(		GmaxGLOBAL, nSolve				)	)
		allocate(	En_temp(				Gmax					)	)
		!
		qi_loc 	= 1
		loc_minBound = nSolve
		loc_maxBound = -1
		call MPI_BARRIER( MPI_COMM_WORLD, ierr)
		!LOOP KPTS
		do qi_glob = myID*qChunk +1, myID*qChunk + qChunk
			!grep # basis func
			Gsize 	= nGq(qi_loc)
			!
			!read current Gvec
			call read_gVec(qi_glob, 		Gvec_qi)
			!
			!SETUP HAM
			call populateH(Gsize,Gvec_qi, Hmat) 
			!
			!SOLVE HAM
			ck_temp	= dcmplx(0.0_dp)
			En_temp	= 0.0_dp
			call eigSolverPART(Hmat(1:Gsize,1:Gsize),En_temp(1:Gsize), ck_temp(1:Gsize,:), found)
			!
			!WRITE FILEs
			!call writeABiN_basVect(qi_glob, Gvec(:,:,qi_loc))
			call writeABiN_basCoeff(qi_glob, ck_temp)
			call writeABiN_energy(qi_glob, En_temp(1:nSolve))
			!get derived quantities & write to file
			call calcAmatANA(	Gsize, qi_glob, Gvec_qi, ck_temp)
			call calcVeloGrad(	Gsize, qi_glob, Gvec_qi, ck_temp)
			!w90 plot preparation
			if( do_w90plot ) call writeUNKs(qi_glob, Gsize, ck_temp(:,1:nBands), Gvec_qi(:,:))
			!
			!
			!count insulating states
			!boundStates = count_negative_eigenvalues(loc_minBound, loc_maxBound, En_temp)
			!
			!GET BANDGAP
			if( qi_loc==1 ) then
				loc_min_gap = En_temp(Gmax)-En_temp(1)
				max_valence	= En_temp(nWfs)
				min_conduct	= En_temp(nWfs+1)
			else
				loc_min_Gap = min(loc_min_Gap, get_band_gap(qi_loc, qi_glob, found, Gsize, nWfs, En_temp) )
				max_valence	= max(max_valence, En_temp(nWfs)	)
				min_conduct	= min(min_conduct, En_temp(nWfs+1)	)
			end if
			!
			!goto next qpt
			qi_loc = qi_loc + 1		
		end do
		!
		!ANALYZE BANDS
		call band_analyzer(loc_min_Gap, max_valence, min_conduct)
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
		integer							::	qi, nn, q_nn, nG_qi, nG_nn, qi_loc
		!
		allocate(	nGq_glob(nQ)					)
		allocate(	Gvec_qi(dim, GmaxGLOBAL)				)
		allocate(	Gvec_nn(dim, GmaxGLOBAL)				)
		allocate(	ck_nn(GmaxGLOBAL, nSolve)		)
		allocate(	ck_qi(GmaxGLOBAL, nSolve)		)
		allocate(	Mmn(nBands, nBands, nntot)		)
		!
		call MPI_ALLGATHER( nGq	, qChunk, MPI_INTEGER, nGq_glob		, qChunk, MPI_INTEGER, MPI_COMM_WORLD, ierr)
		!
		!
		qi_loc = 1
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
			!write(*,'(a,i3,a,i5,a,i5,a,i5,a)')		"[#",myID,",calc_Mmat]: wrote M_matrix for qi=",qi,"; task (",qi_loc,"/",qChunk,") done"
			qi_loc = qi_loc + 1
		end do
		!		
		!
		return
	end subroutine


	subroutine prepare_w90_run()
		real(dp),		allocatable		::	EnQ(:,:)
		!
		allocate(	EnQ(nSolve,nQ)	)
		!
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,'(a,i3,a)')		"[#",myID,";prepare_w90_run]: wrote Mmn files, now collect files to write wannier90 input files"
		call write_w90_matrices()
		write(*,'(a,i3,a)')		"[#",myID,";prepare_w90_run]: wrote w90 matrix input files (.amn, .mmn, .eig)"
		call read_energies(EnQ)
		call writeEnTXT(EnQ)
		write(*,'(a,i3,a)')		"[#",myID,";prepare_w90_run]: wrote energies to txt file"
		!
		!
		return
	end subroutine













!POPULATION OF H MATRIX
	subroutine populateH(nG_qi, Gvec, Hmat)
		!populates the Hamiltonian matrix by adding 
		!	1. kinetic energy terms (onSite)
		!	2. potential terms V(qi,i,j)
		!and checks if the resulting matrix is hermitian( for debugging)
		integer,		intent(in)	:: 	nG_qi
		real(dp),		intent(in)	::	Gvec(:,:)
		complex(dp), intent(inout) 	:: 	Hmat(:,:)
		integer						::	gi
		!init to zero
		Hmat = dcmplx(0.0_dp) 
		!
		!
		!KINETIC ENERGY
		do gi = 1, nG_qi
			Hmat(gi,gi) = 0.5_dp * dot_product(	Gvec(:,gi),	Gvec(:,gi)	)
		end do
		!
		!POTENTIAL WELLS
		call add_potWell(nG_qi, Gvec, Hmat)
		!
		!ZEEMAN
		if(	doZeeman )		call add_Zeeman(nG_qi, Gvec, Hmat)
		!
		!PEIERLS
		if( doMagHam )	 	call add_magHam(nG_qi, Gvec, Hmat)
		!
		!RASHBA	
		if( doRashba )		call add_rashba(nG_qi, Gvec, Hmat)

		!
		!DEBUG
		if(debugHam) then
			if ( .not.	isHermitian(Hmat)	) 	write(*,'(a,i3,a)')	"[#",myID,";populateH]: WARNING Hamiltonian matrix is not Hermitian "
		end if
		!
		!		
		return
	end subroutine



!BAND STRUCTURE TESTS
	integer function count_negative_eigenvalues(loc_minBound, loc_maxBound, energies)
		integer,		intent(inout)		::	loc_minBound, loc_maxBound
		real(dp),		intent(in)			::	energies(:)
		integer								::	boundStates
		!
		!count boundStates
		boundStates = 0
		do while (	energies(boundStates+1) < 0.0_dp .and. boundStates < size(energies)	)
			boundStates = boundStates + 1
		end do
		loc_minBound = min(loc_minBound,boundStates)
		loc_maxBound = max(loc_maxBound,boundStates)
		!
		!return
		count_negative_eigenvalues = boundStates
		return
	end function


	subroutine band_analyzer(loc_minBound, loc_max_valence, loc_min_conduct)
		!band_analyzer(loc_min_Gap, max_valence, min_conduct)
		!
		!	test if all kpts have same amount of negative energy eigenvalues (i.e. no states become metallic in parts of BZ)
		!
		real(dp),		intent(in)		::	loc_min_gap, loc_max_valence, loc_min_conduct
		real(dp),		intent(in)		::	glob_max_valence, glob_min_conduct
		integer							::	glob_minBound	,	glob_maxBound 
		!
		call MPI_REDUCE(loc_max_valence	,	glob_max_valence		, 	1, 	MPI_DOUBLE_PRECISION, 		MPI_MAX,	root, MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(loc_min_conduct	, 	glob_min_conduct		, 	1, 	MPI_DOUBLE_PRECISION, 		MPI_MIN,	root, MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(loc_min_gap		,		glob_min_gap		,	1,	MPI_DOUBLE_PRECISION, 		MPI_MIN,	root, MPI_COMM_WORLD, ierr)
		!
		glob_min_gap_indirect = glob_min_conduct - glob_max_valence
		!
		if( myID == root ) then
			write(*,*)	"*"
			write(*,*)	"*"
			write(*,*)	"*"
			write(*,'(a,i3,a)')	"[#",myID,", solveHam]:********band structure analysis:****************************************"
			!
			!TEST FOR METALLIC BANDS
			write(*,'(a,i3,a,f12.4,a)')	"[#",myID,", solveHam]: minimum direct band gap   =",	glob_min_gap*aUtoEv				," (eV)	"
			write(*,'(a,i3,a,f12.4,a)')	"[#",myID,", solveHam]: minimum indirect band gap =",	glob_min_gap_indirect*aUtoEv		," (eV)	"
			!
			!
			!
			!METALLIC WARNINGS
			if( glob_min_gap <= 0.0_dp	) then
				write(*,'(a,i3,a)')	"[#",myID,", solveHam]: WARNING: direct bandgap does not exist"
			end if
			!
			if( glob_min_gap_indirect <= 0.0_dp	) then
				write(*,'(a,i3,a)')	"[#",myID,", solveHam]: WARNING: indirect bandgap does not exist"
			end if
			write(*,*)	"*"
		end if
		!
		!
		return
	end subroutine



	real(dp) function get_band_gap(qi_loc, qi_glob , found, Gsize, boundStates, En_temp)
		integer,		intent(in)		::	qi_loc, qi_glob, found, Gsize, boundStates
		real(dp),		intent(in)		::	En_temp(:)
		real(dp)						::	bandGap
		!
		if( found /= nSolve )	write(*,'(a,i3,a,i5,a,i5)'	)	"[#",myID,";solveHam]:WARNING only found ",found," bands of required ",nSolve
		if( nBands > found	)	stop	"[solveHam]: ERROR did not find required amount of bands"
		if( Gsize < nSolve	) 	stop	"[solveHam]: cutoff to small to get bands! only get basis functions"
		if( Gsize > Gmax	)	stop	"[solveHam]: critical error in solveHam please contact developer. (nobody but dev will ever read this^^)"
		!
		bandGap =	( En_temp(boundStates+1)-En_temp(boundStates)  ) 
		!
		write(*,'(a,i3,a,i5,a,f6.2,a,f6.2,a,f6.2,a,a,i3,a,i5,a,i5,a)')"[#",myID,", solveHam]: qi=",qi_glob," wann energies= [",En_temp(1)*aUtoEv," : ",&
														En_temp(nWfs)*aUtoEv,"](eV); band gap=",bandGap* aUtoEv," (eV).",&
														" insulating states: #",boundStates,". done tasks=(",qi_loc,"/",qChunk,")"
		!if( boundStates < nWfs) write(*,'(a,i3,a,f8.3,a,f8.3,a)') "[#",myID,", solveHam]: WARNING not enough bound states at qpt=(",qpts(1,qi),",",qpts(2,qi),")."
		!
		get_band_gap = bandGap
		!
		return
	end function







end module ham_Solver
