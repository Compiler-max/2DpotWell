module potWellModel
	!this modules sets up the Hamiltonian matrix and solves it for each k point
	!	in the process the wannier functions are generated aswell with routines from wannGen module
	use mpi
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, machineP, myExp, myLeviCivita, &
								eigSolverPART, isUnit, isHermitian
	use sysPara				
	use planeWave,		only:	calcVeloGrad, calcAmatANA, calcMmat
	use basisIO,		only:	writeABiN_basVect, writeABiN_energy, writeABiN_basCoeff, writeABiN_velo, writeABiN_Amn, writeABiN_Mmn, &
								read_coeff, read_gVec
	use w90Interface,	only:	setup_w90, write_w90_matrices
	implicit none	
	!#include "mpif.h"

	private
	public ::			potWellMethod			



		

	


	contains
!public:
	subroutine potWellMethod()   !call solveHam(wnF, unk, EnW, VeloBwf)
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
		call solveHam()
		!
		!get FD scheme from w90
		if( myID == root )  then
			write(*,'(a,i3,a)')	"[#",myID,";potWellMethod]: Hamiltonian solved, setup w90 now"
			call setup_w90(nntot, nnlist, nncell)
			write(*,'(a,i3,a)')	"[#",myID,";potWellMethod]: wannier setup done, now calc M matrix"
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



		!call w90 interface to write input files & .mmn etc. files
		if( myID == root ) then
			call write_w90_matrices()
			write(*,'(a,i3,a)')		"[#",myID,";solveHam]: wrote w90 matrix input files (.win, .amn, .mmn, .eig, _geninterp.kpt )"
		end if


		!
		return
	end subroutine




!private:
	subroutine solveHam()
		!			solve Ham, write results and derived quantites														
		complex(dp),	allocatable		::	Hmat(:,:) , ck_temp(:,:), Amn_temp(:,:), velo_temp(:,:,:)
		real(dp),		allocatable		::	En_temp(:)
		integer							:: 	qi, qLoc, found, Gsize
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
			write(*,'(a,i3,a)'	)	"[#",myID,";solveHam]: I solved the Ham"
			!
			!get derived quantities
			call calcAmatANA(qLoc, ck_temp, Amn_temp)
			write(*,'(a,i3,a)'	)	"[#",myID,";solveHam]: I calculated A mat"
			call calcVeloGrad(qLoc, ck_temp, velo_temp)
			write(*,'(a,i3,a)'	)	"[#",myID,";solveHam]: I calculated velocities"

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
			write(*,'(a,i3,a,i5,a,f6.2,a,i5,a,i5,a)')"[#",myID,", solveHam]: done for qi=",qi," lowest energy=",En_temp(1),"[Hatree] done tasks=(",qLoc,"/",qChunk,")"
			qLoc = qLoc + 1		
		end do
		!
		!
		return
	end subroutine


	subroutine calc_Mmat(nntot, nnlist, nncell)
		integer,		intent(in)		::	nntot, nnlist(:,:), nncell(:,:,:)
		integer,		allocatable		::	nGq_glob(:), nG_qi, nG_nn, q_nn
		complex(dp),	allocatable		::	ck_qi(:,:), cK_nn(:,:), Mmn(:,:,:)
		real(dp),		allocatable		::	Gvec_qi(:,:), Gvec_nn(:,:)
		real(dp)						::	gShift(2)
		integer							::	qi, nn
		!
		allocate(	nGq_glob(nQ)				)
		allocate(	Gvec_qi(dim, nG)			)
		allocate(	Gvec_nn(dim, nG)			)
		allocate(	ck_nn(GmaxGLOBAL, nSolve)	)
		allocate(	ck_qi(GmaxGLOBAL, nSolve)	)
		allocate(	Mmn(nBands, nBands, nntot)	)
		!
		call MPI_GATHER( nGq	, qChunk, MPI_INTEGER, nGq_glob		, qChunk, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
		if( myID == root ) write(*,*)	"[calc_Mmat]: gathered nGq"
		!
		if(myID == root) then
			write(*,*)	"qChunk=",qChunk
			write(*,*)	"fd sheme:"
			write(*,*)	"nntot=",nntot
			write(*,*)	nnlist
			write(*,*)	nncell
		end if
		!
		do qi = myID*qChunk +1, myID*qChunk + qChunk
			!
			do nn = 1, nntot
				!calc overlap of unks
				if( nncell(3,qi,nn)/= 0 ) stop '[w90prepMmat]: out of plane nearest neighbour found. '
				
				q_nn		=	nnlist(qi,nn)
				nG_qi		= 	nGq_glob(qi)
				nG_nn		= 	nGq_glob(q_nn)
				gShift(1)	= 	nncell(1,qi,nn) * 2.0_dp * PI_dp / aX
				gShift(2)	= 	nncell(2,qi,nn) * 2.0_dp * PI_dp / aY
				!
				!read basis coefficients
				call read_coeff(qi,	ck_qi)
				call read_coeff(q_nn, ck_nn)
				write(*,'(a,i3,a)')		"[#",myID,",calc_Mmat]: read coff"
				!
				!read Gvec
				call read_gVec(qi, 		Gvec_qi)
				call read_gVec(q_nn,	Gvec_nn)
				write(*,'(a,i3,a)')		"[#",myID,",calc_Mmat]: read gvec"
				!
				!
				call calcMmat(qi, q_nn, gShift, nG_qi, nG_nn, Gvec_qi, Gvec_nn, ck_qi, ck_nn, Mmn(:,:,nn)	)
				write(*,'(a,i3,a)')		"[#",myID,",calc_Mmat]: calulated Mmat"
			end do
			!
			!write result to file
			call writeABiN_Mmn(qi, Mmn)
			write(*,'(a,i3,a,i5,a,i3)')		"[#",myID,",calc_Mmat]: done setting up M_matrix for qi=",qi," nntot=",nntot
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
		integer,		intent(in)	:: qLoc
		complex(dp), intent(inout) 	:: Hmat(:,:)
		complex(dp)					:: onSite
		integer						:: i, j
		!init to zero
		Hmat = dcmplx(0.0_dp) 
		!
		do j = 1, nGq(qLoc)
			do i = 1, nGq(qLoc)
				!KINETIC ENERGY + POTENTIAL WELLS
				if(i == j )	then	
					onSite	= 0.5_dp * 	dot_product(Gvec(:,i,qLoc),Gvec(:,i,qLoc))
					Hmat(i,j)	=	V(qLoc, i,j)	+	onSite
				else
					Hmat(i,j)	=	V(qLoc, i,j)
				end if
			end do
		end do
		!
		!ADD PEIERLS
		if( doMagHam )	 call addMagHam( qLoc, Hmat)
		!
		!DEBUG
		if(debugHam) then
			if ( .not.	isHermitian(Hmat)	) 	write(*,'(a,i3,a,i3)')	"[#",myID,";populateH]: WARNING Hamiltonian matrix is not Hermitian at qLoc=",qLoc
		end if
		!
		!		
		return
	end subroutine




!POTENTIAL WELLS	
	complex(dp)	function V(qLoc, i,j)
		integer,	intent(in)	::	qLoc, i, j
		!
		if( doVdesc ) then
			V = Vdesc(qLoc, i,j)
		else	
			V = Vconst(qLoc, i,j)
		end if
		!
		!
		return
	end function


	complex(dp) function Vconst(qLoc, i,j)
		!calc potential matrix elements
		!the integrals were solved analytical and are hard coded in this function
		integer,	intent(in)	::	qLoc,  i, j
		integer					::	at
		complex(dp)				::	Vpot
		real(dp)				::  xL, yL, xR, yR, dGx, dGy
		!
		Vconst 	= dcmplx(0.0_dp)
		dGx		= Gvec(1,j,qLoc) - Gvec(1,i,qLoc) 
		dGy		= Gvec(2,j,qLoc) - Gvec(2,i,qLoc) 
		!
		do at = 1, nAt
			Vpot	=	dcmplx(atPot(at))
			xL	=	atPos(1,at) - atR(1,at)
			xR	=	atPos(1,at) + atR(1,at) 
			yL	=	atPos(2,at) - atR(2,at)
			yR	=	atPos(2,at) + atR(2,at) 
			!
			if( i == j) then		
				Vconst	= Vconst + Vpot 			*	( xR - xL ) * 	( yR - yL )			 			/ vol
			else if( abs(dGx) < machineP ) then	
				Vconst	= Vconst + Vpot  * i_dp 	* 	( xR - xL ) * ( myExp(dGy*yL) - myExp(dGy*yR) )	/ ( vol * dGy )
			else if( abs(dGy) < machineP ) then
				Vconst	= Vconst + Vpot * i_dp	 	* 	( yR - yL) 	* ( myExp(dGx*xL) - myExp(dGx*xR) ) / (vol * dGx )
			else
				Vconst	= Vconst -  Vpot 	 * ( myExp(dGx*xL) - myExp(dGx*xR) ) * ( myExp(dGy*yL) - myExp(dGy*yR) ) / (vol * dGx * dGy )
			end if
		end do
		!
		return
	end function


	complex(dp) function Vdesc(qLoc, i,j)
		!potential integration for a well linear descending in x direction
		!it starts from V0 at xL till V0-dV at xR
		integer,	intent(in)	::	qLoc, i, j
		integer					::	at
		complex(dp)				::	Vpot
		real(dp)				::  xL, yL, xR, yR, dGx, dGy, dV, fact
		!
		Vdesc 	= dcmplx(0.0_dp)
		dGx		= Gvec(1,j,qLoc) - Gvec(1,i,qLoc) 
		dGy		= Gvec(2,j,qLoc) - Gvec(2,i,qLoc) 
		!
		do at = 1, nAt
			Vpot=	dcmplx(atPot(at))
			dV	=	dVpot(at)
			xL	=	atPos(1,at) - atR(1,at)
			xR	=	atPos(1,at) + atR(1,at) 
			yL	=	atPos(2,at) - atR(2,at)
			yR	=	atPos(2,at) + atR(2,at) 
			!
			if( i == j) then
				fact	= (2.0_dp*Vpot - dV)	 / (2.0_dP*vol)	
				Vdesc	= Vdesc 	+			fact *	(xL-xR) * 	(yL-yR)			
				!
				!
				!
			else if( abs(dGx) < machineP ) then	
				fact	= (2.0_dp*Vpot - dV) / ( 2.0_dP* vol * dGy )
				Vdesc	= Vdesc 	- 	i_dp  * fact * (xL-xR) * ( myExp(dGy*yL) - myExp(dGy*yR) )
				!
				!
				!
			else if( abs(dGy) < machineP ) then
				fact	=  1.0_dp / ( dGx**2 * vol * (xL-xR) )
				Vdesc	= Vdesc 	+ 	i_dp * fact * (yR-yL) * myExp(dGx*xL) * (dGx * Vpot			* (xL-xR) + i_dp * dV) 	
				Vdesc	= Vdesc 	-	i_dp * fact * (yR-yL) * myExp(dGx*xR) * (dGx * (Vpot-dV) * (xL-xR) + i_dp * dV)	
				!
				!
				!
			else
				fact	= -1.0_dp * ( myExp(dGy*yL)-myExp(dGy*yR) ) 	/ ( dGx**2 * dGy * vol * (xL-xR) )
				Vdesc	= Vdesc		+			fact * 			myExp(dGx*xL) * (dGx * Vpot			* (xL-xR) + i_dp * dV) 
				Vdesc	= Vdesc		-			fact *			myExp(dGx*xR) * (dGx * (Vpot-dV)	* (xL-xR) + i_dp * dV)
				!
				!
				!
			end if
		end do
		!
		!
		return
	end function








!EXTERNAL MAGNETIC FIELD ( OSCILLATING )
	subroutine addMagHam( qLoc, Hmat)
		!Performs Peierls substitution on plane wave basis.
		!
		!neglects terms which are second order in the external field
		!
		!only contribtutions if i/=j, i.e. to off-diagonal terms
		! no contribution for dGx = 0
		integer,		intent(in)		::	qLoc
		complex(dp),	intent(inout)	::	Hmat(:,:)
		integer							::	i, j
		real(dp)						::	H_prefact, qX_Period, dGx, dGy
		!
		!period of oscillating B field
		qX_period	=	2.0_dp * PI_dp / aX 
		H_prefact 	= 	0.5_dp * Bext(3) / qX_period
		!
		!Add to Hamiltonian
		do j = 1, nGq(qLoc)
			do i = 1, nGq(qLoc)
				if( i/=j ) then
						dGx		= Gvec(1,j,qLoc) - Gvec(1,i,qLoc) 
						dGy		= Gvec(2,j,qLoc) - Gvec(2,i,qLoc) 
						if( abs(dGx) > machineP  ) then
							if( abs(dGy) < machineP ) then
								Hmat(i,j)	= Hmat(i,j) + H_prefact * h1_gyZero( dGx, qX_period)
							else 
								Hmat(i,j)	= Hmat(i,j) + H_prefact * h1_gyFull( dGx, dGy, qX_period)
							end if 
						end if
				end if
			end do
		end do
		!
		!DEBUG
		if( myID == root .and. qLoc == 1) 	write(*,'(a,i3,a)')	"[#",myID,";addMagHam]: hello there"		
		if( norm2(Bext(1:2)) > machineP ) 	write(*,'(a,i3,a)') "[#",myID,";addMagHam]: WARNING found non zero x or y component. External field should be along z direction! "
		!
		return
	end subroutine



	complex(dp)	function h1_gyZero( dGx, qX )
		!matrix element for dGy==0
		real(dp),		intent(in)		:: dGx, qX
		!
		h1_gyZero	= 	qX 		*	( myExp(dGx*aX) - 1.0_dp )		*	aY
		!
		return
	end function

	complex(dp)	function h1_gyFull( dGx, dGy, qX)
		!matrix element for dGy/=0
		real(dp),		intent(in)		:: dGx, dGy, qX
		!
		h1_gyFull	=	qX		*	( myExp(dGx*aX) - 1.0_dp )		* 	( myExp(dGy*aY) - 1.0_dp )
		!
		return
	end function


end module potWellModel
