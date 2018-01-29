module potWellModel
	!this modules sets up the Hamiltonian matrix and solves it for each k point
	!	in the process the wannier functions are generated aswell with routines from wannGen module
	use mpi
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, machineP, myExp, myLeviCivita, &
								eigSolverPART, isUnit, isHermitian
	use sysPara				
	use basisIO,		only:	writeABiN_energy, writeABiN_basis, writeABiN_basCoeff
	implicit none	
	
	private
	public ::					solveHam



		

	


	contains
!public:
	subroutine solveHam()   !call solveHam(wnF, unk, EnW, VeloBwf)
		!solves Hamiltonian at each k point
		!also generates the Wannier functions on the fly (sum over all k)
		!																
		complex(dp),	allocatable		::	Hmat(:,:) , ck_temp(:,:)
		real(dp),		allocatable		::	En_loc(:,:), En_temp(:), En_glob(:,:)
		integer							:: 	qi, qLoc, found, Gsize
		!
		!
		if( myID==root .and. debugHam )		write(*,*)	"[solveHam]: debugging ON. Will do additional tests of the results"
		!
		allocate(	Hmat(				Gmax,	Gmax				)	)
		allocate(	ck_temp(		GmaxGLOBAL, nSolve				)	)
		allocate(	En_temp(				Gmax					)	)	
		allocate(	En_loc(						nSolve	, 	qChunk	)	)
		!
		if( myID == root ) 		allocate(	En_glob(		nSolve, 			nQ		)	)
		if( myID /= root )		allocate(	En_glob(		0,					0		)	)
		!
		qLoc = 1
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
			!COPY TO TARGET ARRAYS
			En_loc(1:nSolve,qLoc)			= En_temp(1:nSolve)
			!
			!DEBUG TESTS
			if( found /= nSolve )	write(*,'(a,i3,a,i5,a,i5)'	)	"[#",myID,";solveHam]: only found ",found," bands of required ",nSolve
			if( nBands > found	)	write(*,'(a,i3,a)'			)	"[#",myID,";solveHam]: warning did not found required amount of bands"
			if( Gsize < nSolve	) 	write(*,'(a,i3,a,i5,a,i5,a)')	"[#",myID,";solveHam]: cutoff to small to get ",nSolve,&
																	" bands! only get",Gsize," basis functions"
			if( Gsize > Gmax	)	write(*,'(a,i3,a)'			)	"[#",myID,";solveHam]: critical error in solveHam, ",&
																	"please contact developer. (nobody but dev will ever read this^^)"
			!
			!WRITE COEFF TO FILE
			call writeABiN_basCoeff(qi, ck_temp)
			!FINALIZE
			write(*,'(a,i3,a,i5,a,i5,a,i5,a)')"[#",myID,", solveHam]: done for qi=",qi," done tasks=(",qLoc,"/",qChunk,")"
			qLoc = qLoc + 1		
		end do
		!		
		!GATHER ON ROOT
		call MPI_GATHER(	En_loc,	nSolve*qChunk,	MPI_DOUBLE_PRECISION,	En_glob,	nSolve*qChunk,	MPI_DOUBLE_PRECISION,	root,	MPI_COMM_WORLD,	ierr)
		!WRITE TO FILE
		if( myID == root ) then
			write(*,*) "root received data, will write to file now"
			call writeABiN_energy(En_glob)
		end if
		!
		!
		return
	end subroutine




!private:
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
			if ( .not.	isHermitian(Hmat)	) then
				write(*,*)"[populateH]: Hamiltonian matrix is not Hermitian :"
			end if
		end if
		!
		!		
		return
	end subroutine




!POTENTIAL WELLS	
	complex(dp)	function V(qLoc, i,j)
		integer,	intent(in)	::	qLoc, i, j
		!
		V	= dcmplx(0.0_dp)
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
		if( norm2(Bext(1:2)) > machineP ) 	write(*,'(a,i3,a)') "[#",myID,";addMagHam]: warning found non zero x or y component. External field should be along z direction! "
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
