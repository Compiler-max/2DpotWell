module potWellModel
	!this modules sets up the Hamiltonian matrix and solves it for each k point
	!	in the process the wannier functions are generated aswell with routines from wannGen module
	use mpi
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, machineP, myExp, myLeviCivita, &
								eigSolverPART, isUnit, isHermitian
	use sysPara				
	implicit none	
	
	private
	public ::					solveHam



		

	


	contains
!public:
	subroutine solveHam(ck, En)   !call solveHam(wnF, unk, EnW, VeloBwf)
		!solves Hamiltonian at each k point
		!also generates the Wannier functions on the fly (sum over all k)
		complex(dp),	intent(out)		::	ck(:,:,:)
		real(dp),		intent(out)		::	En(:,:)	
		!																
		complex(dp),	allocatable		::	Hmat(:,:) , ctemp(:,:)
		real(dp),		allocatable		::	EnT(:)
		integer							:: 	qi, qLoc, found, Gsize
		!									qi: qpts grid index. qLoc: lcoal storage adress for alll quanties at qpts=qi
		!
		!
		if( myID==root .and. debugHam )		write(*,*)	"[solveHam]: debugging ON. Will do additional tests of the results"
			

		allocate(	Hmat(	Gmax,	Gmax	)			)
		allocate(	ctemp(	Gmax, nSolve	)			)
		allocate(	EnT(	Gmax			)			)	

		qLoc = 1
		do qi = myID*qChunk +1, myID*qChunk + qChunk
			!
			!!SETUP HAM
			call populateH(qLoc, Hmat) 
			!!
			!SOLVE HAM
			ctemp	= dcmplx(0.0_dp)
			EnT		= 0.0_dp
			Gsize 	= nGq(qLoc)
			if( Gsize > Gmax	)	write(*,'(a,i3,a)'			)	"[#",myID,";solveHam]: critical error in solveHam, ",&
																		"please contact developer. (nobody but dev will ever read this^^)"
			call eigSolverPART(Hmat(1:Gsize,1:Gsize),EnT(1:Gsize), ctemp(1:Gsize,:), found)
			!
			!
			!COPY TO TARGET ARRAYS
			ck(1:Gmax,1:nSolve,qLoc)	= ctemp(1:Gmax,1:nSolve)
			En(1:nSolve,qLoc)			= EnT(1:nSolve)

			!DEBUG TESTS
			if( found /= nSolve )	write(*,'(a,i3,a,i5,a,i5)'	)	"[#",myID,";solveHam]: only found ",found," bands of required ",nSolve
			if( nBands > found	)	write(*,'(a,i3,a)'			)	"[#",myID,";solveHam]: warning did not found required amount of bands"
			if( Gsize < nSolve	) 	write(*,'(a,i3,a,i5,a,i5,a)')	"[#",myID,";solveHam]: cutoff to small to get ",nSolve,&
																	" bands! only get",Gsize," basis functions"
			!
			!FINALIZE
			qLoc = qLoc + 1
			write(*,'(a,i3,a,i5)')"[#",myID,", solveHam]: done for qi=",qi
		end do
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
					if(i == j )	then	
						onSite	= 0.5_dp * 	dot_product(Gvec(:,i,qLoc),Gvec(:,i,qLoc))
						Hmat(i,j)	=	V(qLoc, i,j)	+	onSite
					else
						Hmat(i,j)	=	V(qLoc, i,j)
					end if
				!end if
			end do
		end do
		!
		if(debugHam) then
			if ( .not.	isHermitian(Hmat)	) then
				write(*,*)"[populateH]: Hamiltonian matrix is not Hermitian :"
			end if
		end if

		return
	end subroutine

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


	integer function countBandsSubZero(EnT)
		real(dp),		intent(in)		:: EnT(:,:)
		integer							:: gammaP, cnt, n
		!
		cnt	= 0
		gammaP	= getGammaPoint()
		do n = 1, size(EnT,1)
			if( EnT(n,gammaP) < 0.0_dp )	then
				cnt	= cnt + 1
			end if
		end do
		countBandsSubZero	= cnt
		!
		return
	end function



end module potWellModel
