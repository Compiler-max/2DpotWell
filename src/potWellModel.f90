module potWellModel
	!this modules sets up the Hamiltonian matrix and solves it for each k point
	!	in the process the wannier functions are generated aswell with routines from wannGen module
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, machineP, myExp, myLeviCivita, &
								eigSolver, eigSolver2, nIntegrate, isUnit, isHermitian
	use sysPara				
	use output,			only:	writeEnAndCK
	implicit none	
	
	private
	public ::					solveHam



		

	


	contains
!public:
	subroutine solveHam(ck, En)   !call solveHam(wnF, unk, EnW, VeloBwf)
		!solves Hamiltonian at each k point
		!also generates the Wannier functions on the fly (sum over all k)
		complex(dp),	intent(out)		::	ck(:,:,:)
																				!wnF( nR	, nSupC,	nWfs	)	
																				!unkW(nR	,	nWfs,  nKpts	)
																				!veloBwf(nR,nK,2*nG)
		real(dp),		intent(out)		::	En(:,:)																	
		complex(dp),	allocatable		::	Hmat(:,:) , ctemp(:,:)
		real(dp),		allocatable		::	EnT(:,:)
		integer							:: 	qi, found, Gmax
		!
		allocate(	Hmat(	nG,	nG		)			)
		allocate(	ctemp(	nG,	nSolve	)			) !could be smaller, just for testing resized
		allocate(	EnT(	nG, nQ	)			)	
		!
		!
		if(debugHam) then
			write(*,*)	"[solveHam]: debugging ON. Will do additional tests of the results"
		end if
		!
		!
		!MPI scatter : ck(:,:,qi), EnT(:,qi)
		do qi = 1, nQ
			!
			!ELECTRONIC STRUCTURE
			call setUpBasis(qi)
			Gmax 	= nGq(qi)
			

			call populateH(qi, Hmat) 	!omp
			
			!old
			!call eigSolver(Hmat, EnT(:,qi))	!mkl
			!ck(1:nG,1:nBands,qi)	= Hmat(1:nG,1:nBands)

			!new
			call eigSolver2(Hmat(1:Gmax,1:Gmax),EnT(1:Gmax,qi), ctemp(1:Gmax,:), found)!a, w ,z, m
			
			ck(1:nG,1:nBands,qi)	= ctemp(1:nG,1:nBands)
			!DEBUG TESTS
			if( debugHam ) then
				if(found /= nSolve )write(*,*)"[solveHam]: only found ",found," bands of required ",nSolve
				if(nBands > found) write(*,*)"[solveHam]: warning did not found required amount of bands"
				if( Gmax < nSolve ) write(*,*)"[solveHam]: cutoff to small to get ",nSolve," bands! only get",Gmax," basis functions"
			end if
			!
			!
			write(*,*)"[solveHam]: done for qi=",qi
			!
		end do
		!MPI gather, ck(:,:,qi), EnT(:,qi)
		!
		!
		!COPY & WRITE ENERGIES/BWFs
		do qi = 1, size(En,2)
			En(:,qi)	= EnT(1:nBands,qi)
		end do
		write(*,*)			"[solveHam]: copied eigenvalues"
		write(*,*)			"[solveHam]: found ", countBandsSubZero(EnT(1:nSolve,:))," bands at the gamma point beneath zero"
		call writeEnAndCK(EnT, ck)
		!
		!
		!
		return
	end subroutine














!private:
	!BASIS SETUP
	subroutine setUpBasis(qi)
		integer,	intent(in)		:: qi
		integer						:: gi
		real(dp)					:: kg(2)
		!
		nGq(qi)	= 0
		do gi = 1, nG
			kg(:)	= qpts(:,qi) + Gtest(:,gi)
			if( norm2(kg) < Gcut ) then
				nGq(qi) = nGq(qi) + 1
				Gvec(:,nGq(qi),qi) = kg(:)
			end if
		end do
		!
		if(nGq(qi) > nG) write(*,'(a,i4,a,i6)')	"[setUpBasis]: warning, somehow counted more basis functions at qi=",qi," limit nG=",nG	
		write(*,'(a,i6,a,i4)')	"[setUpBasis]: using ",nGq(qi), "basis functions at qi=",qi
		!
		!
		return
	end subroutine




	!POPULATION OF H MATRIX
	subroutine populateH(qi, Hmat)
		!populates the Hamiltonian matrix by adding 
		!	1. qinetic energy terms
		!	2. potential terms
		!and checks if the resulting matrix is still hermitian( for debugging)
		integer,		intent(in)	:: qi
		complex(dp), intent(inout) 	:: Hmat(:,:)
		complex(dp)					:: onSite
		integer						:: i, j
		!init to zero
		Hmat = dcmplx(0.0_dp) 
		!
		!GET CUTOFF

	

		!!!$OMP PARALLEL DO SCHEDULE(DYNAMIC,nG/8) DEFAULT(SHARED) FIRSTPRIVATE(q) PRIVATE(j, i, kgi, kgj, onSite)
		do j = 1, nGq(qi)
			do i = 1, nGq(qi)
					if(i == j )	then	
						onSite	= 0.5_dp * 	dot_product(Gvec(:,i,qi),Gvec(:,i,qi))
						Hmat(i,j)	=	V(qi, i,j)	+	onSite
					else
						Hmat(i,j)	=	V(qi, i,j)
					end if
				!end if
			end do
		end do
		!!!$OMP END PARALLEL DO

		!
		if(debugHam) then
			if ( .not.	isHermitian(Hmat)	) then
				write(*,*)"[populateH]: Hamiltonian matrix is not Hermitian :"
			end if
		end if

		return
	end subroutine





	complex(dp)	function V(qi, i,j)
		integer,	intent(in)	::	qi, i, j
		!
		V	= dcmplx(0.0_dp)
		if( doVdesc ) then
			V = Vdesc(qi, i,j)
		else
			V = Vconst(qi, i,j)
		end if
		!
		!
		return
	end function


	complex(dp) function Vconst(qi, i,j)
		!calc potential matrix elements
		!the integrals were solved analytical and are hard coded in this function
		integer,	intent(in)	::	qi, i, j
		integer					::	at
		complex(dp)				::	Vpot
		real(dp)				::  xL, yL, xR, yR, dGx, dGy
		!
		Vconst 	= dcmplx(0.0_dp)
		dGx		= Gvec(1,j,qi) - Gvec(1,i,qi) 
		dGy		= Gvec(2,j,qi) - Gvec(2,i,qi) 
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


	complex(dp) function Vdesc(qi, i,j)
		!potential integration for a well linear descending in x direction
		!it starts from V0 at xL till V0-dV at xR
		integer,	intent(in)	::	qi, i, j
		integer					::	at
		complex(dp)				::	Vpot
		real(dp)				::  xL, yL, xR, yR, dGx, dGy, dV, fact
		!
		Vdesc 	= dcmplx(0.0_dp)
		dGx		= Gvec(1,j,qi) - Gvec(1,i,qi) 
		dGy		= Gvec(2,j,qi) - Gvec(2,i,qi) 
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
