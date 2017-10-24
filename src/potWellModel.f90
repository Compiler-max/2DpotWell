module potWellModel
	!this modules sets up the Hamiltonian matrix and solves it for each k point
	!	in the process the wannier functions are generated aswell with routines from wannGen module
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, machineP, myExp, myLeviCivita, eigSolver, nIntegrate, isUnit, isHermitian
	use sysPara
	use blochWf,		only:	genBwfVelo, testNormUNK							
	use output,			only:	printMat, 	writeEnAndCK
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
		complex(dp),	allocatable		::	Hmat(:,:), unkT(:,:)
		real(dp),		allocatable		::	EnT(:,:)
		integer							:: 	qi, gi, n, nCut
		real(dp)						::	kVal(2)
		!
		allocate(	Hmat(	nG,	nG		)			)
		allocate(	EnT(	nG	, nQ	)			)	
		!
		!
		if(debugHam) then
			write(*,*)	"[solveHam]: debugging ON. Will do additional tests of the results"
		end if
		!
		!
		!MPI scatter : unk, EnT
		do qi = 1, nQ
			kVal	=	qpts(:,qi)
			!
			!ELECTRONIC STRUCTURE

			call populateH(kVal, Hmat) 	!omp
			call eigSolver(Hmat, EnT(:,qi))	!mkl
			ck(1:nG,1:nBands,qi)	= Hmat(1:nG,1:nBands)

			!
			!BLOCH WAVEFUNCTIONS
			!call gaugeCoeff(kVal, Hmat)
			!call genBwfVelo(qi, Hmat, unkT)	!omp
			!unk(:,1:nBands,qi)	= unkT(:,1:nBands)
			!
			!DEBUG TESTS
			if( debugHam ) then
				if( .not. isUnit(Hmat) ) write(*,'(a,i4)')	"[solveHam]: basis coeff not unitary for qi=",qi
			end if
			write(*,*)"[solveHam]: done for qi=",qi
			!
		end do
		!
		!
		!COPY & WRITE ENERGIES/BWFs
		do qi = 1, size(En,2)
			En(:,qi)	= EnT(1:nBands,qi)
		end do
		write(*,*)			"[solveHam]: copied eigenvalues"
		write(*,*)			"[solveHam]: found ", countBandsSubZero(EnT)," bands at the gamma point beneath zero"
		call writeEnAndCK(EnT, ck)
		!
		!
		!
		return
	end subroutine














!private:
	!POPULATION OF H MATRIX
	subroutine populateH(q, Hmat)
		!populates the Hamiltonian matrix by adding 
		!	1. qinetic energy terms
		!	2. potential terms
		!and checks if the resulting matrix is still hermitian( for debugging)
		real(dp)   , intent(in)    	:: q(2)
		complex(dp), intent(inout) 	:: Hmat(:,:)
		real(dp)					:: kgi(2), kgj(2)
		complex(dp)					:: onSite
		integer						:: i, j
		!init to zero
		Hmat = dcmplx(0.0_dp) 
		!
		!GET CUTOFF



		!!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j, i, kgi, kgj, onSite)
		do j = 1, nG
			do i = 1, nG
				kgi(:) 	= q(:) + Gvec(:,i)
				kgj(:)	= q(:) + Gvec(:,j)
				if( norm2(kgi) < Gcut .and. norm2(kgj) < Gcut ) then
					if(i == j )	then	
						onSite	= 0.5_dp * 	dot_product(kgi,kgi)
						Hmat(i,j)	=	V(i,j)	+	onSite
					else
						Hmat(i,j)	=	V(i,j)
					end if
				end if
			end do
		end do
		!!!$OMP END PARALLEL DO

		!
		if(debugHam) then
			if ( .not.	isHermitian(Hmat)	) then
				write(*,*)"[populateH]: Hamiltonian matrix is not Hermitian :"
				!call printMat(nG, Hmat)
			end if
		end if


		return
	end subroutine

	complex(dp)	function V(i,j)
		integer,	intent(in)	::	i, j
		!
		V	= dcmplx(0.0_dp)
		if( doVdesc ) then
			V = Vdesc(i,j)
		else
			V = Vconst(i,j)
		end if
		!
		!
		return
	end function


	complex(dp) function Vconst(i,j)
		!calc potential matrix elements
		!the integrals were solved analytical and are hard coded in this function
		integer,	intent(in)	::	i, j
		integer					::	at
		complex(dp)				::	Vpot
		real(dp)				::  xL, yL, xR, yR, dGx, dGy
		!
		Vconst 	= dcmplx(0.0_dp)
		dGx		= Gvec(1,j) - Gvec(1,i)
		dGy		= Gvec(2,j) - Gvec(2,i)
		!
		do at = 1, nAt
			Vpot	=	dcmplx(atPot(at))
			xL	=	atPos(1,at) - atR(1,at)
			xR	=	atPos(1,at) + atR(1,at) 
			yL	=	atPos(2,at) - atR(2,at)
			yR	=	atPos(2,at) + atR(2,at) 
			!
			if( i == j) then		
				Vconst	= Vconst + Vpot 			*	 ( xR - xL ) * 	( yR - yL )			 / vol
			else if( abs(dGx) < machineP ) then	
				Vconst	= Vconst + Vpot  * i_dp 	* 	( xR - xL ) *( myExp(dGy*yL) - myExp(dGy*yR) )/ ( vol * dGy )
			else if( abs(dGy) < machineP ) then
				Vconst	= Vconst + Vpot * i_dp	 * ( myExp(dGx*xL) - myExp(dGx*xR) ) 	* 	( yR - yL) 		/ (vol * dGx )
			else
				Vconst	= Vconst -  Vpot 	 * ( myExp(dGx*xL) - myExp(dGx*xR) ) * ( myExp(dGy*yL) - myExp(dGy*yR) ) / (vol * dGx * dGy )
			end if
		end do
		!
		return
	end function


	complex(dp) function Vdesc(i,j)
		!potential integration for a well linear descending in x direction
		!it starts from V0 at xL till V0-dV at xR
		integer,	intent(in)	::	i, j
		integer					::	at
		complex(dp)				::	Vpot
		real(dp)				::  xL, yL, xR, yR, dGx, dGy, dV, fact
		!
		Vdesc 	= dcmplx(0.0_dp)
		dGx		= Gvec(1,j) - Gvec(1,i)
		dGy		= Gvec(2,j) - Gvec(2,i)
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
		do n = 1, nG
			if( EnT(n,gammaP) < 0.0_dp )	then
				cnt	= cnt + 1
			end if
		end do
		countBandsSubZero	= cnt
		!
		return
	end function



end module potWellModel


















	!GAUGING BASIS COEFFICIENTS
	!subroutine gaugeCoeff(k, basCoeff)
	!	!method gauges the basis coefficients directly
	!	!	controled via the gaugeSwitch from the input file
	!	!	this should only be used for testing & understanding purposes,
	!	!	since this gauge trafo is not saved in the U matrix (U matrix rotates between ham & wann gauge)
	!	real(dp),		intent(in)		:: k(2)
	!	complex(dp),	intent(inout) 	:: basCoeff(:,:)
	!	!
	!	select case (gaugeSwitch)
	!		case (1)
	!			write(*,*) "[gaugeCoeff]: gauging with G=0 phase"
	!			call gaugeCoeff1(basCoeff)
	!		!case (2)
	!		!	write(*,*) "[gaugeCoeff]: gauging with G(nG0-1) for k<=0, and G(nG0+1) for k>0. (G(nG0)=0)"
	!		!	call gaugeCoeff2(k,basCoeff)
	!		case (0)
	!			!write(*,*) "[gaugeCoeff]: use gauge convention from solver"
	!			! do nothing
	!		case default
	!			write(*,*) "[gaugeCoeff]: unknown gaugeSwitch, eigVecs wont be gauged"
	!	end select
	!	!
	!	return
	!end subroutine


	!subroutine gaugeCoeff1(basCoeff)
	!	!Gauges basCoeff with the G=0 component, effectively removes k-dependence of basCoeff gauge
	!	complex(dp)	, intent(inout) :: basCoeff(:,:)
	!	integer						:: i1,i2
	!	complex(dp)				  	:: phase0
	!	real(dp)				   	:: delta
	!	!
	!	delta=1E-14_dp
	!	!
	!	do i2=1,nG !loop eigVecs
	!		if(	abs( basCoeff(nG0,i2) )	> delta	) then
	!			phase0 = basCoeff(nG0,i2) / abs(basCoeff(nG0,i2)) 
	!			do i1=1,nG !gauge each vector
	!					basCoeff(i1,i2) = basCoeff(i1,i2)/ phase0
	!			end do
	!		endif
	!	end do
	!	!
	!	return
	!end subroutine

	!subroutine gaugeCoeff2(k, basCoeff)
	!	!what does k < 0 mean here ?!, makes only sense in 1D
	!	!for negative k gauge with Gmin, for positive with Gmax
	!	real(dp),		intent(in)		:: k(dim)
	!	complex(dp), 	intent(inout) 	:: basCoeff(:,:)
	!	integer 				   		:: i1,i2
	!	complex(dp)				   		:: phase0
	!	!
	!	do i2=1,nG
	!		if( k <= 0.0_dp) then
	!			phase0 = basCoeff(nG0-1 ,i2) / abs(basCoeff(nG0-1 ,i2))
	!		else
	!			phase0 = basCoeff(nG0+1,i2) / abs(basCoeff(nG0+1,i2))
	!		end if
	!		do i1=1,nG
	!			basCoeff(i1,i2) = basCoeff(i1,i2)/ phase0
	!		end do
	!	end do
	!	!
	!	return
	!end

