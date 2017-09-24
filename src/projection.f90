module projection
	!module contains all the functions to generate the wannier functions
	!	this includes the calculation of the orthonormalized projected bloch wavefunctions
	!	and the FT of the bwfs to calculate the wannier functions
	use omp_lib
	use mathematics, only: dp, PI_dp, acc, myExp, nIntegrate, myMatInvSqrt, isUnit, isIdentity
	use sysPara
	implicit none	
	
	private
	public ::	projectBwf		














	contains
!public:
	subroutine projectBwf(qi, bWf, loBwf, U, failCount, smin, smax)	!projectBwf(qi, bWf, loBwf, U(:,:), failCount, smin, smax)
		!does the projection onto Loewdin-orthonormalized Bloch-like states
		!see Marzari, Vanderbilt PRB 56, 12847 (1997) Sec.IV.G.1 for detailed description of the method 
		integer		,	intent(in)		:: qi
		complex(dp)	,	intent(in)		:: bWf(:,:)   ! bWf(nR,nG)
		complex(dp)	,	intent(out)		:: loBwf(:,:), U(:,:)   !U(nWfs,nWfs)
		integer,		intent(inout)	:: failCount
		real(dp),		intent(inout)	:: smin, smax
		complex(dp)	,	allocatable		:: gnr(:,:), A(:,:), S(:,:), phi(:,:)
		complex(dp)						:: alpha, beta
		integer							:: m, n, k, lda, ldb, ldc, xi
		character*1						:: transa, transb
		!
		allocate(	gnr(nR,nWfs)		)
		allocate( 	S(nWfs,nWfs)		)
		allocate( 	A(nWfs,nWfs)		)
		!
		!SET UP ROTATION MATRIX U
		call genTrialOrb(gnr)
		call calcAmat(bwf,gnr, A)
		write(*,*)	"[projectBwf]: A matrix set up done"
		call calcInvSmat(A, S, smin, smax)
		write(*,*)	"[projectBwf]: S matrix calculation done"
		!U = matmul(S, A)
		transa	= 'n' 
		transb	= 'n'
		m		= size(S,1)
		n		= size(A,2)
		k		= size(S,2)
		alpha	= dcmplx(1.0_dp)
		lda		= size(S,1)
		ldb		= size(A,1)
		beta	= dcmplx(0.0_dp)
		ldc		= size(U,1)
		call zgemm(transa, transb, m, n, k, alpha, S, lda, A, ldb, beta, U, ldc)
		write(*,*)	"[projectBwf]: U matrix set up done"
		!
		!
		!ROTATE BLOCH STATES
		!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(xi)
		do xi = 1, nR
			loBwf(xi,:) = matmul(	U , bWf(xi,1:nWfs)	)
		end do
		!$OMP END PARALLEL DO
		write(*,*)	"[projectBwf]: projections done"
		!
		!
		!DEBUGGING
		if(debugHam) then
			if( .not. isUnit(U) 			) then
				write(*,'(a,i7,a)')"[projectBwf]: qi=",qi," U matrix is not a unitary matrix !"
			end if
			if(	 .not. isOrthonorm(lobWf)	) then
				write(*,'(a,f16.12)')"[projectBwf]: loBwf is NOT a orthonormal basis set, accuracy=",acc
				failCount = failCount + 1
			end if
		end if
		!
		!
		return
	end subroutine
		

	!subroutine genWannF(qi, bWf, wnF)
	!	! generates wannier functions from (projected) bloch wavefunctions
	!	!	DEPRECATED
	!	integer,		intent(in)		:: qi
	!	complex(dp), 	intent(in)  	:: bWf(:,:) ! lobWf(nRpts,nWfs)	
	!	complex(dp), 	intent(inout) 	:: wnF(:,:,:) ! wnF( 	nR, nSC, nWfs		)	
	!	integer 						:: n, Ri, xi
	!	complex(dp)						:: phase
	!	!
	!	!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) DEFAULT(SHARED) PRIVATE(n, Ri, xi, phase) 
	!	do n = 1, nWfs
	!		do Ri = 1, nSC
	!			do xi = 1, nR
	!				phase		 = myExp(	-1.0_dp * dot_product(	qpts(:,qi) , Rcell(:,Ri)	) 	 )
	!				wnF(xi,Ri,n) = wnF(xi,Ri,n) + bWf(xi,n) * phase / real(nQ,dp)
	!			end do
	!		end do
	!	end do
	!	!$OMP END PARALLEL DO
	!	!
	!	!
	!	return
	!end subroutine


	













!privat:
	subroutine genTrialOrb(gnr)
		complex(dp),	intent(out)	:: gnr(:,:)       !gnr( nR, nWfs)	
		real(dp)					:: posX, posY, xc, k, L
		complex(dp)					:: A
		integer						:: n, ri, at
		!
		gnr = dcmplx(0.0_dp)
		do n = 1, nWfs
			do ri = 1, nR
				do at 	= 1, nAt
					if(	insideAt( at , rpts(:,ri) )	) then
						gnr(ri,n) = gVal(at, n, ri)
						!write(*,'(a,i2,a,f6.3,a,f6.3,a,f6.3,a,f6.3,a,f6.3,a,f6.3)')	"[genTrialOrb]: gnr(at=",&
						!			at," r=(",rpts(1,xi),",",rpts(2,xi),&
						!		")= (",posX,",",posY,"))= ",dreal(gnr(xi,n)),"+i*",dimag(gnr(xi,n))
					end if
				end do
				!
			end do
		end do
		!
		return
	end subroutine


	complex(dp) function gVal(at, n, xi)
		!wrapper for calling different trial Orbitals according to input file var trialOrbSw
		integer, 	intent(in)		:: at, n, xi
		!
		select case(trialOrbSw)
			case (1)
				gVal = infPotWell(at, n, xi)
			case default
				gVal = dcmplx(	trialOrbVAL(at) /sum( trialOrbVAL )		)
		end select
		!
		return
	end function


	complex(dp) function infPotWell(at, n, xi)
		integer, 	intent(in)		:: at, n, xi
		real(dp)					:: xc(dim), xrel(dim), k(dim), L(dim)
		complex(dp)					:: A
		!
		!vAvg = abs( sum(vPot) ) * 1.0_dp/ nWells
		L 		= 2.0_dp * atR(:,at)
		xc 		= atPos(:,at)
		k 		= n*PI_dp / L(:)
		A 		= dsqrt(	2.0_dp / product(L)		)
		xrel	= rpts(:,xi) - xc + L / 2.0_dp
		!
		infPotWell = A * dcmplx(		 dsin(	dot_product( k , xrel)		) 	)
		!
		return
	end function



	subroutine calcAmat(bWf,gnr, A)
		!calculates the inner product matrix of bwfs and trial orbital gnr
		complex(dp),	intent(in)	:: bWf(:,:), gnr(:,:) 
		complex(dp),	intent(out)	:: A(:,:) !A(nWfs,nWfs)
		complex(dp),	allocatable	:: f(:)
		integer						:: m,n, xi
		!
	
		A = dcmplx(0.0_dp)
		!
		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, m, f, xi)
		allocate(	f(nR)	)
		!$OMP DO COLLAPSE(2) SCHEDULE(STATIC) 
		do n = 1, nWfs
			do m = 1, nWfs
				f = dcmplx(0.0_dp)
				do xi = 1, nR				
					f(xi) = dconjg(	bWf(xi,m) ) * gnr(xi,n)
				end do
				A(m,n) = nIntegrate(nR, nRx, nRy, dx, dy, f)		
			end do
		end do
		!$OMP END DO
		deallocate(	f )
		!$OMP END PARALLEL
		!
		!
		return
	end subroutine



	subroutine calcInvSmat(A, S, smin, smax)
		!calculates the sqrt inv. of the overlap matrix S
		!
		!	S	= A^dagger A
		complex(dp),	intent(in)		:: A(:,:)
		complex(dp),	intent(out)		:: S(:,:)
		real(dp),		intent(inout)	:: smin, smax
		complex(dp),	allocatable		:: Sold(:,:), Ssqr(:,:)
		integer							:: m,n,k,lda,ldb,ldc, i,j
		real(dp)						:: mi, ma
		complex(dp)						:: alpha, beta
		character*1						:: transa, transb
		!
		allocate(	Sold(nWfs,nWfs)	)
		!CALCULATE S FROM A 
		transa	= 'c'
		transb	= 'n'
		m		= size(A,1)
		n		= size(A,2)
		k		= size(A,2)
		alpha	= dcmplx(1.0_dp)
		beta	= dcmplx(0.0_dp)
		lda		= size(A,1)
		ldb		= size(A,1)
		ldc		= size(A,1)
		!call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
		 call zgemm(transa, transb, m, n, k, alpha, A, lda, A, ldb, beta, S, ldc)
		
		!
		!
		!CALCULATE INVERSE SQRT OF S
		Sold = S
		call myMatInvSqrt(S, mi, ma)

		if( mi < smin ) then
			smin = mi
		end if
		if( ma > smax ) then
			smax = ma
		end if
		!
		!DEBUG
		if(debugHam) then
			write(*,*)"[calcInvSmat]: all done, start debugging test( test if S * (S^-0.5)^2 == I"
			allocate(	Ssqr(nWfs,nWfs)	)
			transa	= 'n'
			transb	= 'n'
			alpha	= dcmplx(1.0_dp)
			beta	= dcmplx(0.0_dp)
			call zgemm(transa, transb, m, n, k, alpha, S   , lda, S   , ldb, beta, Ssqr, ldc)
			call zgemm(transa, transb, m, n, k, alpha, Sold, lda, Ssqr, ldb, beta, Sold, ldc)
			if( .not. isIdentity(Sold) ) then
				write(*,*)"[calcInvSmat]: problem with mat inversion, seems to be not inverse square root"
			end if
		end if
		!
		!
		return
	end subroutine


	logical function isOrthonorm(loBwf)
		!cheks wether loBwf is orthonormal by calculating the overlap matrix elements
		!	<psi_n,k|psi_n,k'> = nSC * \delta(n,m)
		!	WARNING: DOES NOT CHECK ORTHONORMALIZATION BETWEEN K POINTS
		complex(dp),	intent(in)		:: lobWf(:,:)   !lobWf(nR,nWfs)
		complex(dp),	allocatable		:: f(:)
		complex(dp)						:: oLap
		
		integer							:: n, m, ri
		logical							:: realist, imaginist
		allocate( f(nR)		)
		!
		isOrthonorm = .true.
		n 			= 1
		!
		do  while(n <= nWfs .and. isOrthonorm) 
			!DIAGONAL ELEMENTS
			!
			!
			!PREPARE INTEGRATION ARRAY
			do ri = 1, nR
				f(ri) = dconjg( lobWf(ri,n) ) * lobWf(ri,n)
			end do
			oLap	= nIntegrate(nR, nRx, nRy, dx, dy, f)
			oLap	= oLap / dcmplx( nSC	)
			!CHECK CONDITIONS
			realist 	= abs(			abs(dreal(oLap))		-		1.0_dp		) 	> acc
			imaginist	= abs(			dimag(oLap)								)	> acc
			if( realist .or. imaginist ) then
				write(*,'(a,i2,a,f15.10,a,f15.10)')"[isOrthonorm]: overlap(n,n=",n,") = ",dreal(oLap),"+i*",dimag(oLap)
				isOrthonorm = .false.
			end if
			!
			!
			!OFF DIAGONAL
			m = 1
			do while(m <= nWfs	.and. isOrthonorm)
				if(n /= m) then
					!PREPARE INTEGRATION ARRAY
					do ri = 1, nR
						f(ri) = dconjg( lobWf(ri,n) ) * lobWf(ri,m)
					end do
					oLap 	= nIntegrate(nR, nRx, nRy, dx, dy, f)
					oLap	= oLap / dcmplx( nSC	)
					!CHECK CONDITION
					if(abs(oLap) > acc) then
						write(*,'(a,i2,a,i2,a,f15.10,a,f15.10)')"[isOrthonorm]: overlap(",n,",",m,") = ",dreal(oLap),"+i*",dimag(oLap)
						isOrthonorm = .false.
					end if
				end if
				m = m +1
			end do
			n = n + 1
		end do
		!
		return
	end function
		


end module projection