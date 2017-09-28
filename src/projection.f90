module projection
	!module contains all the functions to generate the wannier functions
	!	this includes the calculation of the orthonormalized projected bloch wavefunctions
	!	and the FT of the bwfs to calculate the wannier functions
	use omp_lib
	use mathematics, 	only: 	dp, PI_dp, acc, myExp, nIntegrate,  mySVD, myMatInvSqrt, isUnit, isIdentity
	use blochWf,		only:	genUnk, testNormUNK
	use sysPara
	implicit none	
	
	private
	public ::	projectUnk		














	contains
!public:

	subroutine projectUnk(En, unk, EnP, unkP,tHopp)	!projectBwf(qi, bWf, loBwf, U(:,:), failCount, smin, smax)
		!does the projection onto Loewdin-orthonormalized Bloch-like states
		!see Marzari, Vanderbilt PRB 56, 12847 (1997) Sec.IV.G.1 for detailed description of the method 
		real(dp),		intent(in)		:: En(:,:)		!	En(nBands,nQ)
		complex(dp)	,	intent(in)		:: unk(:,:,:)   ! unk(nR,nG,nQ)
		real(dp),		intent(out)		:: EnP(:,:)		! EnP(nWfs,nQ)
		complex(dp)	,	intent(out)		:: unkP(:,:,:), tHopp(:,:,:)	! unk(nR,nWfs,nQ) , tHopp(nWfs, nWfs nSC)
		complex(dp)	,	allocatable		:: loBwf(:,:), gnr(:,:), A(:,:), U(:,:)
		complex(dp)						:: phase
		integer							:: qi, xi
		!
		allocate(	loBwf(nR,nWfs)		)
		allocate(	gnr(nR,nWfs)		)
		allocate( 	A(nBands,nWfs)		)
		allocate(	U(nBands,nWfs)		)
		!
		if( debugProj ) then
			write(*,*)	"[projectUnk]: debug mode, will test unk normalization after debuging"
		end if

		if ( doProj ) then
			write(*,*)	"[projectUnk]: start with projections"
			!SET UP ROTATION MATRIX U
			call genTrialOrb(gnr)
			!PROJECT
			do qi = 1, nQ
				call calcAmat(qi,unk(:,:,qi),gnr, A) 
				!
				!USE SVD TO CALC U
				call calcUmat(A, U)
		
				!ROTATE BLOCH STATES
				!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(xi, phase)
				do xi = 1, nR
					phase		= myExp(	dot_product(qpts(:,qi),rpts(:,xi))		)
					loBwf(xi,:) = phase * matmul(	unk(xi,1:nBands,qi), U	)
				end do
				!$OMP END PARALLEL DO
				!
				!OVERWRITE UNKS WITH PROJECTED UNKs
				call genUnk(qi, lobWf, unkP(:,:,qi))

				!ADD CURRENT Q TO HOPPING MATRIX
				call addThopp(qi, U,En, tHopp)
			end do


			call calcEnergies(tHopp,EnP)
			!


		else
			write(*,*)	"[projectUnk]: projections disabled, use initial unk files"
			unkP(:,:,:)	= unk(:,1:nWfs,:)
			EnP(:,:)	= En(1:nWfs,:)
		end if
		write(*,*)	"[projectUnk]: done with projections at each k point"	
		!
		!
		!DEBUG
		if( debugProj ) then
			write(*,*)		"[projectUnk]: start test normalization of projected unks"
			if(.not. testNormUNK(unkP)	) then
				write(*,*)	"[projectUnk]: found normalization issues for projected unks"
			else
				write(*,*)	"[projectUnk]: no issues detected"
			end if
		end if
		write(*,*)			"[projectUnk]: finished debuging."
		!
		!
		return
	end subroutine
		



	subroutine addThopp(qi, U,En, tHopp)
		integer,		intent(in)		:: qi
		complex(dp),	intent(in)		:: U(:,:)	!	U(nBands,nWfs)
		real(dp),		intent(in)		:: En(:,:)	!	En(nBands,nQ)
		complex(dp),	intent(inout)	:: tHopp(:,:,:)	!tHopp(nWfs, nWfs nSC)
		complex(dp)						:: phase
		complex(dp),	allocatable		:: EU(:)
		integer							:: R,

		allocate(	EU(nWfs) )
		do R = 1, nSC
			phase			= myExp( - dot_product( qpts(:,qi), Rcell(:,R) ) 		)
			EU				= matmul(dcmplx(En(:,qi)), U(:,:))

			tHopp(:,:,R)	= tHopp(:,:,R) + phase * 
		end do
		return
	end subroutine






!privat:
	subroutine genTrialOrb(gnr)
		complex(dp),	intent(out)	:: gnr(:,:)       !gnr( nR, nWfs)	
		real(dp)					:: posX, posY, xc, k, L
		complex(dp)					:: A
		integer						:: n, ri, at
		!
		if(nWfs > nBands) then
			write(*,*)"[genTrialOrb]: warning, nWfs larger then nBands! No propper subspace..."
		end if

		gnr = dcmplx(0.0_dp)
		
		do ri = 1, nR
			!ONLY IN HOME UNIT CELL
			if(  rpts(1,ri) < aX .and. rpts(2,ri) < aY) then
				if( 	insideAt(1, rpts(:,ri)))	then
					gnr(ri,1)	= gVal(1,1,ri) 
					!gnr(ri,3)	= dcmplx(1.0_dp)
				end if
				!
				if( 	insideAt(2, rpts(:,ri))		)	then
					gnr(ri,2)	= gVal(1,1,ri)
					!gnr(ri,3)	= dcmplx(1.0_dp)
				end if
				!if( .not. insideAt(1, rpts(:,ri))	 .and.  .not. insideAt(2, rpts(:,ri))	) then
				!	gnr(ri,3)	= dcmplx(1.0_dp)
				!end if
			end if
		end do


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
				gVal = dcmplx(1.0_dp)!dcmplx(	trialOrbVAL(at) /sum( trialOrbVAL )		)
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



	subroutine calcAmat(qi, unk,gnr, A)
		!calculates the inner product matrix of bwfs and trial orbital gnr
		integer,		intent(in)	:: qi
		complex(dp),	intent(in)	:: unk(:,:), gnr(:,:) 
		complex(dp),	intent(out)	:: A(:,:) !A(nBands,nWfs)
		complex(dp),	allocatable	:: f(:)
		complex(dp)					:: phase
		integer						:: m,n, xi
		!
		A = dcmplx(0.0_dp)
		!
		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, m, f, xi)
		allocate(	f(nR)	)
		!$OMP DO COLLAPSE(2) SCHEDULE(STATIC) 
		do n = 1, nWfs
			do m = 1, nBands
				f = dcmplx(0.0_dp)
				do xi = 1, nR
					phase	= myExp( 	dot_product( qpts(:,qi), rpts(:,xi))		)				
					f(xi)	= dconjg(	phase * unk(xi,m) ) * gnr(xi,n)
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


	subroutine calcUmat(A, U)
		!U = A S^-0.5 = Z I V
		!	where A= Z d V (singular value decomposition) 
		complex(dp),	intent(in)		:: A(:,:)	!A(nBands,nWfs)
		complex(dp),	intent(out)		:: U(:,:)	!U(nBands,nWfs)
		complex(dp),	allocatable		:: Z(:,:), I(:,:), V(:,:), TMP(:,:)
		real(dp),		allocatable		:: d(:)
		integer							:: m, n, k, lda, ldb, ldc
		complex(dp)						:: alpha, beta
		character*1						:: transa, transb
		!
		allocate(	Z(	nBands, 	nBands	)	)
		allocate(	I(	nBands,		nWfs	)	)
		allocate(	V(	nWfs,		nWfs	)	)
		allocate(	TMP(nBands,		nWfs	)	)
		allocate(	d(	nWfs				)	)
		!
		!SET UP IDENTITY MATRIX
		I	= dcmplx(0.0_dp)
		do m = 1, nBands
			do n = 1, nWfs
				if( n == m) then
					I(m,n)	= dcmplx(1.0_dp)
				end if
			end do
		end do
		!
		!CALC SVD OF A
		call mySVD(A,Z,d,V)
		!
		!CALC U
		!tmp=matmul(I,V)
		transa	= 'n'
		transb	= 'n'
		alpha	= dcmplx(1.0_dp)
		beta	= dcmplx(0.0_dp)
		m		= size(I,1)
		n		= size(V,2)
		k		= size(I,2)
		lda		= m
		ldb		= k
		ldc		= m
		call zgemm(transa, transb, m, n, k, alpha, I, lda, V, ldb, beta, TMP, ldc)
		!
		!U=matmul(Z,tmp)
		transa	= 'n'
		transb	= 'n'
		alpha	= dcmplx(1.0_dp)
		beta	= dcmplx(0.0_dp)
		m		= size(Z,1)
		n		= size(TMP,2)
		k		= size(Z,2)
		lda		= m
		ldb		= k
		ldc		= m
		call zgemm(transa, transb, m, n, k, alpha, Z, lda, TMP, ldb, beta, U, ldc)
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
		!
		!
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
		!
		if( mi < smin ) then
			smin = mi
		end if
		if( ma > smax ) then
			smax = ma
		end if
		!
		!
		!DEBUG
		if(debugHam) then
			! test if S * (S^-0.5)^2 == I
			write(*,*)"[calcInvSmat]: all done, start debugging test"
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


	!logical function isOrthonorm(loBwf)
	!	!cheks wether loBwf is orthonormal by calculating the overlap matrix elements
	!	!	<psi_n,k|psi_n,k'> = nSC * \delta(n,m)
	!	!	WARNING: DOES NOT CHECK ORTHONORMALIZATION BETWEEN K POINTS
	!	complex(dp),	intent(in)		:: lobWf(:,:)   !lobWf(nR,nWfs)
	!	complex(dp),	allocatable		:: f(:)
	!	complex(dp)						:: oLap
	!	
	!	integer							:: n, m, ri
	!	logical							:: realist, imaginist
	!	allocate( f(nR)		)
	!	!
	!	isOrthonorm = .true.
	!	n 			= 1
	!	!
	!	do  while(n <= nWfs .and. isOrthonorm) 
	!		!DIAGONAL ELEMENTS
	!		!
	!		!
	!		!PREPARE INTEGRATION ARRAY
	!		do ri = 1, nR
	!			f(ri) = dconjg( lobWf(ri,n) ) * lobWf(ri,n)
	!		end do
	!		oLap	= nIntegrate(nR, nRx, nRy, dx, dy, f)
	!		oLap	= oLap / dcmplx( nSC	)
	!		!CHECK CONDITIONS
	!		realist 	= abs(			abs(dreal(oLap))		-		1.0_dp		) 	> acc
	!		imaginist	= abs(			dimag(oLap)								)	> acc
	!		if( realist .or. imaginist ) then
	!			write(*,'(a,i2,a,f15.10,a,f15.10)')"[isOrthonorm]: overlap(n,n=",n,") = ",dreal(oLap),"+i*",dimag(oLap)
	!			isOrthonorm = .false.
	!		end if
	!		!
	!		!
	!		!OFF DIAGONAL
	!		m = 1
	!		do while(m <= nWfs	.and. isOrthonorm)
	!			if(n /= m) then
	!				!PREPARE INTEGRATION ARRAY
	!				do ri = 1, nR
	!					f(ri) = dconjg( lobWf(ri,n) ) * lobWf(ri,m)
	!				end do
	!				oLap 	= nIntegrate(nR, nRx, nRy, dx, dy, f)
	!				oLap	= oLap / dcmplx( nSC	)
	!				!CHECK CONDITION
	!				if(abs(oLap) > acc) then
	!					write(*,'(a,i2,a,i2,a,f15.10,a,f15.10)')"[isOrthonorm]: overlap(",n,",",m,") = ",dreal(oLap),"+i*",dimag(oLap)
	!					isOrthonorm = .false.
	!				end if
	!			end if
	!			m = m +1
	!		end do
	!		n = n + 1
	!	end do
	!	!
	!	return
	!end function
		


end module projection