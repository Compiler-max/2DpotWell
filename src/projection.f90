module projection
	!module contains all the functions to generate the wannier functions
	!	this includes the calculation of the orthonormalized projected bloch wavefunctions
	!	and the FT of the bwfs to calculate the wannier functions
	use omp_lib
	use mathematics, 	only: 	dp, PI_dp, acc, myExp, nIntegrate, eigSolver,  mySVD, myMatInvSqrt, isUnit, isIdentity, isHermitian
	use sysPara
	use blochWf,		only:	calcBasis, genUnk, testNormUNK
	use output,			only:	writeInterpBands
	implicit none	
	
	private
	public ::	projectUnk		














	contains
!public:

	subroutine projectUnk(ckH, ckW, Uq)	!projectBwf(qi, bWf, loBwf, U(:,:), failCount, smin, smax)
		!does the projection onto Loewdin-orthonormalized Bloch-like states
		!see Marzari, Vanderbilt PRB 56, 12847 (1997) Sec.IV.G.1 for detailed description of the method 
		complex(dp)	,	intent(in)		:: ckH(:,:,:)  ! unk(nR,nG,nQ)
		complex(dp)	,	intent(out)		:: ckW(:,:,:), Uq(:,:,:)	! Uq(nBands	,nWfs, 	nQ	)	
		real(dp),		allocatable		:: EnP(:,:)	
		complex(dp)	,	allocatable		:: loBwf(:,:), gnr(:,:), A(:,:), Ham(:,:), tmp(:,:)
		
		complex(dp)						:: phase
		integer							:: qi, xi, n, m, R
		!
		allocate(	loBwf(nR,nWfs)		)
		allocate(	gnr(nR,nWfs)		)
		allocate( 	A(nBands,nWfs)		)
		allocate(	Ham(nWfs,nWfs)		)
		allocate(	tmp(nBands,nWfs))
		allocate(	EnP(nWfs,nQ)		)
		!
		if( debugProj ) then
			write(*,*)	"[projectUnk]: debug mode, will test unk normalization after debuging"
		end if
		!
		!
		call genTrialOrb(gnr)
		!
		!PROJECT
		do qi = 1, nQ
			if( doProj ) then 
				Uq(:,:,qi)	= dcmplx(0.0_dp)
				call calcAmat(qi,ckH(:,:,qi) ,gnr, A) 
				call calcUmat(A, Uq(:,:,qi))
			else
				write(*,*)	"[projectUNK]: projection disabled, U= Identity"
				Uq(:,:,qi) = dcmplx(0.0_dp)
				do n = 1, size(Uq,1)
					do m = 1, size(Uq,2)
						if(n==m) then
							Uq(n,m,qi)	= dcmplx(1.0_dp)
						end if
					end do
				end do
			end if
			!
			!ROTATE BLOCH STATES
			ckW(:,:,qi)	= matmul( ckH(:,:,qi) , Uq(:,:,qi)	 ) !	


			!!!!!!ROTATE BLOCH STATES
			!!!!!!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(xi, phase)
			!!!!!do xi = 1, nR
			!!!!!	phase		= myExp(	dot_product(qpts(:,qi),rpts(:,xi))		)
			!!!!!	loBwf(xi,:) = matmul(	phase * unk(xi,1:nBands,qi) / dsqrt(real(nSC,dp)), Uq(:,:,qi)	)
			!!!!!end do
			!!!!!!$OMP END PARALLEL DO
			!!!!!!
			!!!!!!OVERWRITE UNKS WITH PROJECTED UNKs
			!!!!!call genUnk(qi, lobWf, unkP(:,:,qi))
		end do
		!

		write(*,*)	"[projectUnk]: done with projections at each k point"	
		!
		!
		!DEBUG
		!if( debugProj ) then
		!	write(*,*)		"[projectUnk]: start test normalization of projected unks"
		!	!if(.not. testNormUNK(unkP)	) then
		!	!	write(*,*)	"[projectUnk]: found normalization issues for projected unks"
		!	!else
		!	!	write(*,*)	"[projectUnk]: no issues detected"
		!	!end if
		!	!
		!	!
		!	if( .not. doProj ) then 
		!		write(*,*)		"[projectUnk]: start unk comparisson"
		!		call compareUNKs(unk,unkP)
		!	end if 
		!end if
		!write(*,*)			"[projectUnk]: finished debuging."
		!
		!
		return
	end subroutine
		



	subroutine addThopp(qi, U,En, tHopp)
		!	<0|H|R> = (1/sqrt(N)) \sum_q Vq^dagger Eq Vq
		!	where Eq is diagonal matrix (energies in Hamiltonian gauge)
		integer,		intent(in)		:: qi
		complex(dp),	intent(in)		:: U(:,:)	!	U(nBands,nWfs)
		real(dp),		intent(in)		:: En(:,:)	!	En(nBands,nQ)
		complex(dp),	intent(inout)	:: tHopp(:,:,:)	!tHopp(nWfs, nWfs nSC)
		complex(dp)						:: phase
		complex(dp),	allocatable		:: tmp(:,:),EV(:,:), Udag(:,:)
		integer							:: R, m, n
		!
		allocate(	tmp(nBands,nWfs)	)
		allocate(	EV(nWfs,nWfs) 		)
		allocate(	Udag(nWfs,nBands)	)
		!		
		!
		do R = 1, nSC
			!tmp = Eq Vq
			do m = 1, nWfs
				do n = 1, nBands
					tmp(n,m)	= dcmplx(En(n,qi)) * U(n,m)
				end do 
			end do
			! EU = Vq^dagger tmp
			Udag	= dconjg( transpose(U)	)
			EV		= matmul(Udag, tmp)
			!ADD results
			phase			= myExp( - dot_product( qpts(:,qi), Rcell(:,R) ) 		)
			tHopp(:,:,R)	= tHopp(:,:,R) + phase * EV(:,:) /  dcmplx(real(nSC,dp))
		end do
		!
		!tHopp	= tHopp    
		!
		return
	end subroutine


	




!privat:
	subroutine genTrialOrb(gnr)
		complex(dp),	intent(out)	:: gnr(:,:)       !gnr( nR, nWfs)	
		integer						:: n, ri, at, gammaP
		real(dp)					:: drX(2), drY(2)
		!
		if(nWfs > nBands) then
			write(*,*)"[genTrialOrb]: warning, nWfs larger then nBands! No propper subspace..."
		end if

		gammaP	= getGammaPoint()
		gnr 	= dcmplx(0.0_dp)
		

	


		!SINGLE ATOM
		!do n = 1, nWfs
		!	do ri = 1, nR
		!		if( 0.0_dp < rpts(1,ri) .and.  rpts(1,ri) < aX .and. 0.0_dp < rpts(2,ri) .and. rpts(2,ri) < aY) then
		!			if( 	insideAt(1, rpts(:,ri)))	then
		!				gnr(ri,n)	= potWell(1,n,ri)
		!			end if
		!		end if
		!	end do
		!end do

		!TWO ATOMS
		do n = 1, nWfs-1, 2
			do ri = 1, nR
				if( 0.0_dp < rpts(1,ri) .and.  rpts(1,ri) < aX .and. 0.0_dp < rpts(2,ri) .and. rpts(2,ri) < aY) then
					if( 	insideAt(1, rpts(:,ri)))	then
						gnr(ri,n)	= potWell(1,n,ri)
					else if( insideAt(2, rpts(:,ri)) ) then
						gnr(ri,n+1) = potWell(2,n,ri)
					end if
				end if
			end do
		end do		


		return
	end subroutine


	real(dp) function gaussian(rpt,cent, alpha, width)
		real(dp),		intent(in)		:: rpt(2), cent(2), alpha, width
		real(dp)						:: relP(2), val
		!
		relP(:)	= rpt(:) - cent(:)
		gaussian= alpha * exp(	dot_product(relP, relP)	  / (2.0_dp*width)	)
		!
		return
	end function


	logical function insideBond( xpt )
		real(dp),		intent(in)		:: xpt
		!
		!insideBond =  (atPos(1,1) - atR(1,1) <= xpt) .and. (xpt <= atPos(1,2)  + atR(1,2))
		insideBond =  (atPos(1,1)  <= xpt) .and. (xpt <= atPos(1,2) )
		!
		return
	end function




	logical function insideABond( xpt ) 
		real(dp),		intent(in)		:: xpt
		logical							:: left, right
		!
		left	= (  (atPos(1,1) - 1.0_dp*  atR(1,1) ) <=  xpt)	 .and. ( xpt <= 		atPos(1,1)					)
		right	= ( 		 atPos(1,2) 		<= 	xpt) .and. ( xpt <= ( atPos(1,2) + 1.0_dp* atR(1,2) )			)	
		!
		insideABond = left .or. right
		!
		return
	end function 



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



	complex(dp) function potWell(at, n, ri)
		integer, 	intent(in)		:: at, n, ri
		real(dp)					:: x,y, Lx, Ly, xc, yc, kx, ky, A
		!
		x	= rpts(1,ri)
		y	= rpts(2,ri)
		Lx 	= 2.0_dp * atR(1,at)
		Ly 	= 2.0_dp * atR(2,at)
		xc	= atPos(1,at)
		yc	= atPos(2,at)
		kx	= PI_dp/ Lx
		ky	= PI_dp/ Ly
		A	= dsqrt(4.0_dp / (Lx * Ly)	)
		!
		if( n==1 ) then
			potWell = A * dsin( kx*(x -xc + 0.5_dp * Lx) ) * dsin( ky*(y -yc + 0.5_dp * Ly) )
		else if( n==2 ) then
			potWell = A * dcos( kx*(x -xc + 0.5_dp * Lx) ) * dsin( ky*(y -yc + 0.5_dp * Ly) )
		else if( n==3 ) then
			potWell = A * dsin( kx*(x -xc + 0.5_dp * Lx) ) * dcos( ky*(y -yc + 0.5_dp * Ly) )
		end if 


		return
	end function



	complex(dp) function infPotWell(at, n, ri)
		integer, 	intent(in)		:: at, n, ri
		real(dp)					:: x,y, Lx, Ly, xc, yc, kx, ky, A
		!
		x	= rpts(1,ri)
		y	= rpts(2,ri)
		Lx 	= 2.0_dp * atR(1,at)
		Ly 	= 2.0_dp * atR(2,at)
		xc	= atPos(1,at)
		yc	= atPos(2,at)
		kx	= n*PI_dp/ Lx
		ky	= n*PI_dp/ Ly
		!
		A	= dsqrt(4.0_dp / (Lx * Ly)	)
		!
		infPotWell = A * dsin( kx*(x -xc + 0.5_dp * Lx) ) * dsin( ky*(y -yc + 0.5_dp * Ly) )
		!
		return
	end function






	subroutine calcAmat(qi,ckH ,gnr, A)
		!calculates the inner product matrix of bwfs and trial orbital gnr
		integer,		intent(in)	:: qi
		complex(dp),	intent(in)	:: ckH(:,:), gnr(:,:) 
		complex(dp),	intent(out)	:: A(:,:) !A(nBands,nWfs)
		complex(dp),	allocatable	:: psi(:), basVec(:)
		complex(dp)					:: phase
		integer						:: m,n, xi
		!
		allocate(	basVec(nG)	)
		allocate(	psi(nBands)	)
		!
		!
		!INTEGRATION OVER REAL SPACE
		A = dcmplx(0.0_dp)
		do xi = 1, nR
			!GET BASIS VECTOR
			call calcBasis(qi,xi, basVec)
			psi(:)	= matmul( basVec, ckH ) / dsqrt(vol)
			!CALC OVERLAP AT CURRENT GRID POINT
			do n = 1, nWfs
				do m = 1 , nBands 
					A(m,n)	= A(m,n) + psi(m) * gnr(xi,n)
				end do 
			end do
			!phase	= myExp( 	dot_product( qpts(:,qi), rpts(:,xi))		)				
			!f(xi)	= dconjg(	phase * unk(xi,m) / dsqrt(real(nSC,dp))  ) * gnr(xi,n)
		end do
		!NORMALIZE REAL SPACE INTEGRATION
		A	= A / real(nR,dp)
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






	subroutine compareUNKs(unk,unkP)
		complex(dp),	intent(in)		:: unk(:,:,:), unkP(:,:,:)
		integer							:: ri, n, qi, fCount, tCount
		real(dp)						:: avg
		!
		avg		= 0.0_dp
		fCount 	= 0
		tCount	= 0
		!
		do qi = 1, size(unkP,3)
			do n = 1, size(unkP,2)
				do ri = 1, size(unkP,1)
					if( abs( unkP(ri,n,qi)-unk(ri,n,qi)) > acc 	) then
						fCount 	= fCount + 1
						avg 	= avg + abs( unkP(ri,n,qi)-unk(ri,n,qi) )
					end if
					tCount = tCount +1
				end do
			end do
		end do
		!
		avg = avg / real(fCount,dp)
		write(*,*)"[compareUNKs]: comparisson of unks before and after unks where projected with Identity mat"
		write(*,'(a,i8,a,i12,a,f10.8)')"[compareUNKs]: ", fCount, " of ", tCount, " tests failed, average difference =",avg


		return
	end subroutine

















	!subroutine calcInvSmat(A, S, smin, smax)
	!	!calculates the sqrt inv. of the overlap matrix S
	!	!
	!	!	S	= A^dagger A
	!	complex(dp),	intent(in)		:: A(:,:)
	!	complex(dp),	intent(out)		:: S(:,:)
	!	real(dp),		intent(inout)	:: smin, smax
	!	complex(dp),	allocatable		:: Sold(:,:), Ssqr(:,:)
	!	integer							:: m,n,k,lda,ldb,ldc, i,j
	!	real(dp)						:: mi, ma
	!	complex(dp)						:: alpha, beta
	!	character*1						:: transa, transb
	!	!
	!	allocate(	Sold(nWfs,nWfs)	)
	!	!
	!	!
	!	!CALCULATE S FROM A 
	!	transa	= 'c'
	!	transb	= 'n'
	!	m		= size(A,1)
	!	n		= size(A,2)
	!	k		= size(A,2)
	!	alpha	= dcmplx(1.0_dp)
	!	beta	= dcmplx(0.0_dp)
	!	lda		= size(A,1)
	!	ldb		= size(A,1)
	!	ldc		= size(A,1)
	!	!call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
	!	 call zgemm(transa, transb, m, n, k, alpha, A, lda, A, ldb, beta, S, ldc)
	!	!
	!	!
	!	!CALCULATE INVERSE SQRT OF S
	!	Sold = S
!
!	!	call myMatInvSqrt(S, mi, ma)
!	!	!
!	!	if( mi < smin ) then
!	!		smin = mi
!	!	end if
!	!	if( ma > smax ) then
!	!		smax = ma
!	!	end if
!	!	!
!	!	!
!	!	!DEBUG
!	!	if(debugHam) then
!	!		! test if S * (S^-0.5)^2 == I
!	!		write(*,*)"[calcInvSmat]: all done, start debugging test"
!	!		allocate(	Ssqr(nWfs,nWfs)	)
!	!		transa	= 'n'
!	!		transb	= 'n'
!	!		alpha	= dcmplx(1.0_dp)
!	!		beta	= dcmplx(0.0_dp)
!	!		call zgemm(transa, transb, m, n, k, alpha, S   , lda, S   , ldb, beta, Ssqr, ldc)
!	!		call zgemm(transa, transb, m, n, k, alpha, Sold, lda, Ssqr, ldb, beta, Sold, ldc)
!	!		if( .not. isIdentity(Sold) ) then
!	!			write(*,*)"[calcInvSmat]: problem with mat inversion, seems to be not inverse square root"
!	!		end if
!	!	end if
!	!	!
!	!	!
!	!	return
!	!end subroutine
!

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