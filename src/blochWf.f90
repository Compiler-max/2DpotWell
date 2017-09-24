module blochWf
	!generates bloch and lattice periodidc functions given a basCoeff matrix
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, acc, myExp, myLeviCivita, eigSolver, nIntegrate
	use sysPara

	implicit none

	private
	public	::	genBwfVelo, genUnk


	contains







!public
	subroutine genBwfVelo(qi,basCoeff, bWf, velobWf)
		!generates the bloch wavefunctions, with  the basCoeff from eigSolver, using
		!	call zgemm(transa, transb, m, n, k, alpha, a	  , lda, b		, ldb, beta, c , ldc)
		!			c = alpha * op(a) *op(b) + beta * c
		integer		, intent(in)	:: qi
		complex(dp)	, intent(in)	:: basCoeff(:,:)
		complex(dp)	, intent(out)	:: bWf(:,:,:), velobWf(:,:,:,:)	!bWf(nRpts,nG)			
		complex(dp)	, allocatable	:: basVec(:), veloBasX(:), veloBasY(:), tmp(:)
		integer 				 	:: xi, lda, ldb, ldc, m, n, k
		character*1					:: transa, transb
		complex(dp)					:: alpha, beta
		!
		
		!allocate(	tmp(		nG)		)
		!
		transa	= 'n'
		transb	= 'n'
		m		= 1
		n		= nG
		k 		= nG
		alpha	= dcmplx(1.0_dp)
		beta	= dcmplx(0.0_dp)
		lda		= size(basVec)
		ldb		= size(basCoeff,2)
		ldc		= size(tmp)
		!
		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(xi, basVec, veloBasX, veloBasY)
		allocate(	basVec(		nG)		)
		allocate(	veloBasX(	nG)		)
		allocate(	veloBasY(	nG)		)
		!$OMP DO SCHEDULE(STATIC) 
		do xi = 1, nR
			!GET BASIS
			call calcBasis(qi,xi, basVec, veloBasX, veloBasY)
			!
			!WAVE FUNCTIONS
			bwf(xi,:,qi)	= matmul(basVec,basCoeff) / dsqrt(vol)
			!
			!!VELOCITIES
			velobWf(1,xi,:,qi)	= matmul(veloBasX,basCoeff) / dsqrt(vol)
			velobWf(2,xi,:,qi)	= matmul(veloBasY,basCoeff) / dsqrt(vol)
		end do
		!$OMP END DO
		deallocate(	basVec		)
		deallocate(	veloBasX	)
		deallocate(	veloBasY	)
		!$OMP END PARALLEL
		!
		return 
	end subroutine




	subroutine genUnk(qi, bWf, unk)
		! generates the lattice periodic part from given bloch wave functions
		integer,		intent(in)		:: qi
		complex(dp),	intent(in)		:: bWf(:,:) !lobWf(	nR, nWfs)
		complex(dp),	intent(out)		:: unk(:,:,:)   !unk(	nR, nWfs)
		integer							:: xi, n
		complex(dp)						:: phase
		!
		!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2) DEFAULT(SHARED) PRIVATE(n, xi, phase) 
		do n = 1, nWfs
			do xi = 1, nR
				phase		 = myExp( -1.0_dp 	*	 dot_product( qpts(:,qi) , rpts(:,xi)	) 			)
				unk(xi,n,qi) = phase * bWf(xi,n)
			end do
		end do
		!$OMP END PARALLEL DO

		if(debugHam) then
			write(*,*)"[genUnk]: start normalization test"
			if( .not. testNormal(unk)	) then
				write(*,*)	"[genUnk]: unks have normalization issues"
			end if
		end if
		!
		return
	end subroutine


	logical function testNormal(unk)
		complex(dp),	intent(in)		:: unk(:,:,:)
		integer							:: q1, q2, m, n, ri, fcount
		complex(dp),	allocatable		:: f(:)
		complex(dp)						:: oLap
		logical							:: isNorm
		!
		fcount	= 0

		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(q1, q2, m, n, ri, f, oLap) REDUCTION(.AND.:isNorm) REDUCTION(+:fcount)
		allocate( f(nR) 	)
		!$OMP DO COLLAPSE(4), SCHEDULE(STATIC) 
		do q2 = 1, nQ
			do q1 = 1, nQ
				do m = 1, nWfs
					do n = 1, nWfs
						!INTEGRATE
						do ri = 1, nR
							f(ri)	= dconjg( unk(ri,m,q1)	)	* unk(ri,n,q2)
						end do
						oLap	= nIntegrate(nR, nRx,nRy, dx,dy, f)
						!ADJUST IF NONE ZERO
						if( n==m .and. q1==q2	) then
							oLap = oLap - dcmplx(nSC)
						end if
						!CHECK CONDITION
						if( abs(oLap) > acc ) then
							isNorm	= .false.
							fcount	= fcount + 1
						else
							isNorm	= .true.
						end if
					end do
				end do
			end do
		end do
		!$OMP END DO
		deallocate(	f	)
		!$OMP END PARALLEL

		!
		testNormal	= isNorm
		if( .not. testNormal) then
			write(*,'(a,i8,a,i8,a)')	"[testNormal]: ",fcount," of ",nQ**2+nG**2," tests of unk normalization failed"
		end if
		!
		return
	end function

	!subroutine testNormal(bwf)
	!	! <Y_nk1|Y_mk2> = N * \delta_n,m * \delta_k1,k2
	!	complex(dp),	intent(in)		:: bwf(:,:,:)
	!	complex(dp),	allocatable		:: f(:)
	!	integer							:: ri, q1,q2, n,m, count, tot
	!	complex(dp)						:: oLap
	!	real(dp)						:: avg
	!	!
	!	allocate(	f(nR)	)
	!	!
	!	count	= 0
	!	avg		= 0.0_dp
	!	tot		= 0
	!	!
	!	do m = 1, nWfs
	!		do n = 1, nWfs
	!			do q1 = 1, nQ
	!				!do q2 = 1, nQ
	!					!FILL INTEGRATION ARRAY
	!					do ri = 1, nR
	!						f(ri)	= dconjg(bwf(ri,q1,n)) * bwf(ri,q1,m)
	!					end do
	!					oLap	= nIntegrate(nR, nRx,nRy, dx,dy, f)
	!					!CHECK CONDITION
	!					if( dimag(oLap) > acc ) then
	!						count	= count + 1
	!						avg		= avg	+ abs(dreal(oLap)-nSC)
	!					else
	!						if(n==m .and. q1==q2) then
	!							if(abs(dreal(oLap)-nSC) > acc )then
	!								count	= count + 1
	!								avg		= avg	+ abs(dreal(oLap)-nSC)
	!							end if
	!						else
	!							if(abs(dreal(oLap)) > acc )then
	!								count	= count + 1
	!								avg		= avg	+ abs(dreal(oLap)-nSC)
	!							end if
	!						end if
	!					end if
	!					!
	!					!
	!					tot	= tot + 1
	!				!end do
	!			end do
	!		end do
	!	end do
!
!	!	avg	= avg / real(tot,dp)
!	!	write(*,*)"[testNormal]: found ",count," points of ",tot," not normalized bwfs, avg diff=",avg
!
!	!	return
	!end subroutine












!privat
	subroutine calcBasis(qi, ri, basVec, veloBasX, veloBasY)
		!calculates the basis vectors e^i(k+G).r
		!	if |k+G| is larger then the cutoff the basis vector is set to zero
		!	the cutoff enforces a symmetric base at each k point
		integer,	 intent(in)		:: qi, ri
		complex(dp), intent(out)	:: basVec(:), veloBasX(:), veloBasY(:)
		real(dp)				 	:: tmp(2)
		complex(dp)				 	:: phase
		integer 				 	::	i 
		!
		do i =1, nG
			tmp(:) = qpts(:,qi) + Gvec(:,i)
			!
			if( norm2(tmp) < Gcut ) then
				phase			= myExp( dot_product( tmp, rpts(:,ri) )		)
				basVec(i) 		= phase
				veloBasX(i) 	= phase * i_dp * (	qpts(1,qi) + Gvec(1,i)	)
				veloBasY(i)		= phase * i_dp * (	qpts(2,qi) + Gvec(2,i)	)
			else
				basVec(i) 		= dcmplx( 0.0_dp )
				veloBasX(i)		= dcmplx( 0.0_dp )
				veloBasY(i)		= dcmplx( 0.0_dp )
			end if
		end do
		!
		!
		return
	end subroutine






end module blochWf 













