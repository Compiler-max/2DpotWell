module blochWf
	!generates bloch and lattice periodidc functions given a basCoeff matrix
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, acc, myExp, myLeviCivita, eigSolver, nIntegrate
	use sysPara

	implicit none

	private
	public	::	genBwfVelo, genUnk, calcBasis


	contains







!public
	subroutine genBwfVelo(qi,basCoeff, unk)
		!generates the bloch wavefunctions, with  the basCoeff from eigSolver, using
		!	call zgemm(transa, transb, m, n, k, alpha, a	  , lda, b		, ldb, beta, c , ldc)
		!			c = alpha * op(a) *op(b) + beta * c
		integer		, intent(in)	:: qi
		complex(dp)	, intent(in)	:: basCoeff(:,:)
		complex(dp)	, intent(out)	:: unk(:,:)
		complex(dp)	, allocatable	:: basVec(:) 
		integer 				 	:: xi
		complex(dp)					:: phase
		!

		!
		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(xi,phase, basVec)
		allocate(	basVec(		nG)			)
		!$OMP DO SCHEDULE(DYNAMIC,nRx) 
		do xi = 1, nR
			!GET BASIS
			call calcBasis(qi,xi, basVec)
			!
			!WAVE FUNCTIONS
			phase			= myExp( -1.0_dp * dot_product( qpts(:,qi), rpts(:,xi) )		)
			unk(xi,:)	= phase * dsqrt(real(nSc,dp))   * matmul(basVec,basCoeff) 
			!unk(xi,:)	= phase   * matmul(basVec,basCoeff) 
			!unk(xi,:)	= matmul(basVec,basCoeff) 
			!
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		!
		return 
	end subroutine




	subroutine genUnk(qi, bWf, unk)
		! generates the lattice periodic part from given bloch wave functions
		integer,		intent(in)		:: qi
		complex(dp),	intent(in)		:: bWf(:,:) !lobWf(	nR, nWfs)
		complex(dp),	intent(out)		:: unk(:,:)   !unk(	nR, nWfs)
		integer							:: xi, n
		complex(dp)						:: phase
		!
		!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2) DEFAULT(SHARED) PRIVATE(n, xi, phase) 
		do n = 1, nWfs
			do xi = 1, nR
				phase		 = myExp( -1.0_dp 	*	 dot_product( qpts(:,qi) , rpts(:,xi)	) 			)
				unk(xi,n) = phase * dsqrt(real(nSC,dp))* bWf(xi,n)
			end do
		end do
		!$OMP END PARALLEL DO
		!
		return
	end subroutine




!privat
	subroutine calcBasis(qi, ri, basVec)
		!calculates the basis vectors e^i(k+G).r
		!	if |k+G| is larger then the cutoff the basis vector is set to zero
		!	the cutoff enforces a symmetric base at each k point
		integer,	 intent(in)		:: qi, ri
		complex(dp), intent(out)	:: basVec(:)
		integer 				 	:: i 
		!
		basVec	= 0.0_dp
		do i =1, nGq(qi)
			basVec(i) 		= myExp( dot_product( Gvec(:,i,qi), rpts(:,ri) )		)  !/ dsqrt(vol)
		end do
		!
		!
		return
	end subroutine



	subroutine testNormalBwf(qi, bwf)
		integer,		intent(in)		:: qi
		complex(dp),	intent(in)		:: bwf(:,:)
		complex(dp),	allocatable		:: f(:)
		integer							:: n, xi,yi,ri, cnt
		complex(dp)						:: oLap
		!
		allocate(	f(100)	)
		write(*,*)"[testNormalBwf]: qi=",qi
		do n = 1, nG
			!INTEGRATE
			cnt = 1
			do xi = 1, 10
				do yi = 1, 10
					ri		= getRindex(xi,yi)
					f(cnt)	= dconjg( bwf(ri,n) ) * bwf(ri,n)
					cnt		= cnt + 1
				end do
			end do
			oLap	= nIntegrate(100, 10, 10, 0.1_dp, 0.1_dp, f)
			!CONDITION
			!
			write(*,'(a,i3,a,f6.4a,f6.4,a)')"[testNormalBwf]: n=",n,", oLap=(",dreal(oLap),"+i*",dimag(oLap),")."		
			!write(*,'(a,i3,a,f6.4a,f6.4,a)')"[testNormalBwf]: n=",n,", nSC*oLap=(",nSC*dreal(oLap),"+i*",nSC*dimag(oLap),")."			
		end do

		return
	end subroutine

end module blochWf 










	!logical function testNormBWF(bwf)
	!	complex(dp),	intent(in)		:: bwf(:,:,:)
	!	integer							:: q1, q2, m, n, ri, fcount, cnt
	!	complex(dp),	allocatable		:: f(:)
	!	complex(dp)						:: oLap, lphase, rphase
	!	logical							:: isNorm
	!	!
	!	fcount	= 0
!
!	!	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(q1, q2, m, n, ri, f, oLap, lphase, rphase) REDUCTION(.AND.:isNorm) REDUCTION(+:fcount, cnt)
!	!	allocate( f(nR) 	)
!	!	!$OMP DO COLLAPSE(4), SCHEDULE(STATIC) 
!	!	do q2 = 1, nQ
!	!		do q1 = 1, nQ
!	!			do m = 1, nWfs
!	!				do n = 1, nWfs
!	!					!INTEGRATE
!	!					do ri = 1, nR
!	!						lphase	= myExp( dot_product( qpts(:,q1), rpts(:,ri))	)
!	!						rphase	= myExp( dot_product( qpts(:,q2), rpts(:,ri))	)
!	!						f(ri)	= dconjg( lphase * unk(ri,m,q1) 	)	* rphase * unk(ri,n,q2)
!	!					end do
!	!					oLap	= nIntegrate(nR, nRx,nRy, dx,dy, f)
!	!					!ADJUST IF NONE ZERO
!	!					!if( n==m .and. q1==q2	) then
!	!					!	oLap = oLap - dcmplx(nSC)
!	!					!end if
!	!					!CHECK CONDITION
!	!					if( abs(oLap) > acc ) then
!	!						isNorm	= .false.
!	!						fcount	= fcount + 1
!	!						write(*,'(a,i2,a,i2,a,i3,a,i3,a,e10.3,a,e10.3)')	"[testNormal]: n=",n,",m=",m,", q1=",q1,", q2=",q2," oLap =",dreal(oLap),"+i*",dimag(oLap)
!	!					else
!	!						isNorm	= .true.
!	!					end if
!	!					cnt = cnt + 1
!	!				end do
!	!			end do
!	!		end do
!	!	end do
!	!	!$OMP END DO
!	!	deallocate(	f	)
!	!	!$OMP END PARALLEL
!
!	!	!
!	!	testNormUNK	= isNorm
!	!	if( .not. testNormUNK) then
!	!		write(*,'(a,i8,a,i8,a)')	"[testNormUNK]: ",fcount," of ",cnt," tests of unk normalization failed"
!	!	end if
!	!	!
!	!	return
	!end function






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






!	subroutine FDvelocities(unk, velo)
!		complex(dp),	intent(in)		:: unk(:,:,:)		!unk(nR,nStates, nQ)
!		complex(dp),	intent(out)		:: velo(:,:,:,:)	!velo(3,nStates,nStates,nQ)
!		complex(dp),	allocatable		:: f(:)
!		integer							:: qi, n, m, xi, yi, ri, rNN
!		!
!		velo = dcmplx(0.0_dp)
!		
!		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(f, qi, n, m , xi, yi, ri, rNN)
!		allocate(	f(nR)	)
!		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!		do qi = 1, size(unk,3)
!			do n = 1, size(velo,2)
!				do m = 1, size(velo,3)
!					!
!					!X DERIVATIVE
!					do xi = 1, nRx
!						do yi = 1, nRy
!							ri 	= getRindex(xi,yi)
!							if(xi /= nRx) then
!								rNN	= getRindex(xi+1,yi) 
!							else
!								rNN	= getRindex(1,yi)
!							end if
!							!write(*,'(a,f6.3,a,f6.3,a,f6.3,a,f6.3,a)')		"[FDvelocities]: ri= (",rpts(1,ri),",",rpts(2,ri),&
!							!						") ,rNN= (",rpts(1,rNN),",",rpts(2,rNN),")."
!							f(ri)	= dconjg(unk(ri,n,qi)) * ( qpts(1,qi) * unk(ri,m,qi) - i_dp * FD(qi,m,ri,rNN, unk) )
!						end do
!					end do
!					velo(1,n,m,qi) = nIntegrate(nR, nRx, nRy, dx, dy, f)
!					!
!					!
!					!Y DERIVATIVE
!					do xi = 1, nRx
!						do yi = 1, nRy
!							ri 	= getRindex(xi,yi)
!							if(yi /= nRy) then
!								rNN	= getRindex(xi,yi+1) 
!							else
!								rNN	= getRindex(xi,1)
!							end if
!							f(ri)	= dconjg(unk(ri,n,qi)) * ( qpts(2,qi) * unk(ri,m,qi) - i_dp * FD(qi,m,ri,rNN, unk) )
!						end do
!					end do
!					velo(2,n,m,qi) = nIntegrate(nR, nRx, nRy, dx, dy, f)
!					!
!				end do 
!			end do				
!		end do
!		!$OMP END DO
!		deallocate(	f	)
!		!$OMP END PARALLEL
!		!
!		!
!		return
!	end subroutine

!	complex(dp) function FD(qi, m, ri, rNN, unk)
!		!finite difference between ri and rNN
!		integer,		intent(in)		:: qi, m, ri, rNN
!		complex(dp),	intent(in)		:: unk(:,:,:)
!		real(dp)						:: h
!		!
!		h	= norm2( rpts(:,ri) - rpts(:,rNN)	)
!		FD	= unk(rNN,m,qi) - unk(ri,m,qi) / dcmplx(h)
!		!
!		return
!	end function









