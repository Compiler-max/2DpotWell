module blochWf
	!generates bloch and lattice periodidc functions given a basCoeff matrix
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, acc, myExp, myLeviCivita, eigSolver, nIntegrate
	use sysPara

	implicit none

	private
	public	::	calcBasis


	contains







!public
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







!privat
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









