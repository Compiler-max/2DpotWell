module blochWf
	!generates bloch and lattice periodidc functions given a basCoeff matrix
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, acc, machineP,& 
								myExp, myLeviCivita, eigSolver, nIntegrate
	use sysPara

	implicit none

	private
	public	::	calcBasis, UNKoverlap


	contains







!public
	complex(dp) function UNKoverlap(n, m, qi, knb, gShift, ck)
		!calculates the overlap between unk at qi and at a neigbhouring k point knb
		!
		integer,		intent(in)		:: n, m, qi, knb
		real(dp),		intent(in)		:: gShift(2)
		complex(dp),	intent(in)		:: ck(:,:,:)  !ck(			nG		,	nBands  	,	nQ	)		
		integer							:: gi, gj, cnt
		real(dp)						:: delta(2)
		logical							:: notFound
		!
		UNKoverlap	= dcmplx(0.0_dp)
		cnt	= 0
		do gi = 1, nGq(qi)
			notFound 	= .true.
			gj			= 1
			do while( gj<= nGq(knb) .and. notFound ) 
				delta(:)	=  ( Gvec(:,gi,qi)-qpts(:,qi) ) 	-  		( Gvec(:,gj,knb)-qpts(:,knb)-gShift(:) )
				if( norm2(delta) < machineP )	then
					UNKoverlap	= UNKoverlap +  dconjg( ck(gi,n,qi) ) * ck(gj,m,knb) 
					cnt = cnt + 1
					notFound = .false.
				end if
				gj = gj + 1
			end do
		end do
		!
		if( cnt > nGq(qi)	)	write(*,'(a,i8,a,i8)')	"[UNKoverlap]: warning, used ",cnt," where nGmax(qi)=",nGq(qi)
		if( cnt < nGq(qi) / 2.0_dp)	write(*,'(a,i8,a,i8)')	"[UNKoverlap]: warning, used  only",cnt," where nGmax(qi)=",nGq(qi)
		!
		return
	end function


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









