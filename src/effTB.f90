module effTB
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP, myExp, myLeviCivita, nIntegrate, isHermitian
	use sysPara
	use blochWf,		only:	UNKoverlap
	use output,			only:	writeInterpU

	implicit none

	private
	public	::					TBviaKspace, calcConnOnCoarse
	contains	


	subroutine TBviaKspace(ckQ, EnQ, Uq, tHopp, rHopp)
		complex(dp),		intent(in)		:: ckQ(:,:,:), Uq(:,:,:)
		real(dp),			intent(in)		:: EnQ(:,:)
		complex(dp),		intent(out)		:: tHopp(:,:,:), rHopp(:,:,:,:)
		complex(dp),		allocatable		:: Atmp(:,:,:,:)
		complex(dp),		allocatable		:: Htmp(:,:,:)
		complex(dp)							:: phase
		integer								:: R, qi
		!
		allocate(	Atmp(	2,		nWfs,	nWfs,	nQ	)			)
		allocate(	Htmp(			nWfs,	nWfs,	nQ	)			)
		
		!SET UP K SPACE QUANTITIES
		call calcConnOnCoarse(ckQ, Atmp)
		call calcHtmp(EnQ, Uq, Htmp)
		!FT TO REAL SPACE
		rHopp	= dcmplx(0.0_dp)
		tHopp	= dcmplx(0.0_dp)
		do R = 1, nSC
			do qi = 1, nQ
				phase			= myExp( -1.0_dp * dot_product(qpts(1:2,qi),Rcell(1:2,R))		)  /dcmplx(size(tHopp,3))!/ dsqrt(real(nSC,dp) )
				!
				!RHOPP
				rHopp(1,:,:,R)	= rHopp(1,:,:,R) 	+ phase * Atmp(1,:,:,qi) 
				rHopp(2,:,:,R)	= rHopp(2,:,:,R) 	+ phase * Atmp(2,:,:,qi) 
				!
				!THOPP
				tHopp(:,:,R)	= tHopp(:,:,R)		+ phase * Htmp(:,:,qi)
			end do
		end do
		!
		!
		return
	end subroutine





!private

	subroutine calcHtmp(EnQ, Uq, Htmp)
		!Slater Koster Interpolation of tight binding paramters
		!	Htmp = Uq^dag Eq Uq (wann review)
		real(dp),		intent(in)		:: EnQ(:,:)
		complex(dp),	intent(in)		:: Uq(:,:,:)
		complex(dp),	intent(out)		:: Htmp(:,:,:)
		complex(dp),	allocatable		:: Udag(:,:), U(:,:), Ediag(:,:), tmp(:,:)
		integer							:: qi, n
		!
		allocate(	U(		nWfs, nwfs	)		)
		allocate(	tmp(	nWfs, nWfs	)		)
		allocate(	Udag(	nWfs, nwfs	)		)
		allocate(	Ediag(	nWfs, nWfs	)		)
		Htmp	= dcmplx(0.0_dp)
		!
		!
		!SET UP
		do qi = 1, nQ
			U(:,:)	= Uq(:,:,qi)
			Udag	= dconjg( transpose( U ) )
			Ediag	= dcmplx(0.0_dp)
			do n = 1, nWfs
				Ediag(n,n)	= EnQ(n,qi)
			end do
			tmp				= matmul(	Ediag	, 	U		)
			Htmp(:,:,qi)	= matmul( 	Udag	, 	tmp		)
			!
			if( .not. isHermitian(Htmp(:,:,qi)) ) 	then
				write(*,'(a,i3)')	"[TBviaKspace]: generated H matrix is not hermitian, at qi=",qi
			end if
		end do
		!
		!
		return
	end subroutine

		
	subroutine calcConnOnCoarse(ck, A)
		!finite difference on lattice periodic unk to calculate the Berry connection A
		!	A_n(k) 	= <u_n(k)|i \nabla_k|u_n(k)>
		!		 	= i  <u_n(k)| \sum_b{ w_b * b * [u_n(k+b)-u_n(k)]}
		!			= i \sum_b{		w_b	 * [  <u_n(k)|u_n(k+b)> -  <u_n(k)|u_n(k)>]		}
		!
		! see Mazari, Vanderbilt PRB.56.12847 (1997), Appendix B
		!
		complex(dp),	intent(in)		:: ck(:,:,:) 	! ckW(nG, nWfs, nQ)		
		complex(dp),	intent(out)		:: A(:,:,:,:)			
		complex(dp)						:: Mxl, Mxr, Myl, Myr, one
		integer							:: n, m, Z, qi, qx, qy, qxl, qxr, qyl, qyr, found, tot, al, be
		real(dp)						:: wbx,wby, bxl(2), bxr(2), byl(2), byr(2),dmax, avg, delta, &
											Gxl(2), Gyl(2), Gxr(2), Gyr(2), zero(2)
		!
		A 		= dcmplx(0.0_dp)
		Z 		= 4	!amount of nearest neighbours( 2 for 2D cubic unit cell)
		wbx 	= 2.0_dp / 		( real(Z,dp) * dqx**2 )
		wby		= wbx
		write(*,*)	"[calcConnOnCoarse]	weight wb=", wbx
		!wby 	= 1.0_dp /		( real(Z,dp) * dqy**2 )
		!b vector two nearest X neighbours:
		bxl(1) 	= -dqx				
		bxl(2)	= 0.0_dp
		bxr(1) 	= +dqx
		bxr(2)	= 0.0_dp
		!b vector two nearest Y neighbours:
		byl(1) 	= 0.0_dp
		byl(2)	= -dqy
		byr(1) 	= 0.0_dp
		byr(2)	= +dqy
		!
		!DEBUG WEIGHTS
		do al = 1, 2
			do be = 1, 2
				delta 	= 0.0_dp
				delta = delta + wbx * bxl(al) * bxl(be)
				delta = delta + wbx * bxr(al) * bxr(be)
				delta = delta + wby * byl(al) * byl(be)
				delta = delta + wby * byr(al) * byr(be)
				if( al==be .and. abs(delta-1.0_dp) > acc ) then
					write(*,'(a,i1,a,i1,a,f6.3)') &
							"[calcConnCoarse]: weights dont fullfill condition for a=",al," b=",be," delta=",delta
				else if ( al/=be .and. abs(delta) > acc ) then
					write(*,'(a,i1,a,i1,a,f6.3)') & 
							"[calcConnCoarse]: weights dont fullfill condition for a=",al," b=",be,"delta=",delta
				end if
			end do
		end do

		!
		found	= 0
		tot		= 0
		dmax	= 0.0_dp
		avg		= 0.0_dp
		!
		write(*,'(a,f6.3,a,f6.3)')	"[calcConnOnCoarse]: dqx=",dqx," dqy=",dqy
		!
		!$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) SCHEDULE(STATIC) &
		!$OMP& PRIVATE(m,n,qx,qy, qxl, qxr, qyl, qyr, qi,one, Mxl, Mxr, Myl, Myr, Gxl, Gxr, Gyl, Gyr, zero)
		do m = 1, nWfs
			do n = 1, nWfs
				do qx = 1, nQx
					do qy = 1, nQy
						!GET NEIGHBOURS
						qxl	= getLeft( qx,nQx)
						qxr	= getRight(qx,nQx)
						qyl	= getLeft(  qy,nQy)
						qyr = getRight( qy,nQy)

						!GET GRID POSITION OF NEIGHBOURS
						qi	= getKindex(qx,qy)
						qxl	= getKindex(qxl,qy)
						qxr	= getKindex(qxr,qy)
						qyl	= getKindex(qx,qyl)
						qyr	= getKindex(qx,qyr)
						!call testNeighB(qi, qxl, qxr, qyl, qyr)
						!
						!SHIFT NEIGHBOURS BACK TO FIRST BZ
						Gxl(:)	= 0.0_dp
						Gxr(:)	= 0.0_dp
						Gyl(:)	= 0.0_dp
						Gyr(:)	= 0.0_dp
						zero(:)	= 0.0_dp
						if( qx == 1 ) 	Gxl(1)	= - 2.0_dp * PI_dp / aX
						if( qx == nQx)	Gxr(1)	= + 2.0_dp * PI_dp / aX
						if( qy == 1 ) 	Gyl(2)	= - 2.0_dp * PI_dp / aY
						if( qy == nQy)	Gyr(2)	= + 2.0_dp * PI_dp / aY

						!write(*,*)"*"
						!write(*,*)"*"
						!write(*,*)"*"
						!write(*,'(a,f6.3,a,f6.3,a)')	"[calcConnOnCoarse]: q_i=(",qpts(1,qi) ,", ",qpts(2,qi) ,")"
						!write(*,'(a,f6.3,a,f6.3,a)')	"[calcConnOnCoarse]: qxl=(",qpts(1,qxl),", ",qpts(2,qxl),")"
						!write(*,'(a,f6.3,a,f6.3,a)')	"[calcConnOnCoarse]: qxr=(",qpts(1,qxr),", ",qpts(2,qxr),")"
						!write(*,'(a,f6.3,a,f6.3,a)')	"[calcConnOnCoarse]: qyl=(",qpts(1,qyl),", ",qpts(2,qyl),")"
						!write(*,'(a,f6.3,a,f6.3,a)')	"[calcConnOnCoarse]: qyr=(",qpts(1,qyr),", ",qpts(2,qyr),")"
						!write(*,*)"*"
						!write(*,'(a,f6.3,a,f6.3,a)')"[calcConnOnCoarse]:  Gxl=",Gxl(1),", ",Gxl(2),")."
						!write(*,'(a,f6.3,a,f6.3,a)')"[calcConnOnCoarse]:  Gxr=",Gxr(1),", ",Gxr(2),")."
						!write(*,'(a,f6.3,a,f6.3,a)')"[calcConnOnCoarse]:  Gyl=",Gyl(1),", ",Gyl(2),")."
						!write(*,'(a,f6.3,a,f6.3,a)')"[calcConnOnCoarse]:  Gyr=",Gyr(1),", ",Gyr(2),")."
						!
						!OVERLAP TO NEAREST NEIGHBOURS
						one	= UNKoverlap(	n,		m,		qi		, 	qi		,	zero	, ck	)
						Mxl	= UNKoverlap(	n,		m, 		qi		,	qxl 	,	Gxl		, ck	) 
						Mxr	= UNKoverlap(	n,		m, 		qi		,	qxr		,	Gxr		, ck	)
						Myl	= UNKoverlap(	n,		m, 		qi		, 	qyl		,	Gyl		, ck	)
						Myr	= UNKoverlap(	n,		m, 		qi		, 	qyr		,	Gyr		, ck	)

						if(		 n==m .and. abs(abs(one)-1.0_dp) > acc 	) then		!.and.			 abs(one-dcmplx(1.0_dp)) > acc ) then
							write(*,*)"[calcConnOnCoarse]: normalization issue detected"
							write(*,'(a,i2,a,f6.3,a,f6.3)') "[calcConnOnCoarse]: warning n=m=",n," one=",dreal(one),"+i*",dimag(one)
							
						else if( n/=m .and. abs(one) > acc ) then
							write(*,*)"[calcConnOnCoarse]: normalization issue detected"
							write(*,'(a,f15.6,a,f15.6)') "[calcConnOnCoarse]:  zero=",dreal(one),"+i*",dimag(one)
						end if


						!FD SUM OVER NEAREST NEIGHBOURS
						!A(1:2,n,m, qi) = A(1:2,n,m, qi) + i_dp * wbx * bxl(1:2) * ( Mxl - one )
						!A(1:2,n,m, qi) = A(1:2,n,m, qi) + i_dp * wbx * bxr(1:2) * ( Mxr - one )
						!A(1:2,n,m, qi) = A(1:2,n,m, qi) + i_dp * wby * byl(1:2) * ( Myl - one )
						!A(1:2,n,m, qi) = A(1:2,n,m, qi) + i_dp * wby * byr(1:2) * ( Myr - one )


						A(1:2,n,m, qi) = A(1:2,n,m, qi) - wbx * bxl(1:2) * dimag( zlog( Mxl ) )
						A(1:2,n,m, qi) = A(1:2,n,m, qi) - wbx * bxr(1:2) * dimag( zlog( Mxr ) )
						A(1:2,n,m, qi) = A(1:2,n,m, qi) - wby * byl(1:2) * dimag( zlog( Myl ) )
						A(1:2,n,m, qi) = A(1:2,n,m, qi) - wby * byr(1:2) * dimag( zlog( Myr ) )
						
					end do
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		!
		!
		return
	end subroutine


	subroutine testNeighB(qi, qxl, qxr, qyl, qyr)
		integer,		intent(in)		:: qi, qxl, qxr, qyl, qyr
		!
		!
		!X LEFT
		if( 	norm2(qpts(:,qi)-qpts(:,qxl)) > dqx+machineP		 ) then
			write(*,'(a,i3,a,f6.3,a,f6.3,a,a,f6.3,a,f6.3,a)')	&
						"[testNeighB]: problem with x left  at qi=",qi,&
						", qi=(",qpts(1,qi),", ",qpts(2,qi),")",&
						", qxl=(",qpts(1,qxl),", ",qpts(2,qxl),")."
		end if
		!
		!X RIGHT
		if( 	norm2(qpts(:,qxr)-qpts(:,qi)) > dqx+machineP 		 ) then
			write(*,'(a,i3,a,f6.3,a,f6.3,a,a,f6.3,a,f6.3,a)')	&
					"[testNeighB]: problem with x right at qi=",qi,&
					", qi=(",qpts(1,qi),", ",qpts(2,qi),")",&
					", qxr=(",qpts(1,qxr),", ",qpts(2,qxr),")."
		end if
		!
		!
		!Y LEFT
		if( norm2(qpts(:,qi)-qpts(:,qyl)) > dqy+machineP  		 ) then
			write(*,'(a,i3,a,f6.3,a,f6.3,a,a,f6.3,a,f6.3,a)')	&
					"[testNeighB]: problem with y left  at qi=",qi,&
					", qi=(",qpts(1,qi),", ",qpts(2,qi),")",&
					", qyl=(",qpts(1,qyl),", ",qpts(2,qyl),")."
		end if
		!
		!Y RIGHT
		if( 	norm2(qpts(:,qyr)-qpts(:,qi)) > dqy+machineP  		 ) then
			write(*,'(a,i3,a,f6.3,a,f6.3,a,a,f6.3,a,f6.3,a)')	& 
					"[testNeighB]: problem with y right at qi=",qi,&
					", qi=(",qpts(1,qi),", ",qpts(2,qi),")",&
					", qyr=(",qpts(1,qyr),", ",qpts(2,qyr),")."
		end if
		write(*,*)"*"
		write(*,*)"*"
		write(*,*)"*"
		!
		!
		return
	end subroutine




	integer function getLeft(i,N)
		!HELPER for calcConn
		!gets left (lower) neighbour, using the periodicity at boundary
		!
		integer,	intent(in)	:: i,N
		if(i.eq.1) then
			getLeft = N
		else
			getLeft = i-1
		end if
		!
		return
	end function


	integer function getRight(i,N)
		!HELPER for calcConn
		!gets right (upper) neighbour, using the periodicity at boundary
		!
		integer,	intent(in)	:: i,N
		if(i.eq.N) then
			getRight = 1
		else
			getRight = i+1
		end if
		!
		return
	end function



	





	subroutine testIfReal(rImag, tImag)
		complex(dp),	intent(in)		:: rImag(:,:,:,:), tImag(:,:,:)  !rImag(2,nWfs,nwfs,nSC), timag(nWfs,nWfs,nSC))
		integer							:: R,n,m, a, cntR, cntT, tot				
		real(dp)						:: avgR, avgT
		!
		cntT	= 0
		cntR	= 0
		tot		= 0
		avgT	= 0.0_dp
		avgR	= 0.0_dp
		!
		do R = 1, nSC
			do m = 1, nWfs
				do n = 1, nWfs
					if(dimag(tImag(n,m,R)) > machineP ) 	then
						cntT	= cntT + 1
						avgT	= avgT + dimag(tImag(n,m,R))
						write(*,*) 	"[testIfReal]: dimag(tImag)=",dimag(tImag(n,m,R))
					end if
					do a = 1, 2
						if(dimag(rImag(a,n,m,R)) > machineP ) 	then
							cntR	= cntR + 1
							avgR	= avgR + dimag(rImag(a,n,m,R))
							write(*,*) 	"[testIfReal]: dimag(rImag)=",dimag(rImag(a,n,m,R))
						end if
					end do
					tot	= tot + 1
				end do
			end do
		end do
		!
		avgR	= avgR / cntR
		avgT	= avgT / cntT
		!
		write(*,'(a,i4,a,i8,a,f8.4)')	"[testIfReal]: ",cntT," of ",tot," tests for tImag failed, avg=",avgT
		write(*,'(a,i4,a,i8,a,f8.4)')	"[testIfReal]: ",cntR," of ",tot," tests for rImag failed, avg=",avgR
		!
		!
		return
	end subroutine




end module
