module effTB
use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP, myExp, myLeviCivita, nIntegrate, isHermitian
	use sysPara


	implicit none

	private
	public	::					TBviaKspace
	contains	


	subroutine TBviaKspace(unkQ, EnQ, Uq, tHopp, rHopp)
		complex(dp),		intent(in)		:: unkQ(:,:,:), Uq(:,:,:)
		real(dp),			intent(in)		:: EnQ(:,:)
		complex(dp),		intent(out)		:: tHopp(:,:,:), rHopp(:,:,:,:)
		real(dp),			allocatable		:: AconnQ(:,:,:,:)
		complex(dp),		allocatable		:: Htmp(:,:,:)
		complex(dp)							:: phase
		integer								:: R, qi
		!
		allocate(	AconnQ(	2,	nWfs,	nWfs,	nQ		)			)
		allocate(	Htmp(			nWfs,	nWfs,	nQ		)			)
		
		!SET UP K SPACE QUANTITIES
		call calcConnOnCoarse(unkQ, AconnQ)
		call calcHtmp(EnQ, Uq, Htmp)
		!FT TO REAL SPACE
		rHopp	= dcmplx(0.0_dp)
		tHopp	= dcmplx(0.0_dp)
		do R = 1, nSC
			do qi = 1, nQ
				phase			= myExp( -1.0_dp * dot_product(qpts(:,qi),Rcell(:,R))		)  / dsqrt(real(nQ,dp) )
				!
				!RHOPP
				rHopp(1,:,:,R)	= rHopp(1,:,:,R) 	+ phase * AconnQ(1,:,:,qi) 
				rHopp(2,:,:,R)	= rHopp(2,:,:,R) 	+ phase * AconnQ(2,:,:,qi) 
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
		real(dp),		intent(in)		:: EnQ(:,:)
		complex(dp),	intent(in)		:: Uq(:,:,:)
		complex(dp),	intent(out)		:: Htmp(:,:,:)
		complex(dp),	allocatable		:: Udag(:,:), U(:,:), Ediag(:,:)
		integer							:: qi, n
		!
		allocate(	U(		nWfs, nwfs	)		)
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
			Htmp(:,:,qi)	= matmul(	Ediag	, 	U				)
			Htmp(:,:,qi)	= matmul( 	Udag	, 	Htmp(:,:,qi) 	)
			!
			if( .not. isHermitian(Htmp(:,:,qi)) ) 	then
				write(*,'(a,i3)')	"[TBviaKspace]: generated H matrix is not hermitian, at qi=",qi
			end if
		end do
		!
		!
		return
	end subroutine

		
	subroutine calcConnOnCoarse(unk, A)
		!finite difference on lattice periodic unk to calculate the Berry connection A
		!	A_n(k) 	= <u_n(k)|i \nabla_k|u_n(k)>
		!		 	= i  <u_n(k)| \sum_b{ w_b * b * [u_n(k+b)-u_n(k)]}
		!			= i \sum_b{		w_b * b * [  <u_n(k)|u_n(k+b)> -  <u_n(k)|u_n(k)>]		}
		!
		! see Mazari, Vanderbilt PRB.56.12847 (1997), Appendix B
		!
		complex(dp),	intent(in)		:: unk(:,:,:)		
		real(dp),		intent(out)		:: A(:,:,:,:)			
		complex(dp)						:: Mxl, Mxr, Myl, Myr, one
		integer							:: n, m, Z, qi, qx, qy, qxl, qxr, qyl, qyr, found, tot, al, be
		real(dp)						:: wbx,wby, bxl(2), bxr(2), byl(2), byr(2),dmax, avg, delta
		!
		A 		= 0.0_dp
		Z 		= 4	!amount of nearest neighbours( 2 for 2D cubic unit cell)
		wbx 	= 2.0_dp / 		( real(Z,dp) * dqx**2 )
		wby		= wbx
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
						!OVERLAP TO NEAREST NEIGHBOURS
						one	= UNKoverlap(	n,		m,		qi		, 	qi		, unk	)
						Mxl	= UNKoverlap(	n,		m, 		qi		,	qxl 	, unk	) 
						Mxr	= UNKoverlap(	n,		m, 		qi		,	qxr		, unk	)
						Myl	= UNKoverlap(	n,		m, 		qi		, 	qyl		, unk	)
						Myr	= UNKoverlap(	n,		m, 		qi		, 	qyr		, unk	)

						if(		 n==m 	) then		!.and.			 abs(one-dcmplx(1.0_dp)) > acc ) then
							write(*,'(a,i2,a,f6.3,a,f6.3)') "[calcConnOnCoarse]: n=m=",n," one=",dreal(one),"+i*",dimag(one)
						else
							write(*,'(a,f6.3,a,f6.3)') "[calcConnOnCoarse]:  one=",dreal(one),"+i*",dimag(one)
						end if


						!FD SUM OVER NEAREST NEIGHBOURS
						A(1:2,n,m, qi) = A(1:2,n,m, qi) + wbx * bxl(1:2) * dimag( Mxl - one )
						A(1:2,n,m, qi) = A(1:2,n,m, qi) + wbx * bxr(1:2) * dimag( Mxr - one )
						A(1:2,n,m, qi) = A(1:2,n,m, qi) + wby * byl(1:2) * dimag( Myl - one )
						A(1:2,n,m, qi) = A(1:2,n,m, qi) + wby * byr(1:2) * dimag( Myr - one )


						!A(1:2,n,m, qi) = A(1:2,n,m, qi) - wbx * bxl(1:2) * dimag( zlog( Mxl ) )
						!A(1:2,n,m, qi) = A(1:2,n,m, qi) - wbx * bxr(1:2) * dimag( zlog( Mxr ) )
						!A(1:2,n,m, qi) = A(1:2,n,m, qi) - wby * byl(1:2) * dimag( zlog( Myl ) )
						!A(1:2,n,m, qi) = A(1:2,n,m, qi) - wby * byr(1:2) * dimag( zlog( Myr ) )
						
					end do
				end do
			end do
		end do
		!
		!
		return
	end subroutine





	complex(dp) function UNKoverlap(n, m, qi, knb, unk)
		!HELPER for calcConn
		!calculates the overlap between unk at qi and at a neigbhouring k point knb
		!	integration only over the first unit cell
		!
		integer,		intent(in)		:: n, m, qi, knb
		complex(dp),	intent(in)		:: unk(:,:,:)  !unk(	nR, nK, nWfs/nG	)
		complex(dp),	allocatable		:: f(:)
		integer							:: ri
		!
		allocate( f(nR)	)
		do ri = 1, nR
			f(ri)	= dconjg( unk(ri,n,qi) ) * unk(ri,m,knb)
		end do
		!integrate, normalize for integration over only one unit cell
		UNKoverlap = nIntegrate(nR, nRx, nRy, dx, dy, f) / real(nSC,dp)
		!
		!
		return
	end function







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