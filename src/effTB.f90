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
		complex(dp),		allocatable		:: Htmp(:,:,:), rImag(:,:,:,:), tImag(:,:,:)
		complex(dp)							:: phase
		integer								:: R, qi
		!
		allocate(	AconnQ(	2,	nWfs,	nWfs,	nQ		)			)
		allocate(	Htmp(			nWfs,	nWfs,	nQ		)			)
		allocate(	rImag(		2,	nWfs,	nwfs,	nSC		)			)
		allocate(	timag(			nWfs,	nWfs,	nSC		)			)
		rImag	= dcmplx(0.0_dp)
		tImag	= dcmplx(0.0_dp)
		!
		!SET UP K SPACE QUANTITIES
		call calcConnOnCoarse(unkQ, AconnQ)
		call calcHtmp(EnQ, Uq, Htmp)
		!FT TO REAL SPACE
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
		integer							:: n, m, Z, qi, qx, qy, qxl, qxr, qyl, qyr, found, tot
		real(dp)						:: wbx,wby, bxl(2), bxr(2), byl(2), byr(2),dmax, avg, val 
		!
		A 		= 0.0_dp
		Z 		= 4	!amount of nearest neighbours( 2 for 2D cubic unit cell)
		wbx 	= 1.0_dp / 		( real(Z,dp) * dqx**2 )
		wby 	= 1.0_dp /		( real(Z,dp) * dqy**2 )
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
		found	= 0
		tot		= 0
		dmax	= 0.0_dp
		avg		= 0.0_dp
		!
		!!!!!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC),  &
		!!!!!$OMP& DEFAULT(SHARED), PRIVATE(m, n, qx, qxl, qxr, qy, qyl, qyr, qi, one, Mxl, Mxr, Myl, Myr, val),&
		!!!!!$OMP& REDUCTION(+:found,avg,tot), REDUCTION(max:dmax)
		do m = 1, nWfs
			do n = 1, nWfs
				do qx = 1, nQx
					qxl	= getLeft( qx,nQx)
					qxr	= getRight(qx,nQx)
					!
					do qy = 1, nQy
						qyl	= getLeft(  qy,nQy)
						qyr = getRight( qy,nQy)
						qi	= getKindex(qx,qy)
						if(		norm2(qpts(:,qi)-qpts(:,qyl)) > dqx		) then
							write(*,*)"[calcConnOnCoarse]: fd does not get the neighbours"
						end if
						if(		norm2(qpts(:,qyr)-qpts(:,qi)) > dqx		) then
							write(*,*)"[calcConnOnCoarse]: fd does not get the neighbours"
						end if
						!
						!OVERLAP TO NEAREST NEIGHBOURS
						!one = UNKoverlap(	n,		m, 		qi		, 		qi					, unk	)
						Mxl	= UNKoverlap(	n,		m, 		qi		, getKindex( qxl, qy ) 		, unk	) 
						Mxr	= UNKoverlap(	n,		m, 		qi		, getKindex( qxr, qy )		, unk	)
						Myl	= UNKoverlap(	n,		m, 		qi		, getKindex( qx ,qyl )		, unk	)
						Myr	= UNKoverlap(	n,		m, 		qi		, getKindex( qx ,qyr )		, unk	)

						if(n == m) then
							one	= dcmplx(1.0_dp)
						else
							one = dcmplx(0.0_dp)
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
						!DEBUG:
						if(n == m) then
							val	= abs( abs(one) - 1.0_dp )
						else
							val = abs(one)
						end if
						
						if( val > acc ) then
							!write(*,'(a,i2,a,i7,a,f16.8,a,f16.8)')	"[calcConn]: n=",n," unk normalization problem at ki=",ki,&
							!							" one=",dreal(one),"+i*",dimag(one)
							found 	= found + 1
							avg		= avg + val
							if( val > dmax) then
								dmax = val
							end if 
						end if
						tot = tot + 1
					end do
				end do
			end do
		end do
		!!!!!$OMP END PARALLEL DO
		!
		!DEBUG
		avg	= avg / real(found,dp)
		write(*,'(a,i6,a,i6,a,f16.12,a,f16.12)')	"[calcConnOnCoarse]: ",found," of ",tot,&
										" checked unk functions had normalization issues;  max delta=",dmax,&
										" avg diff=",avg
		!
		return
	end subroutine


	integer function getLeft(i,N)
		!HELPER for calcConn
		!gets left (lower) neighbour, using the periodicity at boundary
		!
		integer,	intent(in)	:: i,N
		if(i==1) then
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
		if(i==N) then
			getRight = 1
		else
			getRight = i+1
		end if
		!
		return
	end function



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
		if( dimag(UNKoverlap) > acc ) then
			write(*,*)"[UNKoverlap]: none zero imaginary part:,",dimag(UNKoverlap)
		end if
		!
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