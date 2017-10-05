module berry
	!module contains methods related to the Berry phase theory of modern polarization
	!	i.e. calculation of connection, velocities, curvatures and polarization
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, myLeviCivita, nIntegrate
	use sysPara
	use	wannier,		only:	calcHopping
	use gaugeTrafo,		only:	DoWannInterpol
	use polarization,	only:	calcPolViaA
	use semiClassics,	only:	calcFirstOrdP
	use output,			only:	writeConnCurv
	implicit none

	private
	public	::	berryMethod
	contains




!public
	subroutine berryMethod(unkW, wnf, pBerry, pNiu)
		complex(dp),	intent(in)		:: unkW(:,:,:), wnf(:,:,:)
		real(dp),		intent(out)		:: pBerry(2), pNiu(3)
		real(dp),		allocatable		:: EnH(:,:)
		complex(dp),	allocatable		:: AconnH(:,:,:,:), FcurvH(:,:,:,:), veloH(:,:,:,:),  tHopp(:,:,:), rHopp(:,:,:,:)
		!
		allocate(			EnH(		nWfs	, nK							)			)
		allocate(			AconnH(		3		, nWfs,nWfs			,	nK		)			)
		allocate(			FcurvH(		3		, nWfs, nWfs		,	nK		)			)
		allocate(			veloH(		3		, nWfs,nWfs			,	nK		)			)
		allocate(			tHopp(		nWfs	, 	nWfs	,	nSc				)			)
		allocate(			rHopp(	2	,	nWfs, 	nWfs, 	nSC					)			)						
		!
		write(*,*)	"[berrryMethod]: hello from Berry"
		!GET CONNECTION, CURVATURE & VELOCITY


		call calcHopping(wnf, tHopp, rHopp)
		
		!call calcRhopp(unkW, rHopp)

		!INTERPOLATE CONN,CURV, VELO
		call DoWannInterpol(rHopp, tHopp, EnH, AconnH, FcurvH, veloH)

		write(*,*)	"[berrryMethod]: interpolation done"
		!INTEGRATE CONNECTION
		call calcPolViaA(AconnH,pBerry)
		write(*,*)	"[berrryMethod]: calculated zero order pol"
		




		if(doNiu) then
			write(*,*)	"[berrryMethod]: now calc first order pol"
			call calcFirstOrdP(FcurvH, AconnH, veloH, EnH, pNiu)

		end if

		if(doPei) then

		end if

		write(*,*)	"[berrryMethod]: all done"

		!OUTPUT
		!call writeConnCurv(AconnH, FcurvH)
		!
		!
		return
	end subroutine




	subroutine calcRhopp(unkW, rHopp)
		complex(dp),	intent(in)		:: unkW(:,:,:)
		complex(dp),	intent(out)		:: rHopp(:,:,:,:)
		real(dp),		allocatable		:: AWCoarse(:,:,:,:)
		integer							:: qi, R
		complex(dp)						:: phase
		!
		allocate(	AWcoarse(2,nWfs,nWfs,nQ)	)
		!
		call calcConnOnCoarse(unkW, AWcoarse)
		write(*,*)	"[calcRhopp]: calculated the connection on coarse q mesh"
		do R = 1, nSC
			do qi = 1, nQ
				phase			= myExp( -1.0_dp * dot_product(qpts(:,qi),Rcell(:,R))		)
				rHopp(1,:,:,R)	= rHopp(1,:,:,R) + phase * AWcoarse(1,:,:,qi) !/ dsqrt(real(nQ,dp) )
				rHopp(2,:,:,R)	= rHopp(2,:,:,R) + phase * AWcoarse(2,:,:,qi) !/ dsqrt( real(nQ,dp) )
			end do
		end do
		!
		deallocate(	AWcoarse	)
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













end module berry






	!subroutine berryMethod(unk, En, pBerry, pNiu)
	!	complex(dp),	intent(in)		:: unk(:,:,:)
	!	real(dp),		intent(in)		:: En(:,:)
	!	real(dp),		intent(out)		:: pBerry(2), pNiu(3)
	!	complex(dp),	allocatable		:: velo(:,:,:,:)
	!	real(dp),		allocatable		:: Aw(:,:,:,:), FW(:,:,:)
	!	!
	!	allocate(			Aw(			3		, nWfs,nWfs			,	nQ		)			)
	!	allocate(			velo(		3		, nWfs,nWfs			,	nQ		)			)
	!	allocate(			FW(			3		, nWfs				,	nQ		)			)
!
!	!	call calcConnOnCoarse(unk, Aw)
!	!	call calcPolViaA(Aw,pBerry)
!	!	
!	!	!ToDo: SWITCH FOR USING NIU
!	!	if(doNiu) then
!	!		!Get zero order velocities and curvatures
!	!		!call calcVeloMat(unk, VeloBwf, velo)
!	!		call FDvelocities(unk, velo)
!	!		call calcCurv(En, velo, Fw)
!	!		!use them for calc of first order pol 
!	!		call calcFirstOrdP(Fw, Aw, velo, En, pNiu)
!	!	end if
!
!
!	!	!OUTPUT
!	!	call writeConnCurv(Aw, Fw)
!	!	!
!	!	!
!	!	return
	!end subroutine









	!subroutine calcVeloMat(unk, veloBwf, Velo)
	!	!calculates matrix elements of the velocity operator
	!	!	velocity operator is analytically applied to the plane wave basis 
	!	!	and then weighted by the basCoeff obtained from the solver and stored in veloBwf
	!	!	
	!	complex(dp),		intent(in)		:: unk(:,:,:), veloBwf(:,:,:,:)	!	unk(nR, nK, nWfs) , veloBwf(2,nR,nG,nQ)
	!	complex(dp),		intent(out)		:: Velo(:,:,:,:)   !Velo(3,	nK,	nWfs	, nwFs)		
	!	integer								:: qi, n,m, ri
	!	complex(dp),		allocatable		:: fx(:), fy(:)
	!	complex(dp)							:: BWFphase
	!	!
	!	!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(m, n, qi, ri, BWFphase, fx, fy)
	!	allocate(	fx(nR)	)
	!	allocate(	fy(nR)	)
	!	!
	!	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!	do m = 1, nWfs
	!		do n = 1, nWfs
	!			do qi = 1, nQ
	!				!FILL INTEGRATION ARRAY
	!				do ri = 1, nWfs
	!					BWFphase= myExP( 	dot_product( qpts(:,qi), rpts(:,ri) ))
	!					fx(ri)	= dconjg(	BWFphase * unk(ri,n,qi)	)	* veloBwf(1,ri, m, qi )
	!					fy(ri)	= dconjg(	BWFphase * unk(ri,n,qi)	)	* veloBwf(2,ri, m, qi )
	!				end do
	!				!INTEGRATE
	!				Velo(1,n,m,qi)	= nIntegrate(nR, nRx, nRy, dx, dy, fx)
	!				Velo(2,n,m,qi)	= nIntegrate(nR, nRx, nRy, dx, dy, fy)
	!				Velo(3,n,m,qi) 	= dcmplx(0.0_dp)
	!			end do
	!		end do
	!	end do
	!	!$OMP END DO
	!	!$OMP END PARALLEL	
	!	!
	!	return
	!end subroutine









