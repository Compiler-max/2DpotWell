module gaugeTrafo
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP, myExp, myLeviCivita, nIntegrate, eigSolver, rotMat, myCommutat
	use sysPara
	use berry,			only:	
	use wannier,		only:	genUnkW

	!use 
	implicit none

	private
	public ::					DoGaugeTrafo

	contains

!public
	subroutine DoGaugeTrafo(unkW, tHopp, EnH, AconnH, FcurvH, veloH)
		complex(dp),	intent(in)		:: unkW(:,:,:), tHopp(:,:,:)
		real(dp),		intent(out)		:: EnH(:,:)
		complex(dp),	intent(out)		:: AconnH(:,:,:,:), FcurvH(:,:,:,:)
		complex(dp),	intent(out)		:: veloH(:,:,:,:)
		complex(dp),	allocatable		:: rHopp(:,:,:,:), U(:,:), HW(:,:), HaW(:,:,:), AW(:,:,:), FW(:,:,:,:)
		integer							:: ki
		!
		
		allocate(	rHopp(	2	,	nWfs, 	nWfs, 	nSC		)		)
		allocate(	U(		nWfs, 	nWfs					)		)
		allocate(	HW(		nWfs, 	nWfs					)		)
		allocate(	HaW(	2	,	nWfs, 	nWfs			)		)
		allocate(	AW(		2	,	nWfs, 	nWfs			)		)
		allocate(	AW(		2	,	nWfs, 	nWfs			)		)
		allocate(	Fw(		2	,	2	,	nWfs,	nWfs	)	)

		call calcRhopp(unkW, rHopp)


		do ki = 1, nK
			call interpolateMat(tHopp, rHopp, HW, HaW, AW, FW)
			call gaugeBack(Hw, HaW, AW, FW, EnH(:,ki), U, AconnH(:,:,:,ki), FcurvH(:,:,:,ki), veloH(:,:,:,ki))	
		end do	



		return
	end subroutine



	subroutine gaugeBack(Hw, HaW, AW, FW, EnH, U, AconnH, FcurvH, veloH)
		complex(dp),	intent(in)		:: Hw(:,:)
		complex(dp),	intent(inout)	:: HaW(:,:,:), AW(:,:,:), FW(:,:,:,:)
		real(dp),		intent(out)		:: EnH(:)
		complex(dp),	intent(out)		:: U(:,:), AconnH(:,:,:), FcurvH(:,:,:), veloH(:,:,:)
		complex(dp),	allocatable		:: DH(:,:,:)
		integer							:: ki
		!
		allocate(	DH(2,nWfs,nWfs)	)


		!1. CALC 
		U	= HW
		!GET U MAT & ENERGIES
		call eigSolver(U, EnH)
		!ROTATE WITH u
		call calcBarMat(U, HaW, AW, FW)
		!CONNECTION
		call calcA(EnH, AW, HaW, AconnH, DH)
		!VELOCITIES
		call calcVelo(EnH, AW, HaW, veloH)
		!CURVATURE
		call calcCurv(FW, DH, AW, FcurvH)

		!
		!		
		deallocate(	DH	)
		return
	end subroutine


!gaugeBack HELPERS	
	subroutine calcCurv(FW, DH, AW, FcurvH)
		complex(dp),	intent(in)		:: FW(:,:,:,:), DH(:,:,:), AW(:,:,:)
		complex(dp),	intent(out)		:: FcurvH(:,:,:)
		complex(dp),	allocatable		:: FH(:,:,:,:)
		!
		!ToDO
		FcurvH	= dcmplx(0.0_dp)

		return
	end subroutine


	subroutine calcBarMat(U, HaW, AW, FW)
		!Helper for gaugeBack
		!	 for any quantity O: \bar(O) = U^dag O U
		complex(dp),	intent(in)		:: U(:,:)
		complex(dp),	intent(inout)	:: HaW(:,:,:), Aw(:,:,:), FW(:,:,:,:)
		complex(dp),	allocatable		:: Uc(:,:)
		integer							:: a,b
		!
		allocate(	Uc( size(U,2), size(U,1) )		)
		!	
		Uc	= dconjg( transpose(U)	)
		!
		do a = 1, 2
			HaW(a,:,:) 	= matmul(	HaW(a,:,:)	, 	U				)
			HaW(a,:,:)	= matmul(	Uc			, 	HaW(a,:,:)		)
			!
			AW(a,:,:) 	= matmul(	AW(a,:,:) 	,	U				)
			AW(a,:,:) 	= matmul(	Uc			,	AW(a,:,:)		)
			!
			do b = 1, 2
				FW(a,b,:,:) 	= matmul(	FW(a,b,:,:) 	,	U				)
				FW(a,b,:,:) 	= matmul(	Uc				,	FW(a,b,:,:)		)
			end do
		end do
		!
		!
		deallocate( Uc	)
		return
	end subroutine


	subroutine calcA(EnH, AW, HaW, AconnH, DH)
		! Helper for gaugeBack
		!	A^(H) = \bar{A}^(H) + i D^(H)
		real(dp),		intent(in)		:: EnH(:)
		complex(dp),	intent(in)		:: Aw(:,:,:), HaW(:,:,:)
		complex(dp),	intent(out)		:: AconnH(:,:,:), DH(:,:,:)
		integer							:: m, n 
		!
		!SET UP D MATRIX
		do m = 1, nWfs
			do n = 1, nWfs
				if( n /= m) then
					DH(:,n,m)	= HaW(:,n,m) / ( EnH(m) - EnH(n) + machineP	)
				else
					DH(:,n,m)	= dcmplx(0.0_dp)
				end if
			end do
		end do
		!
		!CALC CONNECTION
		AconnH(:,:,:) = AW(:,:,:) + i_dp * DH(:,:,:)
		!
		!
		return
	end subroutine


	subroutine calcVelo(EnH, AW, HaW, veloH)
		!Helper for gaugeBack
		!	v_nm	= \bar{Ha}_nm - i_dp * (Em-En) * \bar{A}_nm
		real(dp),		intent(in)		:: EnH(:)
		complex(dp),	intent(in)		:: AW(:,:,:), HaW(:,:,:)
		complex(dp),	intent(out)		:: veloH(:,:,:)
		integer							:: n, m
		!
		do m = 1, nWfs
			do n = 1, nWfs
				veloH(:,n,m)	= HaW(:,n,m) - i_dp * dcmplx(	EnH(m) - EnH(n)		) * AW(:,n,m)
			end do
		end do
		!
		return
	end subroutine















!privat
	subroutine interpolateMat(tHopp, rHopp, HW, HaW, AW, FW)
		complex(dp),	intent(in)		:: tHopp(:,:,:), rHopp(:,:,:,:)
		complex(dp),	intent(out)		:: HW(:,:), HaW(:,:,:), AW(:,:,:), FW(:,:,:,:)
		integer							:: ki, R, a, b
		complex(dp)						:: phase
		!
		do R = 1, nSC
			phase			= myExp(	dot_product(kpts(:,ki),Rcell(:,R))	)
			!HAM
			HW(:,:)		= HW(:,:) 		+ phase		 					* tHopp(:,:,R)
			!HAM DERIVATIVE
			HaW(1,:,:)	= HaW(1,:,:)	+ phase * i_dp *  Rcell(1,R) 	* tHopp(:,:,R)
			HaW(2,:,:)	= HaW(2,:,:)	+ phase * i_dp *  Rcell(2,R) 	* tHopp(:,:,R)
			!CONNECTION
			AW(1,:,:)	= AW(1,:,:) 	+ phase 						* rHopp(1,:,:,R)
			AW(2,:,:)	= AW(2,:,:) 	+ phase 						* rHopp(2,:,:,R)
			!CURVATURE
			do a = 1, 2
				do b = 1, 2
					FW(a,b,:,:) = FW(a,b,:,:) + phase * i_dp * Rcell(a,R) * rHopp(b,:,:,R)
					FW(a,b,:,:) = FW(a,b,:,:) - phase * i_dp * Rcell(b,R) * rHopp(a,:,:,R)
				end do
			end do
		end do
		!
		!
		return
	end subroutine

!interpolateMat Helpers:
	subroutine calcRhopp(unkW, rHopp)
		complex(dp),	intent(in)		:: unkW(:,:,:)
		complex(dp),	intent(out)		:: rHopp(:,:,:,:)
		real(dp),		allocatable		:: AWCoarse(:,:,:,:)
		integer							:: qi, R
		complex(dp)						:: phase
		!
		allocate(	AWcoarse(2,nR,nWfs,nQ)	)
		!
		call calcConnOnCoarse(unkW, AWcoarse)
		do R = 1, nSC
			do qi = 1, nQ
				phase			= myExp( -1.0_dp * dot_product(qpts(:,qi),Rcell(:,R))		)
				rHopp(1,:,:,R)	= rHopp(1,:,:,R) + phase * AWcoarse(1,:,:,qi) / nQ
				rHopp(2,:,:,R)	= rHopp(2,:,:,R) + phase * AWcoarse(2,:,:,qi) / nQ
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
		complex(dp),	intent(in)		:: unk(:,:,:)		!unk(	nR, nK/nKw, nWfs/nG	)
		real(dp),		intent(out)		:: A(:,:,:,:)			!Aconn(	3,nK, nWfs)		)	
		complex(dp)						:: Mxl, Mxr, Myl, Myr, one
		integer							:: n, m, Z, qi, qx, qy, qxl, qxr, qyl, qyr, found, tot
		real(dp)						:: wbx,wby, bxl(2), bxr(2), byl(2), byr(2),dmax, avg, val 
		!
		A 		= 0.0_dp
		Z 		= 4	!amount of nearest neighbours( 2 for 2D cubic unit cell)
		wbx 	= 3.0_dp / 		( real(Z,dp) * dqx**2 )
		wby 	= 3.0_dp /		( real(Z,dp) * dqy**2 )
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
		!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC),  &
		!$OMP& DEFAULT(SHARED), PRIVATE(m, n, qx, qxl, qxr, qy, qyl, qyr, qi, one, Mxl, Mxr, Myl, Myr, val),&
		!$OMP& REDUCTION(+:found,avg,tot), REDUCTION(max:dmax)
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
						!
						!OVERLAP TO NEAREST NEIGHBOURS
						one = UNKoverlap(	n,		m, 		qi		, 		qi					, unk	)
						Mxl	= UNKoverlap(	n,		m, 		qi		, getKindex( qxl, qy ) 		, unk	) 
						Mxr	= UNKoverlap(	n,		m, 		qi		, getKindex( qxr, qy )		, unk	)
						Myl	= UNKoverlap(	n,		m, 		qi		, getKindex( qx ,qyl )		, unk	)
						Myr	= UNKoverlap(	n,		m, 		qi		, getKindex( qx ,qyr )		, unk	)

						
						!FD SUM OVER NEAREST NEIGHBOURS
						A(1:2,n,m, qi) = A(1:2,n,m, qi) + wbx * bxl(1:2) * dimag( Mxl - one )
						A(1:2,n,m, qi) = A(1:2,n,m, qi) + wbx * bxr(1:2) * dimag( Mxr - one )
						A(1:2,n,m, qi) = A(1:2,n,m, qi) + wby * byl(1:2) * dimag( Myl - one )
						A(1:2,n,m, qi) = A(1:2,n,m, qi) + wby * byr(1:2) * dimag( Myr - one )
						!DEBUG:
						val	= abs( abs(one) - 1.0_dp )
						if( (val > acc .and. n==m) .or. (abs(val-1.0_dp)>acc .and. n/=m)  ) then
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
		!$OMP END PARALLEL DO
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
	!integrate
	UNKoverlap = nIntegrate(nR, nRx, nRy, dx, dy, f)
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

!!CONVERT OMEGA TENSOR TO VECTOR
!		do m = 1, nWfs
!			do n = 1, nWfs
!				do ki = 1, nK
!					!EVAL CROSS PRODUCT
!					do c = 1,3
!						do b= 1,3
!							do a=1,3
!								if( myLeviCivita(a,b,c) /= 0) then
!									Fh(c,ki,n,m) = Fh(c,ki,n,m) + myLeviCivita(a,b,c) * FhTens(a,b,ki,n,m)
!								end if
!							end do 
!						end do
!					end do
!					!
!				end do
!			end do
!		end do
!
	















end module gaugeTrafo