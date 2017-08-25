module wannier
	!contains subroutines for working with existing wannier functions
	!	center calculations, etc. 
	use mathematics,	only:	dp, PI_dp,i_dp, myExp, nIntegrate
	use sysPara, 		only: 	readInp, insideAt, &
									dim, aX, aY,vol, nAt, atR, atPos, atPot,&
									nG, nG0, Gcut, nK, nWfs, nSC,nR, nRx, nRy, R0, dx, dy, dkx, dky, &
									Gvec, atPos, atR, kpts, rpts,Rcell, gaugeSwitch, trialOrbSw, trialOrbVAL, Zion

	implicit none

	private
	public :: isNormal

	contains


















!public:
	logical function isNormal(wnF)
		!checks if Wannier functions fullfill
		!	<Rn|R'm> = \delta(R,R') \delta(n,m)
		!	CURRENTLY ONLY CHECKING <1n|1m> = \delta(n,m)
		complex(dp),	intent(in)		:: wnF(:,:,:) !wnF( 	nR, nSC, nWfs		)	
		complex(dp),	allocatable		:: f(:)
		complex(dp)						:: oLap
		real(dp)						:: thres
		logical							:: rLog, iLog
		integer							:: n, m, ri, sc
		!
		allocate(	f(nR)	)
		isNormal 	= .true.
		thres		= 1e-2_dp
		sc			= 1

		!
		!
		n 			= 1
		do while ( n<= nWfs .and. isNormal )
			!INTEGRATE OVERLAPS
			do ri = 1, nR
				f(ri) 	= dconjg( wnF(ri,sc,n) ) * wnF(ri,sc,n)
			end do
			oLap 	= nIntegrate(nR, nRx, nRy, dx, dy, f)
			!CHECK CONDITIONS
			rLog 	=	abs(			abs(dreal(oLap))		-		1.0_dp	) 	> thres
			iLog	=	abs(			dimag(oLap)								)	> thres
			if( rLog .or. iLog) then
				write(*,'(a,i2,a,i2,a,i2,a,i2,a,f6.3,a,f6.3)')"[isNormal]: < R=",sc,", n=",n, &
																" |R=",sc,", m=",n," > = ", dreal(oLap),"+i*",dimag(oLap)
				isNormal = .false.
			end if
			!write(*,'(a,i2,a,i2,a,i2,a,i2,a,f6.3,a,f6.3)')"[isNormal]: < R=",sc,", n=",n, &
			!									" |R=",sc,", m=",n," > = ", dreal(oLap),"+i*",dimag(oLap)
			!
			!
			!
			m = 1
			do while (m<= nWfs	.and. isNormal )
				if(m /= n) then
					!INTEGRATE OVERLAPS
					do ri = 1, nR
						f(ri) 	= dconjg( wnF(ri,sc,n) ) * wnF(ri,sc,m)
					end do
					oLap 	= nIntegrate(nR, nRx, nRy, dx, dy, f)
					!CHECK CONDITIONS
					if( abs(oLap) > thres ) then
						write(*,'(a,i2,a,i2,a,i2,a,i2,a,f6.3,a,f6.3)')"[isNormal]: < R=",sc,", n=",n, &
																" |R=",sc,", m=",m," > = ", dreal(oLap),"+i*",dimag(oLap)	
						isNormal = .false.						
					end if
					!write(*,'(a,i2,a,i2,a,i2,a,i2,a,f6.3,a,f6.3)')"[isNormal]: < R=",sc,", n=",n, &
					!										" |R=",sc,", m=",m," > = ", dreal(oLap),"+i*",dimag(oLap)	!
				end if
				m	= m+1
			end do
			!
			!
			!
			n 		= n + 1
		end do
		!
		return
	end




	!CENTERS AND POLARIZATION:
	subroutine calcWcent(wnF, xc)
		complex(dp), intent(in)  :: wnF(:,:,:) !wnF(nRpts,nSupC,nWfs)
		real(dp)   , intent(out) :: xc(:,:)
		integer					 :: n
		!
		write(*,'(a,i2,a,f10.5,a,f10.5,a)')"[calcWcent]: using R(",R0,")= (",Rcell(1,R0),", ",Rcell(2,R0),") as home unit cell"
		do n = 1, nWfs
			call wXw(R0,R0,n,n,wnF, xc(:,n))
			write(*,'(a,i2,a,f10.6,a,f10.6,a)')"[calcWcent]: n=",n," center= (",xc(1,n),",",xc(2,n),")"
		end do
		!
		return
	end


	subroutine calcWsprd(wnF, xc, sprd)
		complex(dp), intent(in)		:: wnF(:,:,:)
		real(dp),	 intent(in)		:: xc(:,:)
		real(dp)   , intent(out)	:: sprd(:,:)
		integer					 	:: n
		real(dp)					:: tmp(2)
		!
		do n = 1, nWfs
			call wXXw(R0,R0,n,n, wnF, tmp)
			sprd(1,n) = abs( tmP(1) - xc(1,n)**2 )
			sprd(2,n) = abs( tmP(2) - xc(2,n)**2 )
		end do
		!
		return
	end


	subroutine calc0ElPol(wCent,pE, pI, pT)
		real(dp), intent(in)	:: wCent(:,:)
		real(dp), intent(out)	:: pE(2), pI(2), pT(2)
		integer 				:: n, at
		real(dp) 				:: vol, cent(2)
		!
		vol 	= aX * aY !1D
		pT		= 0.0_dp
		pE	 	= 0.0_dp
		pI		= 0.0_dp
		!
		do n = 1,size(wCent)
			cent(1) = dmod(wCent(1,n),aX)
			cent(2)	= dmod(wCent(2,n),aY)
		
			write(*,'(a,f8.5,a,f8.5,a,f8.6,a,f8.6,a)')"[calc0ElPol]: Wcent = (",wCent(1,n),", ",wCent(2,n),") modified cent = (", cent(1),", ",cent(2),")"
			pE = pE + cent					!electron contribution
		end do
		do at = 1, nAt
			pI = pI + Zion(at) * atPos(:,at) 	!Ion contribution  
		end do
		!
		!pE = dmod(pE, aLatt)
		pE = -pE / vol
		pI = pI / vol
		pT =pI + pE				!Total											
		!
		return
	end














!privat:
		subroutine wXw(R1,R2,n,m, wnF, res)
		!positional operator expectation value
		integer,		intent(in)	:: R1, R2, n, m
		complex(dp),	intent(in)	:: wnF(:,:,:)	 !wnF( nRpts, nSupC,		nWfs)
		real(dp),		intent(out)	:: res(2)		
		real(dp)				 	:: norm, thres
		real(dp), allocatable	 	:: fx(:),fy(:)
		integer 				 	:: xi
		!
		thres = 1e-10_dp
		allocate(	fx(nR) 	)
		allocate(	fy(nR)	)
		!
		do xi = 1, nR	 
			if(dimag( dconjg(wnF(xi,R1,n)) * wnF(xi,R2,m)) > thres  ) then
				write(*,*)"[wXw]: waning wnf overlap not strictly real"
			end if
			fx(xi)	=  rpts(1,xi)		*	dreal( dconjg(wnF(xi,R1,n)) * wnF(xi,R2,m)  )
			fy(xi)	=  rpts(2,xi)		*	dreal( dconjg(wnF(xi,R1,n)) * wnF(xi,R2,m)  )
		end do
		!
		!call nIntegrate(nIntSwitch, dx, f	, wXw	)
		res(1)		= nIntegrate(nR, nRx, nRy, dx, dy, fx)
		res(2)		= nIntegrate(nR, nRx, nRy, dx, dy, fy)
		!
		!write(*,*)"[wXw]: scaling with ",1.0_dp/norm
		!wXw = wXw  /norm
		return
	end


	subroutine wXXw(R1,R2,n,m, wnF, res)
		!positional operator squared, expectation value
		integer,		intent(in)	:: R1, R2, n, m
		complex(dp),	intent(in)	:: wnF(:,:,:)
		real(dp),		intent(out)	:: res(2)	
		real(dp)				 	:: norm
		real(dp), allocatable	 	:: fx(:),fy(:)
		integer 				 	:: xi
		!
		allocate( 	fx(nR) 	)
		allocate(	fy(nR)	)
		do xi = 1, nR	 !
			fx(xi)  = rpts(1,xi)**2 		*	dreal( dconjg(wnF(xi,R1,n)) * wnF(xi,R2,m)  )
			fy(xi)	= rpts(2,xi)**2			*	dreal( dconjg(wnF(xi,R1,n)) * wnF(xi,R2,m)  )
		end do
		!
		res(1)	= nIntegrate(nR, nRx, nRy, dx, dy, fx)
		res(2)	= nIntegrate(nR, nRx, nRy, dx, dy, fy)
		!call nIntegrate(nIntSwitch, dx, fn	, norm	)
		!write(*,*)"[wXXw]: scaling with ",1.0_dp/norm
		!wXXw = wXXw  / norm
		!
		return
	end




end module wannier