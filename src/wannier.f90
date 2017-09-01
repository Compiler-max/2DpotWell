module wannier
	!contains subroutines for working with existing wannier functions
	!	center calculations, etc. 
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, nIntegrate
	use sysPara, 		only: 	readInp, insideAt, getKindex, getRindex, &
									dim, aX, aY,vol, nAt, atR, atPos, atPot,&
									nG, nG0, Gcut, &
									nK, nKx, nKy, nKw, nKxW, nKyW,  nWfs, nSC, nSCx, nSCy, dkx, dky, &
									nR, nRx, nRy, R0, dx, dy, &
									Gvec, atPos, atR, kpts, kptsW, rpts,Rcell, gaugeSwitch, trialOrbSw, trialOrbVAL, Zion

	implicit none

	private
	public :: isNormal, calcWcent, calcWsprd, calc0ElPol, genTBham, genUnkW, interpConnCurv

	contains


















!public:
	integer function isNormal(wnF)
		!checks if Wannier functions fullfill
		!	<Rn|R'm> = \delta(R,R') \delta(n,m)
		!	CURRENTLY ONLY CHECKING <1n|1m> = \delta(n,m)
		complex(dp),	intent(in)		:: wnF(:,:,:) !wnF( 	nR, nSC, nWfs		)	
		complex(dp),	allocatable		:: f(:)
		complex(dp)						:: oLap
		real(dp)						:: avg, dmax
		logical							:: rLog, iLog
		integer							:: n, m, ri, sc,sc1, sc2, tot, diffSC
		!
		allocate(	f(nR)	)
		isNormal 	= 0
		sc			= 1
		dmax		= 0.0_dp
		avg			= 0.0_dp
		diffSC		= 0
		!
		!
		do sc1 = 1, nSC
			do sc2 = 1, nSC
				do n = 1, nWfs
					!INTEGRATE OVERLAPS
					do ri = 1, nR
						f(ri) 	= dconjg( wnF(ri,sc1,n) ) * wnF(ri,sc2,n)
					end do
					oLap 	= nIntegrate(nR, nRx, nRy, dx, dy, f)
					!CHECK CONDITIONS
					if( sc1 == sc2) then
						rLog 	=	abs(	abs(dreal(oLap))	-		1.0_dp	) 	> acc
						iLog	=	abs(	dimag(oLap)							)	> acc
					else
						rLog 	=	abs(	abs(dreal(oLap))					) 	> acc
						iLog	=	abs(	dimag(oLap)							)	> acc
						if( rLog .or. iLog) then
							diffSC = diffSC + 1
						end if 
					end if
					if( rLog .or. iLog) then
						!write(*,'(a,i2,a,i2,a,i2,a,i2,a,f6.3,a,f6.3)')"[isNormal]: < R=",sc,", n=",n, &
						!												" |R=",sc,", m=",n," > = ", dreal(oLap),"+i*",dimag(oLap)
						isNormal	= isNormal + 1
						avg			= avg + abs(oLap) -1.0_dp
						if( abs(oLap)-1.0_dp > dmax) then
							dmax	= abs(oLap)
						end if
					end if
					!write(*,'(a,i2,a,i2,a,i2,a,i2,a,f6.3,a,f6.3)')"[isNormal]: < R=",sc,", n=",n, &
					!									" |R=",sc,", m=",n," > = ", dreal(oLap),"+i*",dimag(oLap)
					!
					!
					if(sc1 /= sc2) then
						do m = 1, nWfs
							if(m /= n) then
								!INTEGRATE OVERLAPS
								do ri = 1, nR
									f(ri) 	= dconjg( wnF(ri,sc1,n) ) * wnF(ri,sc2,m)
								end do
								oLap 	= nIntegrate(nR, nRx, nRy, dx, dy, f)
								!CHECK CONDITIONS
								if( abs(oLap) > acc ) then
									!write(*,'(a,i2,a,i2,a,i2,a,i2,a,f6.3,a,f6.3)')"[isNormal]: < R=",sc,", n=",n, &
									!										" |R=",sc,", m=",m," > = ", dreal(oLap),"+i*",dimag(oLap)	
									isNormal	= isNormal + 1
									avg			= avg + abs(oLap)
									if( abs(oLap) > dmax) then
										dmax	= abs(oLap)
									end if		
								end if
								!write(*,'(a,i2,a,i2,a,i2,a,i2,a,f6.3,a,f6.3)')"[isNormal]: < R=",sc,", n=",n, &
								!										" |R=",sc,", m=",m," > = ", dreal(oLap),"+i*",dimag(oLap)	!
							end if
							tot = tot + 1
						end do
					end if
					!
					tot = tot + 1
				end do
			end do
		end do
		!
		avg	= avg / real(isNormal,dp)
		write(*,'(a,i5,a,i8,a,f16.12,a,f16.12)')	"[isNormal]: ",isNormal," of ",tot, &
													" are not properly normalized. dmax=",dmax," avg diff=",avg 
		write(*,'(a,i8)')	"[isNormal]: found ",diffSC, "issues between different unit cells "
		!
		return
	end




	!CENTERS AND POLARIZATION:
	subroutine calcWcent(wnF, xc)
		!calculates the centers of the wannier functions, by evaluating the wXw expectaion values
		!	xc(n) = <R0,n|r|R0,n>		, where R0 is the home unit cell specified in the input file
		!
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
		!calculates the spreads of the wannier functions, by evaluating the wXXw expectaion values
		!	sprd(n) = <R0,n|r^2|R0,n>		, where R0 is the home unit cell specified in the input file
		!
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
		!calcuates 0 order polarization from the wannier centers
		!	this is only defined up to a uncertainty quantum of e*\vec{a}/V0, 
		!	where e is electron charge, \vec{a} a bravais lattice vector and V0 the volume of the unit cell
		!	uncertainty is resolved by projecting the wannier centers into the first unit cell at (0,0) 
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
		!ELECTRONIC
		do n = 1,size(wCent,2)
			cent(1) = dmod(wCent(1,n),aX) !get current center by projection into first unit cell
			cent(2)	= dmod(wCent(2,n),aY)
			!!
			!write(*,'(a,f8.5,a,f8.5,a,f8.6,a,f8.6,a)')"[calc0ElPol]: Wcent = (",wCent(1,n),", ",wCent(2,n),") modified cent = (", cent(1),", ",cent(2),")"
			pE = pE + cent				
		end do
		!IONIC
		do at = 1, nAt
			pI = pI + Zion(at) * atPos(:,at) 
		end do
		!NORMALIZE
		!pE = pE / vol
		!pI = pI / vol
		!
		!SHIFT WITH RESPECT TO CENTER OF UNIT CELL
		cent(1)	= aX * 0.5_dp
		cent(2)	= aY * 0.5_dp
		pE = (pE - cent ) / vol
		pI = (pI - cent ) / vol 
		!TOTAL
		pT = pI + pE												
		!
		return
	end


	!TIGHT BINDING MODEL
	subroutine genTBham(ki, wnF, Htb)
		!subroutine generates a tight binding Hamiltonian via
		!
		!	H_n,m = \sum{R} exp{i k.R} <0n|H(k)|Rm>
		!
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: wnF(:,:,:)	!wnF(nR	, nSC, nWfs	)
		complex(dp),	intent(out)		:: Htb(:,:,:) 	!Htb(R, nWfs,nWfs)
		integer							:: n,m, R
		complex(dp)						:: phaseR
		!
		Htb	= dcmplx(0.0_dp)
		!
		!
		do m = 1, nWfs
			do n = 1, nWfs
				!SUM OVER CELLS R
				do R = 1, nSC
					phaseR		= myExp( 	dot_product( kpts(:,ki) , Rcell(:,R) )			)	
					Htb(R,n,m)	= Htb(R,n,m) + phaseR * Hwnf(n,m,ki,R, wnF)
				end do
			end do
		end do
		!
		return
	end





	!INTERPOLATION
	subroutine genUnkW(wnF, unk)
		!generates lattice periodic functions unk from the Wannier functions wnf
		!	currently uses the intial k points again, but could also be interpolated into coarse mesh
		complex(dp),	intent(in)		:: wnF(:,:,:)			!	wnF( 	nR, nSC, nWfs	)	
		complex(dp),	intent(out)		:: unk(:,:,:)
		integer							:: n, R, ki, xi
		real(dp)						:: cellP
		!
		unk = dcmplx(0.0_dp)
		!GENERATE BLOCH LIKE FUNCTIONS
		do n = 1, nWfs
			do R = 1, nSC
				do ki = 1, nKw
					cellP = dot_product(	kptsW(:,ki) , 	Rcell(:,R)	)
					do xi = 1, nR
						unk(xi,ki,n) = unk(xi,ki,n) + myExp(cellP) * wnF(xi,R,n) 	  
					end do
				end do
			end do
		end do
		!EXTRACT LATTICE PERIODIC PART
		do n = 1, nWfs
			do ki = 1, nKw
				do xi = 1, nR
					cellP			= -1.0_dp * dot_product( kpts(:,ki) ,	rpts(:,xi)	)	
					unk(xi,ki,n)	= myExp(cellP) * unk(xi,ki,n)
				end do
			end do
		end do
		!
		!
		return
	end


	subroutine interpConnCurv(wnF, Aconn, Fcurv)
		!calculates the Berry connection from the Wannier centers
		!
		!	A_n(k)		= \sum_R exp^{i k.R} <0n|r|Rn>
		!
		!and the connection via
		!
		! 	F_a,b_n(K)	= d_a A_b - d_b A_a
		!
		!where the derivatives are evaluated analytically
		!	d_a A_b 	= d_a \sum_R exp^{i k.R} <0n|r_b|Rn>
		!				= \sum_R (i R_a) exp^{i k.R} <0n|r_b|Rn>
		!
		complex(dp),	intent(in)		:: wnF(:,:,:)					 !		wnF( 	nR		, 	nSC		, nWfs	)	
		complex(dp),	intent(out)		:: Aconn(:,:,:), Fcurv(:,:,:,:)  !		Aconn(	2	,	nK	, nWfs	)		
		integer							:: n, Ri, ki
		real(dp)						:: cellP, rExpec(2)
		!
		Aconn	= dcmplx(0.0_dp)
		Fcurv	= dcmplx(0.0_dp)
		!
		do n = 1, nWfs
			do Ri = 1, nSC
				call wXw(1,Ri,n,n,wnF, rExpec)
				do ki =1, nKw
					cellP			= dot_product( 	kptsW(:,ki)	,	Rcell(:,Ri)		)		
					Aconn(:,ki,n)	= Aconn(:,ki,n) +  myExp(cellP) * dcmplx( rExpec )
					!xy
					Fcurv(1,2,ki,n)	= Fcurv(1,2,ki,n) + myExp(cellP) * i_dp * dcmplx( Rcell(1,n) * rExpec(2) )
					Fcurv(1,2,ki,n)	= Fcurv(1,2,ki,n) - myExp(cellP) * i_dp * dcmplx( Rcell(2,n) * rExpec(1) )
					!yx
					Fcurv(2,1,ki,n)	= Fcurv(2,1,ki,n) + myExp(cellP) * i_dp * dcmplx( Rcell(2,n) * rExpec(1) )
					Fcurv(2,1,ki,n)	= Fcurv(2,1,ki,n) - myExp(cellP) * i_dp * dcmplx( Rcell(1,n) * rExpec(2) )
				end do
			end do
		end do
		!
		return
	end





















!privat:
		subroutine wXw(R1,R2,n,m, wnF, res)
		!positional operator expectation value
		integer,		intent(in)	:: R1, R2, n, m
		complex(dp),	intent(in)	:: wnF(:,:,:)	 !wnF( nRpts, nSupC,		nWfs)
		real(dp),		intent(out)	:: res(2)		
		real(dp)				 	:: norm
		real(dp), allocatable	 	:: fx(:),fy(:)
		integer 				 	:: xi
		!
		allocate(	fx(nR) 	)
		allocate(	fy(nR)	)
		!
		do xi = 1, nR	 
			if(dimag( dconjg(wnF(xi,R1,n)) * wnF(xi,R2,m)) > acc  ) then
				write(*,*)"[wXw]: waning wnf overlap not strictly real =",dimag( dconjg(wnF(xi,R1,n)) * wnF(xi,R2,m))
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


	complex(dp) function Hwnf(n,m, ki, R, wnF)
		!calculates the Hamiltonian expectation value
		!
		!	<0n|H(ki)|Rm>
		!
		integer,		intent(in)		:: n, m, ki, R
		complex(dp),	intent(in)		:: wnF(:,:,:)	!!wnF(nR	, nSC, nWfs	)
		real(dp)						:: Htmp
		complex(dp),	allocatable		:: f(:)
		integer							:: ri, at
		!
		allocate(	f(nR)	)
		!
		do ri = 1, nR
			!SET UP HAM ON REAL SPACE GRID
			Htmp	= 0.5_dp	* dot_product(kpts(:,ki),kpts(:,ki))			!hbar**2 / (2*me) * k**2 kin.energy
			do at = 1, nAt
				if( insideAt(at,rpts(:,ri))		) then
					Htmp	= Htmp + atPot(at)
				end if
			end do
			!FILL INTEGRATION ARRAY
			f(ri)	= dconjg( wnF(ri,R0,n) )	*	dcmplx( Htmp )		*	wnF(ri,R,m)
		end do
		!
		Hwnf	= nIntegrate(nR, nRx,nRy, dx,dy, f)
		!
		!
		return
	end





end module wannier