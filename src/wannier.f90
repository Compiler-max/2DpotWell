module wannier
	!contains subroutines for working with existing wannier functions
	!	center calculations, etc. 
	use mathematics,	only:	dp, PI_dp, i_dp, myExp, nIntegrate
	use sysPara, 		only: 	readInp, insideAt, getKindex, getRindex, &
									dim, aX, aY,vol, nAt, atR, atPos, atPot,&
									nG, nG0, Gcut, nK, nKx, nKy,  nWfs, nSC, nSCx, nSCy, nR, nRx, nRy, R0, dx, dy, dkx, dky, &
									Gvec, atPos, atR, kpts, rpts,Rcell, gaugeSwitch, trialOrbSw, trialOrbVAL, Zion

	implicit none

	private
	public :: isNormal, calcWcent, calcWsprd, calc0ElPol, genUnkW, calcConn, calcPolViaA 

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
		do n = 1,size(wCent,2)
			cent(1) = dmod(wCent(1,n),aX) !get current center by projection into first unit cell
			cent(2)	= dmod(wCent(2,n),aY)
			!
			write(*,'(a,f8.5,a,f8.5,a,f8.6,a,f8.6,a)')"[calc0ElPol]: Wcent = (",wCent(1,n),", ",wCent(2,n),") modified cent = (", cent(1),", ",cent(2),")"
			pE = pE + cent					!electron contribution
		end do
		do at = 1, nAt
			pI = pI + Zion(at) * atPos(:,at) 	!Ion contribution  
		end do
		!
		pE = -pE / vol
		pI = pI / vol
		pT =pI + pE				!Total											
		!
		return
	end



	!INTERPOLATION
	subroutine genUnkW(wnF,kpts, unk)
		!generates lattice periodic functions unk from the Wannier functions wnf
		!	currently uses the intial k points again, but could also be interpolated into coarse mesh
		complex(dp),	intent(in)		:: wnF(:,:,:)			!	wnF( 	nR, nSC, nWfs	)	
		real(dp),		intent(in)		:: kpts(:,:)
		complex(dp),	intent(out)		:: unk(:,:,:)
		integer							:: n, R, ki, xi
		real(dp)						:: cellP
		!
		unk = dcmplx(0.0_dp)
		!GENERATE BLOCH LIKE FUNCTIONS
		do n = 1, nWfs
			do R = 1, nSC
				do ki = 1, size(unk,2)
					cellP = dot_product(	kpts(:,ki) , 	Rcell(:,R)	)
					do xi = 1, nR
						unk(xi,ki,n) = unk(xi,ki,n) + myExp(cellP) * wnF(xi,R,n) 	  
					end do
				end do
			end do
		end do
		!EXTRACT LATTICE PERIODIC PART
		do n = 1, nWfs
			do ki = 1, nK
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



	subroutine calcConn(unk, A)
		!finite difference on lattice periodic unk to calculate the Berry connection A
		!	A_n(k) = <u_n(k)|i \nabla_k|u_n(k)>
		!
		! see Mazari, Vanderbilt PRB.56.12847 (1997)
		!
		complex(dp),	intent(in)		:: unk(:,:,:)		!unk(	nR, nK, nWfs/nG	)
		complex(dp),	intent(out)		:: A(:,:,:)			!Aconn(	2,nK, nWfs)		)	
		complex(dp)						:: Mxl, Mxr, Myl, Myr, M, one
		integer							:: n, Z, ki, kx, ky, kxl, kxr, kyl, kyr
		real(dp)						:: thres, wbx,wby, bxl(2), bxr(2), byl(2), byr(2) !for nearest neighbours, assuming cubic mesh
		!
		thres	= 1e-3_dp
		A 		= dcmplx(0.0_dp)
		Z 		= 4	!amount of nearest neighbours( 2 for 2D cubic unit cell)
		wbx 	= 3.0_dp / 		( real(Z,dp) * dkx**2 )
		wby 	= 3.0_dp /		( real(Z,dp) * dky**2 )
		!b vector two nearest X neighbours:
		bxl(1) 	= -dkx
		bxl(2)	= 0.0_dp
		bxr(1) 	= +dkx
		bxr(2)	= 0.0_dp
		!b vector two nearest Y neighbours:
		byl(1) 	= 0.0_dp
		byl(2)	= -dky
		byr(1) 	= 0.0_dp
		byr(2)	= +dky


		do n = 1, nWfs
			do kx = 1, nKx
				kxl	= getLeft(kx,nKx)
				kxr	= getRight(kx,nKx)
				!
				do ky = 1, nKy
					kyl	= getLeft(ky,nKy)
					kyr = getRight(ky,nKy)
					ki	= getKindex(kx,ky)
					!
					!OVERLAP TO NEAREST NEIGHBOURS
					one = overlap(	n, 		ki		, 		ki					, unk	)
					Mxl	= overlap(	n, 		ki		, getKindex( kxl, ky ) 		, unk	) 
					Mxr	= overlap(	n, 		ki		, getKindex( kxr, ky )		, unk	)
					Myl	= overlap(	n, 		ki		, getKindex( kx ,kyl )		, unk	)
					Myr	= overlap(	n, 		ki		, getKindex( kx ,kyr )		, unk	)
					!
					!write(*,'(a,f15.12,a,f15.12)')"[calcConn]: Mxl=",dreal(Mxl),"+i*",dimag(Mxl)
					!FD SUM OVER NEAREST NEIGHBOURS
					A(:,ki,n) = A(:,ki,n) + wbx * bxl(:) * i_dp * ( Mxl - one )
					A(:,ki,n) = A(:,ki,n) + wbx * bxr(:) * i_dp * ( Mxr - one )
					A(:,ki,n) = A(:,ki,n) + wby * byl(:) * i_dp * ( Myl - one )
					A(:,ki,n) = A(:,ki,n) + wby * byr(:) * i_dp * ( Myr - one )
					!FD SUM OVER NEAREST NEIGHBOURS
					!A(1,ki,n) = A(1,ki,n) + wbx * bxl(1) * i_dp * ( Mxl - one )
					!A(2,ki,n) = A(2,ki,n) + wbx * bxl(2) * i_dp * ( Mxl - one )
					!A(1,ki,n) = A(1,ki,n) + wbx * bxr(1) * i_dp * ( Mxr - one )
					!A(2,ki,n) = A(2,ki,n) + wbx * bxr(2) * i_dp * ( Mxr - one )
					!A(1,ki,n) = A(1,ki,n) + wby * byl(1) * i_dp * ( Myl - one )
					!A(2,ki,n) = A(2,ki,n) + wby * byl(2) * i_dp * ( Myl - one )
					!A(1,ki,n) = A(1,ki,n) + wby * byr(1) * i_dp * ( Myr - one )
					!A(2,ki,n) = A(2,ki,n) + wby * byr(2) * i_dp * ( Myr - one )

					!
					if(abs( abs(one) - 1.0_dp ) > thres ) then
						write(*,'(a,i2,a,i7,a,f16.8,a,f16.8)')	"[calcConn]: n=",n," unk normalization problem at ki=",ki,&
													" one=",dreal(one),"+i*",dimag(one)
					end if
				end do
			end do
		end do
		!
		!
		return
	end


	complex(dp) function overlap(n, ki, knb, unk)
		!HELPER for calcConn
		!calculates the overlap between unk at ki and at a neigbhouring k point knb
		!	integration only over the first unit cell
		!
		integer,		intent(in)		:: n, ki, knb
		complex(dp),	intent(in)		:: unk(:,:,:)  !unk(	nR, nK, nWfs/nG	)
		complex(dp),	allocatable		:: f(:)
		integer							:: xi,yi,ri,rloc, nRx1, nRy1, nR1
		!
		!Set integration range to first unit cell
		nRx1 	= int(		real(nRx,dp) / real(nSCx,dp)		)
		nRy1 	= int(		real(nRy,dp) / real(nSCy,dp)		)
		nR1		= nRx1 * nRy1 
		allocate(	f(nR1)	)
		!
		!fill integration array
		f 		= dcmplx(0.0_dp)
		do yi = 1, nRy1
			do xi = 1, nRx1
				ri		= getRindex(xi,yi)			!overall index, to get correct position from unk
				rloc 	= (yi-1) * nRx1 + xi		!for mapping to f array
				f(rloc)	= dconjg( unk(ri,ki,n) ) * unk(ri,knb,n)
				!write(*,'(a,f10.6,a,f10.6)')	"[overlap] f=",dreal(f(rloc)),"+i*",dimag(f(rloc))
			end do
		end do
		!
		!integrate
		overlap = nIntegrate(nR1, nRx1, nRy1, dx, dy, f	)
		
		!write(*,'(a,f10.6,a,f10.6)')"[overlap]=",dreal(overlap),"+i*",dimag(overlap)
		!
		!
		return
	end


	subroutine calcPolViaA(A, pElA)
		!calculates the polarization by integrating connection over the brillouin zone
		! r_n 	= <0n|r|0n> 
		!		=V/(2pi)**2 \integrate_BZ <unk|i \nabla_k|unk>
		!		=V/(2pi)**2 \integrate_BZ A(k)
		complex(dp),	intent(in)		:: A(:,:,:)			!A(2,	 nWfs, nK	)	
		real(dp),		intent(out)		:: pElA(2)
		complex(dp)						:: val(2)
		real(dp)						:: thres
		integer							:: n, ki
		!
		val		= dcmplx(0.0_dp)
		pElA	= 0.0_dp
		thres	= 1e-10_dp
		!
		!SUM OVER K SPACE AND OVER STATES
		do n 	= 1, nWfs
			do ki = 1, nK
				val(1) = val(1) + A(1,ki,n)
				val(2) = val(2) + A(2,ki,n)
			end do
			!val(1)	= nIntegrate(nK, nKx, nKy, dkx, dky, A(1,n,:)	)
			!val(2)	= nIntegrate(nK, nKx, nKy, dkx, dky, A(2,n,:)	)
		end do
		!NORMALIZE
		val(1) = val(1) / real(nK,dp)
		val(2) = val(2) / real(nK,dp)
		!
		!HARVEST
		pElA(1) = dreal( val(1)	)
		pElA(2) = dreal( val(2)	)
		!
		!
		!DEBUGGING
		if(		dimag( val(1) ) > thres 	.or. 	dimag( val(2) ) > thres		) then
			write(*,*)"[calcPolViaA]: non zero imaginary part of polarization"
		end if
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
	end

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
	end




end module wannier