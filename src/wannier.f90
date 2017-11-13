module wannier
	!contains subroutines for working with existing wannier functions
	!	center calculations, etc. 
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, nIntegrate, eigSolver
	use sysPara
	use blochWf,		only:	calcBasis
	use polarization,	only:	calcPolWannCent
	use output,			only:	writeWannFiles

	implicit none

	private
	public :: wannMethod, calcHopping

	contains



















!public:
	subroutine wannMethod(ckW,  pWann)
		complex(dp),	intent(in)		:: ckW(:,:,:)		!	unk(	nR 	,	nWfs	, nQ	)
		real(dp),		intent(out)		:: pWann(3)
		complex(dp),	allocatable		:: wnF(:,:,:)
		real(dp),		allocatable		:: wCent(:,:), wSprd(:,:)
		integer							:: normCount
		!
		allocate(			wnF( 		nR		, 	nSC		, nWfs		)				)
		allocate(			wCent(		2		, 	nWfs				)				)
		allocate(			wSprd(		2		, 	nWfs				)				)

		!
		!Generate Wannier functions and calc polarization from wannier centers
		call genWannF(ckW, wnF)
		write(*,*)	"[wannMethod]: generated the wannier functions"
		call calcCentSpread(wnF, wCent, wSprd)
		write(*,*)	"[wannMethod]: calculated centers and spreads"
		call calcPolWannCent(wCent,pWann)
		write(*,*)	"[wannMethod]: calculated polarization via centers"
		!call interpolateBands(wnF)
		!write(*,*)	"[wannMethod]: interpolated bands"
		!
		!write results
		if( writeBin )	call writeWannFiles(wnF, wCent, wSprd)
		write(*,*)	"[wannMethod]: wrote wannier files"

		!DEBUG
		if(debugWann) then
			write(*,*)	"[wannMethod]: all done. start debugging.." 
			normCount = isNormal(wnF)
			if(normCount /= 0) then
				write(*,*)	"[wannMethod]: generated wannier functions have normalization issues"
			end if
			!
			!
			call isReal(wnF)
			write(*,*)	"[wannMethod]: debugging test done"
		end if
		!
		!
		return
	end subroutine






	subroutine calcHopping(wnf, tHopp, rHopp)
		complex(dp),	intent(in)		:: wnf(:,:,:)
		complex(dp),	intent(out)		:: tHopp(:,:,:), rHopp(:,:,:,:) !tHopp(nWfs,nWfs,nSC)
		integer							:: n,m, R
		real(dp)						:: wRw(2)
		!
		write(*,*)"[calcHopping]: hello"
		do R = 1, nSC
			do m = 1, nWfs
				do n = 1, nWfs
					tHopp(n,m,R)		= wHw(R0,R,n,m, wnF)
					write(*,'(a,i3,a,i3,a,i3,a,f10.5)')"[calcHopping]: calculated R=",R,", n=",n,", m=",m,": wHw=",dreal(tHopp(n,m,R))
					!
					call wXw(R0,R,n,m, wnF, wRw)
					rHopp(1:2,n,m,R)	= wRw(1:2)
					write(*,*)"[calcHopping]: calculated wXw"
				end do
			end do
		end do
		write(*,*)"[calcHopping]: by by"
		!
		return
	end subroutine
















!privat:
	subroutine genWannF(ckW, wnF)
		! generates wannier functions from (projected) bloch wavefunctions
		!
		complex(dp), 	intent(in)  	:: ckW(:,:,:) ! ckW(nG, nWfs, nQ)
		complex(dp), 	intent(inout) 	:: wnF(:,:,:) ! wnF( 	nR, nSC, nWfs		)	
		complex(dp),	allocatable		:: basVec(:)
		integer 						:: Ri, xi, qi, gMax
		complex(dp)						:: phase
		real(dp)						:: nQreal
		!
		
		!
		nQreal	= real(nQ,dp)
		wnF		= dcmplx(0.0_dp)
		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Ri, xi, qi, phase, basVec, gMax ) 
		allocate(	basVec(nG)		)
		!$OMP DO SCHEDULE(STATIC) COLLAPSE(1) 
		do Ri = 1, nSC
			do qi = 1 , nQ
				gMax	= nGq(qi)
				do xi = 1, nR
					call calcBasis(qi, xi, basVec)
					phase			= myExp( -1.0_dp *	dot_product( qpts(:,qi) ,  Rcell(:,Ri) ) )	
					phase			= phase 	/ ( dsqrt(real(nSC,dp)) )
					!phase			= phase		/ dsqrt( vol )
					wnF(xi,Ri,:)	= wnF(xi,Ri,:) + phase * matmul( basVec(1:gMax), ckW(1:gMax,:,qi)	)
				end do
			end do
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		!
		!
		return
	end subroutine



	subroutine calcCentSpread(wnF, xc, sprd)
		!calculates the centers of the wannier functions, by evaluating the wXw expectaion values
		!	xc(n) = <R0,n|r|R0,n>		, where R0 is the home unit cell specified in the input file
		!
		complex(dp), intent(in)  	:: wnF(:,:,:) !wnF(nRpts,nSupC,nWfs)
		real(dp)   , intent(out) 	:: xc(:,:), sprd(:,:)
		integer					 	:: n
		real(dp)					:: tmp(2)
		!
		write(*,'(a,i2,a,f10.5,a,f10.5,a)')"[calcWcent]: using R(",R0,")= (",Rcell(1,R0),", ",Rcell(2,R0),") as home unit cell"
		do n = 1, nWfs
			call wXw(R0,R0,n,n,wnF, xc(:,n))
			call wXXw(R0,R0,n,n, wnF, tmp)
			sprd(:,n) = abs( tmp(:) - xc(:,n)**2 )
			write(*,'(a,i2,a,f10.6,a,f10.6,a,e16.9)')"[calcWcent]: n=",n," center= (",xc(1,n),",",xc(2,n),"), norm2(sprd)=",norm2(sprd(:,n))
			!
		end do
		!
		return
	end subroutine



!MATRIX ELEMENTS
	subroutine wXw(R1,R2,n,m, wnF, res)
		!positional operator expectation value
		integer,		intent(in)		:: R1, R2, n, m
		complex(dp),	intent(in)		:: wnF(:,:,:)	 !wnF( nRpts, nSupC,		nWfs)
		real(dp),		intent(out)		:: res(2)		
		real(dp), 		allocatable	 	:: fx(:),fy(:)
		integer 				 		:: xi
		!
		allocate(	fx(nR) 	)
		allocate(	fy(nR)	)
		!
		do xi = 1, nR	 
			if(dimag( dconjg(wnF(xi,R1,n)) * wnF(xi,R2,m)) > acc  ) then
				!write(*,*)"[wXw]: waning wnf overlap not strictly real =",dimag( dconjg(wnF(xi,R1,n)) * wnF(xi,R2,m)	)
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
	end subroutine


	subroutine wXXw(R1,R2,n,m, wnF, res)
		!positional operator squared, expectation value
		integer,		intent(in)	:: R1, R2, n, m
		complex(dp),	intent(in)	:: wnF(:,:,:)
		real(dp),		intent(out)	:: res(2)	
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
	end subroutine


	real(dp) function wHw(R1,R2,n,m, wnF)
		integer,		intent(in)		:: R1, R2, n, m
		complex(dp),	intent(in)		:: wnF(:,:,:)
		real(dp)						:: resV, resX, resY
		real(dp),		allocatable 	:: f(:)
		integer							:: ri, xi, yi, cnt, left, right 
		!
		allocate(	f(nR)	)
		!POTENTIAL TERM
		do ri = 1, nR
			f(ri)	=  dreal(	dconjg(wnf(ri,n,R1)) * getPot(ri) * wnf(ri,m,R2) )
		end do 	
		resV = nIntegrate(nR, nRx, nRy, dx, dy, f)
		!
		!
		!X DERIVATIVE
		cnt = 1
		f	= 0.0_dp
		do xi = 1, nRx
			do yi = 1, nRy
				ri		= getRindex(xi,yi)
				left	= getRleftX(xi,yi)
				right	= getRrightX(xi,yi)
				!
				f(cnt) 	= dreal( 	dconjg(wnf(ri,n,R1)) * ( wnf(right,m,R2) - 2.0_dp * wnf(ri,m,R2) + wnf(left,m,R2) )   ) / dx**2
				cnt 	= cnt + 1
			end do
		end do
		resX = nIntegrate(nR, nRx, nRy, dx, dy, f)
		!
		!
		!Y DERIVATIVE
		cnt = 1
		f	= 0.0_dp
		do xi = 1, nRx
			do yi = 1, nRy
				ri		= getRindex(xi,yi)
				left	= getRleftY(xi,yi)
				right	= getRrightY(xi,yi)
				!
				f(cnt) 	=  dreal(	dconjg(wnf(ri,n,R1)) * ( wnf(right,m,R2) - 2.0_dp * wnf(ri,m,R2) + wnf(left,m,R2) )		) / dy**2 
				cnt 	= cnt + 1
			end do
		end do
		resY = nIntegrate(nR, nRx, nRy, dx, dy, f)
		!
		!
		!SUM RESULTS
		wHw= resV - 0.5_dp * (resX + resY )
		!
		!
		return
	end function




!DEBUGGING:
	integer function isNormal(wnF)
		!checks if Wannier functions fullfill
		!	<Rn|R'm> = \delta(R,R') \delta(n,m)
		!	CURRENTLY ONLY CHECKING <1n|1m> = \delta(n,m)
		complex(dp),	intent(in)		:: wnF(:,:,:) !wnF( 	nR, nSC, nWfs		)	
		real(dp),	allocatable			:: f(:)
		real(dp)						:: oLap,  avg, dmax
		integer							:: n, m, ri, sc,sc1, sc2, tot
		!
		
		
		!
		!!!!$OMP PARALLEL DEFAULT(SHARED)	PRIVATE(n, m, ri, sc, sc1, sc2, rLog, iLog, oLap, f) &
		!!!!$OMP& REDUCTION(+:avg,tot, diffSC, isNormal) REDUCTION(max:dmax) 
		isNormal 	= 0
		sc			= 1
		dmax		= 0.0_dp
		avg			= 0.0_dp
		tot			= 0
		allocate(	f(nR)	)
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC) 
		do sc1 = 1, nSC
			do sc2 = 1, nSC
				do n = 1, nWfs
					do m = 1, nWfs
						!INTEGRATE OVERLAPS
						do ri = 1, nR
							f(ri) 	= abs( dconjg(wnF(ri,sc1,n)) * wnF(ri,sc2,n) )
						end do
						oLap 	= nIntegrate(nR, nRx, nRy, dx, dy, f)
						!APPLY CONDITION
						if( sc1 == sc2 .and. n==m ) then
							oLap = abs(oLap - 1.0_dp)
						end if
						!TEST CONDITION
						if(  abs(oLap) > acc) then
							write(*,'(a,i2,a,i2,a,i2,a,i2,a,f9.4)')	"[isNormal]: sc1=",sc1,", sc2=",sc2,", n=",n,", m=",m,": oLap=",oLap
							isNormal	= isNormal + 1
							avg	= avg + oLap
							if( oLap > dmax) then
								dmax = oLap
							end if 
						end if
						tot = tot +1
					end do
				end do
			end do
		end do
		!!!!$OMP END DO
		!!!!$OMP END PARALLEL
		avg	= avg / real(isNormal,dp)
		write(*,'(a,i5,a,i8,a,f16.7,a,f16.7)')	"[isNormal]: ",isNormal," of ",tot, &
													" are not properly normalized. dmax=",dmax," avg diff=",avg 
		!											
		!
		return
	end function


	subroutine isReal(wnF)
		complex(dp),	intent(in)		:: wnF(:,:,:)
		integer							:: ri, R, n, count, tot
		real(dp)						:: avg
		!
		count	= 0
		avg		= 0.0_dp
		tot 	= 0
		!
		do n = 1, nWfs
			do R = 1, nSC
				do ri = 1, nR
					if(	abs(dimag(wnF(ri,R,n))) > acc  ) then
						count 	= count + 1 
						avg 	= avg + abs(dimag(wnF(ri,R,n)))
					end if
					tot	= tot + 1
				end do
			end do
		end do
		!
		write(*,'(a,i7,a,i8,a,f16.12)')	"[isReal]: ",count," of ",tot,&
										" tested Wannier function grid points had non zero imaginary part; avgerage diff=",avg
		!
		!
		return
	end subroutine







end module wannier