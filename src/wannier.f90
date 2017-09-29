module wannier
	!contains subroutines for working with existing wannier functions
	!	center calculations, etc. 
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, nIntegrate, eigSolver
	use sysPara
	use polarization,	only:	calcPolWannCent
	use output,			only:	writeWannFiles, writeInterpBands

	implicit none

	private
	public :: wannMethod, isNormal, isReal, calcCentSpread, genTBham, genUnkW, calcWannMat

	contains


















!public:
	subroutine wannMethod(unk, pWann)
		complex(dp),	intent(in)		:: unk(:,:,:)		!	unk(	nR 	,	nWfs	, nQ	)
		real(dp),		intent(out)		:: pWann(2)
		complex(dp),	allocatable		:: wnF(:,:,:)		
		real(dp),		allocatable		:: wCent(:,:), wSprd(:,:)
		integer							:: normCount
		!
		allocate(			wnF( 		nR		, 	nSC		, nWfs		)				)
		allocate(			wCent(		2		, 	nWfs				)				)
		allocate(			wSprd(		2		, 	nWfs				)				)
		!
		!Generate Wannier functions and calc polarization from wannier centers
		call genWannF2(unk, wnF)
		write(*,*)	"[wannMethod]: generated the wannier functions"
		call calcCentSpread(wnF, wCent, wSprd)
		write(*,*)	"[wannMethod]: calculated centers and spreads"
		call calcPolWannCent(wCent,pWann)
		write(*,*)	"[wannMethod]: calculated polarization via centers"
		!call interpolateBands(wnF)
		!write(*,*)	"[wannMethod]: interpolated bands"
		!
		!write results
		call writeWannFiles(wnF, wCent, wSprd)
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
		
		
		!
		!$OMP PARALLEL DEFAULT(SHARED)	PRIVATE(n, m, ri, sc, sc1, sc2, rLog, iLog, oLap, f) &
		!$OMP& REDUCTION(+:avg,tot, diffSC, isNormal) REDUCTION(max:dmax) 
		isNormal 	= 0
		sc			= 1
		dmax		= 0.0_dp
		avg			= 0.0_dp
		diffSC		= 0
		allocate(	f(nR)	)
		!$OMP DO COLLAPSE(3) SCHEDULE(STATIC) 
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
		!$OMP END DO
		deallocate(	f )
		!$OMP END PARALLEL
		avg	= avg / real(isNormal,dp)
		write(*,'(a,i5,a,i8,a,f16.12,a,f16.12)')	"[isNormal]: ",isNormal," of ",tot, &
													" are not properly normalized. dmax=",dmax," avg diff=",avg 
		write(*,'(a,i8)')	"[isNormal]: found ",diffSC, "issues between different unit cells "
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




	!CENTERS AND POLARIZATION:
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
			sprd(1,n) = abs( tmP(1) - xc(1,n)**2 )
			sprd(2,n) = abs( tmP(2) - xc(2,n)**2 )
			write(*,'(a,i2,a,f10.6,a,f10.6,a,e16.9)')"[calcWcent]: n=",n," center= (",xc(1,n),",",xc(2,n),"), norm2(sprd)=",norm2(sprd(:,n))
			!
		end do
		!
		return
	end subroutine






	!TIGHT BINDING MODEL
	subroutine interpolateBands(wnF)
		complex(dp),	intent(in)		:: wnF(:,:,:)
		complex(dp),	allocatable		:: Htb(:,:)
		real(dp),		allocatable		:: Ew(:,:)
		integer							:: ki
		!
		allocate(	Htb(nWfs,nWfs)	)
		allocate(	Ew(nWfs, nK)	)
		!
		do ki = 1, nK
			call genTBham(ki,wnF,Htb)
			call eigSolver(Htb, Ew(:,ki))
		end do
		!
		call writeInterpBands(Ew)
		!
		!
		return
	end subroutine

	subroutine genTBham(ki, wnF, Htb)
		!subroutine generates a tight binding Hamiltonian via
		!
		!	H_n,m = \sum{R} exp{i k.R} <0n|H(k)|Rm>
		!
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: wnF(:,:,:)	!wnF(nR	, nSC, nWfs	)
		complex(dp),	intent(out)		:: Htb(:,:) 	!Htb(nWfs,nWfs)
		integer							:: n,m, R
		complex(dp)						:: phaseR, Hexp
		!
		Htb	= dcmplx(0.0_dp)
		!
		!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(m,n, R, phaseR, Hexp) ,&
		!$OMP& COLLAPSE(2)	SCHEDULE(STATIC) 
		do m = 1, nWfs
			do n = 1, nWfs
				!SUM OVER CELLS R
				do R = 1, nSC
					phaseR		= myExp( 	dot_product( kpts(:,ki) , Rcell(:,R) )			)
					call wHw(1,R,n,m, wnF, Hexp)	
					Htb(n,m)	= Htb(n,m) + phaseR * Hexp
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		!
		return
	end subroutine





	
	subroutine genUnkW(wnF, unkW)
		!generates lattice periodic functions unk from the Wannier functions wnf
		!	uses the coarse k point mesh
		complex(dp),	intent(in)		:: wnF(:,:,:)			!	wnF( 	nR, nSC, nWfs	)	
		complex(dp),	intent(out)		:: unkW(:,:,:)
		integer							:: n, R, qi, xi
		real(dp)						:: cellP
		!
		unkW	= dcmplx(0.0_dp)
		!GENERATE BLOCH LIKE FUNCTIONS
		do n = 1, nWfs
			do qi = 1, nQ
				do R = 1, nSC
					do xi = 1, nR
						cellP = -1.0_dp * dot_product(	qpts(:,qi) , 	rpts(:,xi)-Rcell(:,R)	)
						unkW(xi,qi,n) = unkW(xi,qi,n) + myExp(cellP) * wnF(xi,R,n) 	  
					end do
				end do
			end do
		end do
		!
		!
		return
	end subroutine



	!INTERPOLATION
	subroutine calcWannMat(wnF, Hw, Hwa, Aw, Fw)
		!calculates the Matrix elements in wannier gauge
		!	see Wang/Vanderbilt PRB 74, 195118 (2006)
		complex(dp),	intent(in)		:: wnF(:,:,:)		!wnF( 	nR	, nSC, nWfs	)
		complex(dp),	intent(out)		:: Hw(:,:,:)		!Hw(nKi, nWfs, nWfs)
		complex(dp),	intent(out)		:: Hwa(:,:,:,:)		!Hwa(3,nKi, nWfs, nWfs)
		complex(dp),	intent(out)		:: Aw(:,:,:,:)		!Aw(3,nKi, nWfs, nWfs)
		complex(dp),	intent(out)		:: Fw(:,:,:,:,:)	!Fw(3,3,nKi,nWfs,nWfs)
		complex(dp)						:: Rphase, Hexp
		real(dp)						:: rExp(3)
		integer							:: m, n, ki, R, i, j, myID
		!
		Hw			= dcmplx(0.0_dp)
		Hwa			= dcmplx(0.0_dp)
		Aw			= dcmplx(0.0_dp)
		Fw			= dcmplx(0.0_dp)
		rExp(3)		= 0.0_dp
		!



		!
		!
		!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) DEFAULT(SHARED) PRIVATE(Rphase, Hexp, rExp, m, n, ki, R, i, j, myID)
		do m = 1, nWfs
			do n =  1, nWfs
				do ki = 1, nK
					!FOURIER TRAFO (sequential)
					do R = 1, nSC
						Rphase	= myExp( dot_product( 	kpts(:,ki), Rcell(:,R)		) )
						!
						!ENERGY QUANTITIES
						call wHw(R0,R,n,m, wnf, Hexp)
						Hw(ki,n,m) 		= 			Hw(ki,n,m)			+ Rphase 						* Hexp
						Hwa(:,ki,n,m)	= 			Hwa(:,ki,n,m)		+ Rphase * i_dp * Rcell(:,R) 	* Hexp
						!
						!POSITIONAL QUANTITIES
						call wXw(R0,R,n,m, wnF, rExp(1:2))
						Aw(1:2,ki,n,m)	= 			Aw(1:2,ki,n,m)		+ Rphase 						* rExp(1:2)
						!
						do i = 1, 2
							do j = 1, 2
								Fw(i,j,ki,n,m)	= 	Fw(i,j,ki,n,m) 		+ Rphase * i_dp * Rcell(i,R) 	* rExp(j)
								Fw(i,j,ki,n,m)	= 	Fw(i,j,ki,n,m) 		- Rphase * i_dp * Rcell(j,R) 	* rExp(i)
							end do
						end do
					end do
					!
					myID	= OMP_GET_THREAD_NUM()
					write(*,'(a,i3,a,i3,a,i3,a,i3)')		"[calcWannMat,id=",myID,"]: I did m=",m," n=",n," ki=",ki
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		!
		!
		return
	end subroutine














!privat:
	subroutine genWannF2(unk, wnF)
		! generates wannier functions from (projected) bloch wavefunctions
		!
		complex(dp), 	intent(in)  	:: unk(:,:,:) ! lobWf(nRpts,nWfs)	
		complex(dp), 	intent(inout) 	:: wnF(:,:,:) ! wnF( 	nR, nSC, nWfs		)	
		integer 						:: n, Ri, xi, qi
		complex(dp)						:: phase
		real(dp)						:: nQreal
		!
		nQreal	= real(nQ,dp)
		wnF		= dcmplx(0.0_dp)
		!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) DEFAULT(SHARED) PRIVATE(n, Ri, xi, qi, phase) 
		do n = 1, nWfs
			do Ri = 1, nSC
				do xi = 1, nR
					do qi = 1 , nQ
						phase			= myExp(	-1.0_dp * dot_product(	qpts(:,qi) , Rcell(:,Ri)	) 	 )
						phase			= phase * myExp( dot_product( qpts(:,qi) , rpts(:,xi)))
						wnF(xi,Ri,n) = wnF(xi,Ri,n) + unk(xi,n,qi) * phase / dsqrt(nQreal)
					end do
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		!
		!
		return
	end subroutine


	subroutine wHw(R1,R2,n,m, wnF, res)
		!Hamilton operator expectation value wHw
		!
		!	wHw	= <R1,n|H(r)|R2,m>
		!		= -hbar**2/(2me)	<R1,n|\nabla**2 + V(r)|R2,m>
		!		= -0.5			[	<R1,n| dx**2 |R2,m>	 + <R1,n| dy**2 |R2,m>	+ <R1,n|V(r)|R2,m>]
		!
		integer,		intent(in)		:: R1, R2, n, m
		complex(dp),	intent(in)		:: wnF(:,:,:)	 !wnF( nRpts, nSupC,		nWfs)
		complex(dp),	intent(out)		:: res
		complex(dp)						:: resX, resY, resV, Vri
		real(dp),		allocatable		:: fx(:), fy(:), fv(:), fxx(:), fyy(:)
		integer							:: xi, yi, ri, lx, rx, ly, ry, at
		!
		allocate(	fx(nR)	)	
		allocate(	fy(nR)	)	
		allocate(	fv(nR)	)	
		allocate(	fxx(nR)	)
		allocate(	fyy(nR)	)		
		!
		!CALCULATE DERIVATIVES 
		do yi = 1, nRy
			do xi = 1, nRx
				ri 	= getRindex(xi,yi)
				!
				!X deriv
				call getNeighbX(xi,yi,lx,rx)
				fxx(ri)	= dconjg( wnF(ri,R1,n) )	* (		wnF(lx,R2,m) - 2.0_dp*wnF(ri,R2,m) + wnF(rx,R2,m)		)	/ dx**2
				!
				!Y deriv
				call getNeighbY(xi,yi,ly,ry)
				fyy(ri)	= dconjg( wnF(ri,R1,n))		* (		wnF(ly,R2,m) - 2.0_dp*wnF(ri,R2,m) + wnF(ry,R2,m)		)	/ dy**2
			end do
		end do
		!
		!FILL INTEGRATION ARRAY
		do ri = 1, nR
			fx(ri)	= -0.5_dp * dconjg(	wnF(ri,R1,n)	) * fxx(ri)
			fy(ri)	= -0.5_dp * dconjg(	wnF(ri,R1,n)	) * fyy(ri)
			!
			Vri	= dcmplx(0.0_dp)
			do at = 1, nAt
				if( insideAt(at,rpts(:,ri))	) then
					Vri = dcmplx( atPot(at) )
				end if	
			end do
			fv(ri)	= -0.5_dp * dconjg( wnF(ri,R1,n)	) * wnF(ri,R2,m) * Vri
		end do
		!
		resX	= nIntegrate(nR, nRx, nRy, dx, dy, fx)
		resY	= nIntegrate(nR, nRx, nRy, dx, dy, fy)
		resV	= nIntegrate(nR, nRx, nRy, dx, dy, fv)
		!
		res		= resX + resY + resV
		!
		!
		return
	end subroutine


	subroutine getNeighbX(xi,yi, lx, rx)
		!Helper for wHw
		integer,	intent(in)		:: xi,yi
		integer,	intent(out)		:: lx, rx
		!LEFT X NEIGHBOUR
		if(xi == 1) then
			lx	= getRindex( nRx	, yi	)
		else
			lx	= getRindex( xi-1	, yi	)
		end if
		!
		!RIGHT X
		if( xi == nRx ) then
			rx	= getRindex( 1		, yi	)
		else
			rx	= getRindex( xi+1	, yi 	)
		end if
		!
		!
		return
	end subroutine


	subroutine getNeighbY(xi, yi, ly, ry)
		!Helper for wHw
		integer,	intent(in)		:: xi,yi
		integer,	intent(out)		:: ly, ry
		!LEFT X NEIGHBOUR
		if(yi == 1) then
			ly	= getRindex( xi		, nRy	)
		else
			ly	= getRindex( xi		, yi-1	)
		end if
		!
		!RIGHT X
		if( yi == nRy ) then
			ry	= getRindex( xi		, 1	)
		else
			ry = getRindex( xi		, yi+1 	)
		end if
		!
		!
		return
	end subroutine


	subroutine wXw(R1,R2,n,m, wnF, res)
		!positional operator expectation value
		integer,		intent(in)		:: R1, R2, n, m
		complex(dp),	intent(in)		:: wnF(:,:,:)	 !wnF( nRpts, nSupC,		nWfs)
		real(dp),		intent(out)		:: res(2)		
		real(dp)				 		:: norm
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
	end subroutine


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
	end function





end module wannier