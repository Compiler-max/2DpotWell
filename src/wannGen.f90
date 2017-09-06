module wannGen
	!module contains all the functions to generate the wannier functions
	!	this includes the calculation of the orthonormalized projected bloch wavefunctions
	!	and the FT of the bwfs to calculate the wannier functions
	use mathematics, only: dp, PI_dp, acc, myExp, nIntegrate, myMatInvSqrt, isUnit, isIdentity
	use sysPara
	implicit none	
	
	private
	public ::	projectBwf, genWannF			














	contains
!public:
	subroutine projectBwf(qi, bWf, loBwf, U, failCount, smin, smax)	!projectBwf(qi, bWf, loBwf, U(:,:), failCount, smin, smax)
		!does the projection onto Loewdin-orthonormalized Bloch-like states
		!see Marzari, Vanderbilt PRB 56, 12847 (1997) Sec.IV.G.1 for detailed description of the method 
		integer		,	intent(in)		:: qi
		complex(dp)	,	intent(in)		:: bWf(:,:)   ! bWf(nR,nG)
		complex(dp)	,	intent(out)		:: loBwf(:,:), U(:,:)   !U(nWfs,nWfs)
		integer,		intent(inout)	:: failCount
		real(dp),		intent(inout)	:: smin, smax
		complex(dp)	,	allocatable		:: gnr(:,:), A(:,:), S(:,:), phi(:,:)
		integer							:: n, nBands, xi
		!
		allocate(	gnr(nR,nWfs)		)
		allocate( 	S(nWfs,nWfs)		)
		allocate( 	A(nWfs,nWfs)		)
		loBwf 	= dcmplx(0.0_dp)
		U 		= dcmplx(0.0_dp)
		!
		!SET UP ROTATION MATRIX U
		call genTrialOrb(gnr)
		call calcAmat(bwf,gnr, A)
		call calcInvSmat(A, S, smin, smax)
		U = matmul(S, A)
		!
		!ROTATE BLOCH STATES
		do xi = 1, nR
			loBwf(xi,:) = matmul(	U , bWf(xi,1:nWfs)	)
		end do
		!
		!DEBUGGING
		if( .not. isUnit(U) 			) then
			write(*,*)"[projectBwf]: U matrix is not a unitary matrix !"
		end if
		if(	 .not. isOrthonorm(lobWf)	) then
			write(*,'(a,f16.12)')"[projectBwf]: loBwf is NOT a orthonormal basis set, accuracy=",acc
			failCount = failCount + 1
		end if
		!
		return
	end subroutine
		

	subroutine genWannF(qi, bWf, wnF)
		! generates wannier functions from (projected) bloch wavefunctions
		!
		integer,		intent(in)		:: qi
		complex(dp), 	intent(in)  	:: bWf(:,:) ! lobWf(nRpts,nWfs)	
		complex(dp), 	intent(inout) 	:: wnF(:,:,:) ! wnF( 	nR, nSC, nWfs		)	
		integer 				 :: n,xi,Ri
		real(dp)				 :: cellP
		!
		do n = 1, nWfs
			do Ri = 1, nSC
				cellP = -1.0_dp * dot_product(	qpts(:,qi) , Rcell(:,Ri)	) 
				do xi = 1, nR		
						wnF(xi,Ri,n) = wnF(xi,Ri,n) + bWf(xi,n) * myExp(cellP) / real(nK,dp)
				end do
			end do
		end do
		!
		return
	end subroutine

















!privat:
	subroutine genTrialOrb(gnr)
		complex(dp),	intent(out)	:: gnr(:,:)       !gnr( nR, nWfs)	
		real(dp)					:: posX, posY, xc, k, L
		complex(dp)					:: A
		integer						:: n, ri, at
		!
		gnr = dcmplx(0.0_dp)
		do n = 1, nWfs
			do ri = 1, nR
				do at 	= 1, nAt
					if(	insideAt( at , rpts(:,ri) )	) then
						gnr(ri,n) = gVal(at, n, ri)
						!write(*,'(a,i2,a,f6.3,a,f6.3,a,f6.3,a,f6.3,a,f6.3,a,f6.3)')	"[genTrialOrb]: gnr(at=",&
						!			at," r=(",rpts(1,xi),",",rpts(2,xi),&
						!		")= (",posX,",",posY,"))= ",dreal(gnr(xi,n)),"+i*",dimag(gnr(xi,n))
					end if
				end do
				!
			end do
		end do
		!
		return
	end subroutine


	complex(dp) function gVal(at, n, xi)
		!wrapper for calling different trial Orbitals according to input file var trialOrbSw
		integer, 	intent(in)		:: at, n, xi
		!
		select case(trialOrbSw)
			case (1)
				gVal = infPotWell(at, n, xi)
			case default
				gVal = dcmplx(	trialOrbVAL(at) /sum( trialOrbVAL )		)
		end select
		!
		return
	end function


	complex(dp) function infPotWell(at, n, xi)
		integer, 	intent(in)		:: at, n, xi
		real(dp)					:: xc(dim), xrel(dim), k(dim), L(dim)
		complex(dp)					:: A
		!
		!vAvg = abs( sum(vPot) ) * 1.0_dp/ nWells
		L 		= 2.0_dp * atR(:,at)
		xc 		= atPos(:,at)
		k 		= n*PI_dp / L(:)
		A 		= dsqrt(	2.0_dp / product(L)		)
		xrel	= rpts(:,xi) - xc + L / 2.0_dp
		!
		infPotWell = A * dcmplx(		 dsin(	dot_product( k , xrel)		) 	)
		!
		return
	end function



	subroutine calcAmat(bWf,gnr, A)
		!calculates the inner product matrix of bwfs and trial orbital gnr
		complex(dp),	intent(in)	:: bWf(:,:), gnr(:,:) 
		complex(dp),	intent(out)	:: A(:,:) !A(nWfs,nWfs)
		complex(dp),	allocatable	:: f(:)
		complex(dp)					:: overlap
		integer						:: m,n, xi
		allocate(	f(nR)	)
		!
		A = dcmplx(0.0_dp)
	
		do n = 1, nWfs
			!do m = 1, nWfs 
				!if band m is to contribute to wannier function n then calculate overlap
				!if( minBand(n) <= m .and. m <= maxBand(n)) then
			f = dcmplx(0.0_dp)
			do xi = 1, nR
				!write(*,'(a,i3,a,i4,a,f10.5,a,f10.5,a)')"g(n=", n,", xi=",xi,")= " , dreal( gnr(xi,n) ), "+ i* ",dimag( gnr(xi,n) ),")"
				f(xi) = dconjg(	bWf(xi,n) ) * gnr(xi,n)
				!write(*,'(a,i3,a,i4,a,f10.5,a,f10.5,a)')"f(n=", n,", xi=",xi,")= " , dreal( f(xi) ), "+ i* ",dimag( f(xi) ),")"
			end do
			A(n,n) = nIntegrate(nR, nRx, nRy, dx, dy, f)
			
			!write(*,'(a,i2,a,f15.10,a,f15.10)')	"[calcAmat]: n=,",n," A(n,n)=",dreal(A(n,n))," +i*",dimag(A(n,n))
				!else
				!	A(m,n) = dcmplx(0.0_dp)
				!end if
				!
			!end do
		end do
		!
		return
	end subroutine



	subroutine calcInvSmat(A, S, smin, smax)
		!calculates the sqrt inv. of the overlap matrix S
		complex(dp),	intent(in)		:: A(:,:)
		complex(dp),	intent(out)		:: S(:,:)
		real(dp),		intent(inout)	:: smin, smax
		complex(dp),	allocatable		:: Acon(:,:),Sold(:,:), Ssqr(:,:)
		integer							:: m,n, i,j
		real(dp)						:: mi, ma
		m = size(A,1)
		n = size(A,2)
		allocate(	Acon(n,m)	)
		allocate(	Sold(n,n)	)
		allocate(	Ssqr(n,n)	)
		!
		S = dcmplx(0.0_dp)
		Acon = transpose(A)
		Acon = dconjg(Acon)
		S = matmul(Acon , A)
		Sold = S
		!
		call myMatInvSqrt(S, mi, ma)

		if( mi < smin ) then
			smin = mi
		end if
		if( ma > smax ) then
			smax = ma
		end if

		!DEBUG: Check if S is really the inverse sqrt
		Ssqr = matmul(S,S)
		Sold = matmul(Sold,Ssqr)
		
		if( .not. isUnit(Sold) ) then
			write(*,*)"[calcInvSmat]: problem with mat inversion, seems to be not inverse square root"
		end if
		!
		!
		return
	end subroutine


	logical function isOrthonorm(loBwf)
		!cheks wether loBwf is orthonormal by calculating the overlap matrix elements
		!	<psi_n,k|psi_n,k'> = nSC * \delta(n,m)
		!	WARNING: DOES NOT CHECK ORTHONORMALIZATION BETWEEN K POINTS
		complex(dp),	intent(in)		:: lobWf(:,:)   !lobWf(nR,nWfs)
		complex(dp),	allocatable		:: f(:)
		complex(dp)						:: oLap
		
		integer							:: n, m, ri
		logical							:: realist, imaginist
		allocate( f(nR)		)
		!
		isOrthonorm = .true.
		n 			= 1
		!
		do  while(n <= nWfs .and. isOrthonorm) 
			!DIAGONAL ELEMENTS
			!
			!
			!PREPARE INTEGRATION ARRAY
			do ri = 1, nR
				f(ri) = dconjg( lobWf(ri,n) ) * lobWf(ri,n)
			end do
			oLap 		= nIntegrate(nR, nRx, nRy, dx, dy, f)
			oLap			= oLap / dcmplx( nSC	)
			!CHECK CONDITIONS
			realist 	= abs(			abs(dreal(oLap))		-		1.0_dp		) 	> acc
			imaginist	= abs(			dimag(oLap)								)	> acc
			if( realist .or. imaginist ) then
				write(*,'(a,i2,a,f15.10,a,f15.10)')"[isOrthonorm]: overlap(n,n=",n,") = ",dreal(oLap),"+i*",dimag(oLap)
				isOrthonorm = .false.
			end if
			!
			!
			!OFF DIAGONAL
			m = 1
			do while(m <= nWfs	.and. isOrthonorm)
				if(n /= m) then
					!PREPARE INTEGRATION ARRAY
					do ri = 1, nR
						f(ri) = dconjg( lobWf(ri,n) ) * lobWf(ri,m)
					end do
					oLap 	= nIntegrate(nR, nRx, nRy, dx, dy, f)
					oLap	= oLap / dcmplx( nSC	)
					!CHECK CONDITION
					if(abs(oLap) > acc) then
						write(*,'(a,i2,a,i2,a,f15.10,a,f15.10)')"[isOrthonorm]: overlap(",n,",",m,") = ",dreal(oLap),"+i*",dimag(oLap)
						isOrthonorm = .false.
					end if
				end if
				m = m +1
			end do
			n = n + 1
		end do
		!
		return
	end function
		


end module wannGen