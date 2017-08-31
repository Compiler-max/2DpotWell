module berry
	!module contains methods related to the Berry phase theory of modern polarization
	!	i.e. calculation of connection, velocities, curvatures and polarization
	use mathematics,	only:	dp, PI_dp, i_dp, myExp, myLeviCivita, nIntegrate
	use sysPara,		only: 	readInp, getKindex, getRindex, &
									dim, aX, aY,vol, nAt, atR, atPos, atPot,&
									nG, nG0, Gcut, nK, nKx, nKy, nWfs, nSC, nSCx, nSCy, nR, nRx, nRy, dx, dy, dkx, dky, &
									Gvec, atPos, atR, kpts, rpts, gaugeSwitch
	implicit none

	private
	public	::	calcConn, calcCurv, calcVeloMat, calcPolViaA

	contains








!public
		subroutine calcConn(unk,nxk, nyk, A)
		!finite difference on lattice periodic unk to calculate the Berry connection A
		!	A_n(k) 	= <u_n(k)|i \nabla_k|u_n(k)>
		!		 	= i  <u_n(k)| \sum_b{ w_b * b * [u_n(k+b)-u_n(k)]}
		!			= i \sum_b{		w_b * b * [  <u_n(k)|u_n(k+b)> -  <u_n(k)|u_n(k)>]		}
		!
		! see Mazari, Vanderbilt PRB.56.12847 (1997), Appendix B
		!
		complex(dp),	intent(in)		:: unk(:,:,:)		!unk(	nR, nK/nKw, nWfs/nG	)
		integer,		intent(in)		:: nxk, nyk
		real(dp),		intent(out)		:: A(:,:,:)			!Aconn(	3,nK, nWfs)		)	
		complex(dp)						:: Mxl, Mxr, Myl, Myr, M, one
		integer							:: n, Z, ki, kx, ky, kxl, kxr, kyl, kyr
		real(dp)						:: thres, wbx,wby, bxl(2), bxr(2), byl(2), byr(2) !for nearest neighbours, assuming cubic mesh
		!
		thres	= 1e-3_dp
		A 		= 0.0_dp
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
			do kx = 1, nxk
				kxl	= getLeft(kx,nxk)
				kxr	= getRight(kx,nxk)
				!
				do ky = 1, nyk
					kyl	= getLeft(ky,nyk)
					kyr = getRight(ky,nyk)
					ki	= getKindex(kx,ky)
					!
					!OVERLAP TO NEAREST NEIGHBOURS
					one = UNKoverlap(	n, 		ki		, 		ki					, unk	)
					Mxl	= UNKoverlap(	n, 		ki		, getKindex( kxl, ky ) 		, unk	) 
					Mxr	= UNKoverlap(	n, 		ki		, getKindex( kxr, ky )		, unk	)
					Myl	= UNKoverlap(	n, 		ki		, getKindex( kx ,kyl )		, unk	)
					Myr	= UNKoverlap(	n, 		ki		, getKindex( kx ,kyr )		, unk	)
					!
					!write(*,'(a,f15.12,a,f15.12)')"[calcConn]: Mxl=",dreal(Mxl),"+i*",dimag(Mxl)
					!FD SUM OVER NEAREST NEIGHBOURS
					A(1:2,ki,n) = A(1:2,ki,n) + wbx * bxl(:) * dimag( Mxl - one )
					A(1:2,ki,n) = A(1:2,ki,n) + wbx * bxr(:) * dimag( Mxr - one )
					A(1:2,ki,n) = A(1:2,ki,n) + wby * byl(:) * dimag( Myl - one )
					A(1:2,ki,n) = A(1:2,ki,n) + wby * byr(:) * dimag( Myr - one )
					!FD SUM OVER NEAREST NEIGHBOURS
					!A(:,ki,n) = A(:,ki,n) + wbx * bxl(:) * dimag(	log( Mxl ) )
					!A(:,ki,n) = A(:,ki,n) + wbx * bxr(:) * dimag(	log( Mxr ) )
					!A(:,ki,n) = A(:,ki,n) + wby * byl(:) * dimag(	log( Myl ) )
					!A(:,ki,n) = A(:,ki,n) + wby * byr(:) * dimag(	log( Myr ) )
					!
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




	subroutine calcCurv(En, Velo, Fcurv)
	!Calculates the connection via Kubo formula on matrix elements of velocity operator
	!see Wang/Vanderbilt PRB 74, 195118 (2006) eq.(5)
	!
	!	F_n,c(k) = \sum{a,b} leviCivi(a,b,c) F_n,{a,b}(k)
	!
	!
		real(dp),		intent(in)		:: En(:,:)		!En(nK,nWfs)
		complex(dp),	intent(in)		:: Velo(:,:,:,:) !Velo(3,nK,nWfs,nWfs)
		real(dp),		intent(out)		:: Fcurv(:,:,:)  !Fcurv(3,nK,nWfs)
		integer							:: n, ki, a,b,c
		!
		Fcurv = dcmplx(0.0_dp)
		do n = 1, nWfs
			do ki = 1, nK
				do c = 1,3
					do b= 1,3
						do a=1,3
							if( myLeviCivita(a,b,c) /= 0) then
								Fcurv(c,ki,n) = Fcurv(c,ki,n) + myLeviCivita(a,b,c) * omega(n,a,b,ki,Velo, En)
							end if
						end do 
					end do
				end do
			end do
		end do
		!
		!
		return
	end



	subroutine calcVeloMat(unk, veloBwf, Velo)
		!calculates matrix elements of the velocity operator
		!	velocity operator is analytically applied to the plane wave basis 
		!	and then weighted by the basCoeff obtained from the solver and stored in veloBwf
		!	
		complex(dp),		intent(in)		:: unk(:,:,:), veloBwf(:,:,:)	!	unk(nR, nK, nWfs) , veloBwf(nR,nK ,2*nWfs)
		complex(dp),		intent(out)		:: Velo(:,:,:,:)   !Velo(3,	nK,	nWfs	, nwFs)		
		integer								:: ki, n,m, ri
		complex(dp),		allocatable		:: fx(:), fy(:)

		allocate(	fx(nR)	)
		allocate(	fy(nR)	)

		
		do m = 1, nWfs
			do n = 1, nWfs
				do ki = 1, nK
					!FILL INTEGRATION ARRAY
					do ri = 1, nWfs
						fx(ri)	= myExP( 	dot_product( kpts(:,ki), rpts(:,ri) )		)	*unk(ri,ki,n)	* veloBwf(ri, ki,	m		)
						fy(ri)	= myExP( 	dot_product( kpts(:,ki), rpts(:,ri) )		)	*unk(ri,ki,n)	* veloBwf(ri, ki,	nWfs+m	)
					end do
					!INTEGRATE
					Velo(1,ki,n,m)	= nIntegrate(nR, nRx, nRy, dx, dy, fx)
					Velo(2,ki,n,m)	= nIntegrate(nR, nRx, nRy, dx, dy, fy)
					Velo(3,ki,n,m) 	= dcmplx(0.0_dp)
				end do
			end do
		end do
		!	
		!
		return
	end


	subroutine calcPolViaA(A, pElA)
		!calculates the polarization by integrating connection over the brillouin zone
		! r_n 	= <0n|r|0n> 
		!		=V/(2pi)**2 \integrate_BZ <unk|i \nabla_k|unk>
		!		=V/(2pi)**2 \integrate_BZ A(k)
		real(dp),		intent(in)		:: A(:,:,:)			!A(2,	 nK, nWfs	)	
		real(dp),		intent(out)		:: pElA(:)
		complex(dp)	,	allocatable		:: val(:)
		real(dp)						:: thres
		integer							:: n, ki
		!
		allocate(	val( size(A,1) )	)
		val		= dcmplx(0.0_dp)
		pElA	= 0.0_dp
		thres	= 1e-10_dp
		!
		!SUM OVER K SPACE AND OVER STATES
		do n 	= 1, size(A,3)
			do ki = 1, size(A,2)
				val(:)	= val(:) + A(:,ki,n)
			end do
		end do
		!
		!NORMALIZE
		val		= val / real(size(A,2),dp)
		!
		!HARVEST
		pElA	= val !/ (aX*aY)
		!
		!MOD TO FIRST UNIT CELL
		pElA(1)	= mod( pElA(1), aX )
		pElA(2) = mod( pElA(2), aY )
		!
		!DEBUGGING
		if(		dimag( val(1) ) > thres 	.or. 	dimag( val(2) ) > thres		) then
			write(*,*)"[calcPolViaA]: non zero imaginary part of polarization"
		end if
		return
	end







!privat
	real(dp) function omega(n,a,b,ki,Velo,En)
		!returns curvature tensor element a,b
		!see Wang/Vanderbilt PRB 74, 195118 (2006) eq.(9)
		!
		integer,		intent(in)		:: n, a, b, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)!Velo(3,nK,nWfs,nWfs)
		real(dp),		intent(in)		:: En(:,:)!En(nK,nWfs)
		integer							:: m
		!
		omega	= 0.0_dp
		do m = 1, nWfs
			if(m /= n) then
				omega = omega +		dimag(	Velo(a,ki,n,m) * Velo(b,ki,m,n)	) 	/ 	( En(ki,m) - En(ki,n) )**2
			end if 
		end do		
		omega	= -2.0_dp * omega
		!
		!
		return
	end


	!CONNECTION HELPERS
	complex(dp) function UNKoverlap(n, ki, knb, unk)
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
		UNKoverlap = nIntegrate(nR1, nRx1, nRy1, dx, dy, f	)
		
		!write(*,'(a,f10.6,a,f10.6)')"[overlap]=",dreal(overlap),"+i*",dimag(overlap)
		!
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





end module berry