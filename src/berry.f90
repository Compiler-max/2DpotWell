module berry
	!module contains methods related to the Berry phase theory of modern polarization
	!	i.e. calculation of connection, velocities, curvatures and polarization
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, myLeviCivita, nIntegrate
	use sysPara
	implicit none

	private
	public	::	isKperiodic, calcWaveMat, calcConnOnCoarse, calcVeloMat, calcCurv
	contains
















!public
	integer function isKperiodic(unk)
		!	Test the condition u_nk = e^{iGr} u_nk+G
		!	see	King-Smith Vanderbilt PRB 48.7 4442 (1993)
		!
		!	the values at each x direction boundary (kx,1) is compared with (kx,nKy)
		!	analog for y direction
		!
		complex(dp),	intent(in)		:: unk(:,:,:)
		integer							:: n, qix, qiy, ri, totC, isX, isY
		real(dp)						:: Gx(2), Gy(2), val,dmax, avg
		complex(dp)						:: phase, u0, u1
		!
		Gx(1)		= 2.0_dp * PI_dp / aX
		Gx(2)		= 0.0_dp
		Gy(1)		= 0.0_dp
		Gy(2)		= 2.0_dp * PI_dp / aX
		isKperiodic = 0
		isX			= 0
		isY			= 0
		totC		= 0
		dmax		= 0.0_dp
		avg			= 0.0_dp
		!
		!
		do 	n = 1, nWfs 
			!CHECK IN X DIRECTION
			do  qix = 1, nQx 
				do ri = 1, nR
					phase 	= myExp( 	dot_product( Gx(:) , rpts(:,ri) )		)
					u0 		= unk( ri, getKindex(qix,1  ), n)
					u1		= unk( ri, getKindex(qix,nQy), n)
					val 	= abs(u0-phase*u1)
					if(  val > acc	) then 
						isKperiodic = isKperiodic 	+ 1
						isX			= isX			+ 1
						avg			= avg + val
						if( val > dmax) then
							dmax = val
						end if 
					end if
					totC = totC + 1
				end do	
			end do
			!
			!CHECK IN Y DIRECTION
			do  qiy = 1, nQy
				do ri = 1, nR
					phase 	= myExp( 	dot_product( Gy(:) , rpts(:,ri) )		)
					u0 		= unk( ri, getKindex(1  ,qiy), n)
					u1		= unk( ri, getKindex(nQx,qiy), n)
					val		= abs(u0-phase*u1)
					if( val > acc	) then 
						isKperiodic = isKperiodic	+ 1
						isY			= isY			+ 1
						avg			= avg + val
						if(val > dmax) then
							dmax = val
						end if 
					end if
					totC = totC + 1
				end do
			end do
		end do
		!
		!
		avg	= avg/ real(isKperiodic,dp)
		if(isKperiodic /= 0) then
			write(*,'(a,i8,a,i8,a,f16.12,a,f16.12)')	"[isKperiodic]: ",isKperiodic," of ", totC,&
														" unk test didnt pass; max delta=",dmax, &
														" average diff=",avg
			write(*,'(a,i8,a,i8)')	"[isKperiodic]:  fails along x boundary =",isX, "; along y=", isY 
		end if
		!
		return
	end function






	subroutine calcWaveMat(unk, Hw, Hwa, Aw, Fw)
		!calculates the matrices described in Wang/Vanderbilt PRB 74, 195118 (2006)
		!	via integration k space
		complex(dp),	intent(in)		:: unk(:,:,:)
		complex(dp),	intent(out)		:: Hw(:,:,:), Hwa(:,:,:,:), Aw(:,:,:,:), Fw(:,:,:,:,:)		!Hw(nKi, nWfs, nWfs)
		complex(dp)						:: Hnmq, qphase, HwTmp, HwaTmp, AwTmp, FwTmp
		real(dp),	allocatable		:: AwCoarse(:,:,:,:)
		integer							:: m, n, ki, qi, R, a,b, myID, nThreads, chunk
		!
		allocate(	AwCoarse(3,nQ,nWfs, nWfs)	)
		call calcConnOnCoarse(unk, AwCoarse)
		!
		Hw	= dcmplx(0.0_dp)
		Hwa	= dcmplx(0.0_dp)
		Aw	= dcmplx(0.0_dp)
		Fw	= dcmplx(0.0_dp)
		!

		
		!myID	= OMP_GET_THREAD_NUM()
		!nThreads= OMP_GET_NUM_THREADS()
		!write(*,'(a,i3,a,i3,a)')	"[calcWaveMat]: hello from omp thread ",myID,"of ",nThreads," threads"
		
		!

		!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) DEFAULT(SHARED) PRIVATE(m, n, ki, R, qi, Hnmq, qphase)
		do m = 1, nWfs
			do n =  1, nWfs
				do ki = 1, nK	
					!SUM OVER CELLS & COARSE K MESH (SEQUENTIAL)
					do R = 1, nSC
						do qi = 1, nQ
							qphase			= myExp( 	dot_product(	kpts(:,ki) - qpts(:,qi)	, Rcell(:,R)	)		)	
							!HAMILTONIAN QUANTITIES
							Hnmq 			= unHum(n,m,qi, unk)
							Hw(ki,n,m) 		= Hw(ki,n,m) 		+ 						qphase * Hnmq 				/ nQ
							Hwa(:,ki,n,m)	= Hwa(:,ki,n,m) 	+ i_dp * Rcell(:,R)	*	qphase * Hnmq 				/ nQ
							!POSITIONAL QUANTITIES
							Aw(:,ki,n,m)	= Aw(:,ki,n,m)		+ 						qphase * AwCoarse(:,qi,n,m) / nQ
							do b = 1, 3
								do a = 1, 3
									Fw(a,b,ki,n,m)	= Fw(a,b,ki,n,m)	+ i_dp * Rcell(a,R)	*	qphase * AwCoarse(b,qi,n,m)	/ nQ 
									Fw(a,b,ki,n,m)	= Fw(a,b,ki,n,m)	- i_dp * Rcell(b,R)	*	qphase * AwCoarse(a,qi,n,m)	/ nQ
								end do
							end do
						end do	
					end do
					!
					myID	= OMP_GET_THREAD_NUM()
					write(*,'(a,i3,a,i3,a,i3,a,i3)')		"[calcWaveMat,id=",myID,"]: I did m=",m," n=",n," ki=",ki
				end do
			end do
		end do
		!$OMP END PARALLEL DO
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
		!
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

						if(n==1 .and. m==1) then
							!write(*,*)	"*"
							!write(*,'(a,f16.12,a,f16.12,a,f16.12,a,f16.12,a)')"[calcConnOnCoarse]: qCe=(",qpts(1, qi					),	",",&	
							!										qpts(2,	qi						),"). ",	"oLap=(",dreal(one),", ",dimag(one),")."
							!!
							!write(*,'(a,f16.12,a,f16.12,a,f16.12,a,f16.12,a)')"[calcConnOnCoarse]: qxl=(",qpts(1, getKindex( qxl, qy )	),	",",&	
							!										qpts(2,	getKindex( qxl, qy )	),"). ",	"oLap=(",dreal(Mxl),", ",dimag(Mxl),")."
							!!
							!write(*,'(a,f16.12,a,f16.12,a,f16.12,a,f16.12,a)')"[calcConnOnCoarse]: qxr=(",qpts(1, getKindex( qxr, qy )	),	",",&	
							!										qpts(2,	getKindex( qxr, qy )	),"). ",	"oLap=(",dreal(Mxr),", ",dimag(Mxr),")."
							!!
							!write(*,'(a,f16.12,a,f16.12,a,f16.12,a,f16.12,a)')"[calcConnOnCoarse]: qyl=(",qpts(1, getKindex( qx ,qyl )	),	",",&	
							!										qpts(2,	getKindex( qx ,qyl )	),"). ",	"oLap=(",dreal(Myl),", ",dimag(Myl),")."
							!!
							!write(*,'(a,f16.12,a,f16.12,a,f16.12,a,f16.12,a)')"[calcConnOnCoarse]: qyr=(",qpts(1, getKindex( qx ,qyr )	),	",",&	
							!										qpts(2,	getKindex( qx ,qyr )	),"). ",	"oLap=(",dreal(Myr),", ",dimag(Myr),")."
							
							!write(*,*)	"*"
							!write(*,*)"[calcConnOnCoarse]: qCe=(",qpts(1, qi					),	",",&	
							!										qpts(2,	qi						),"). ",	"oLap=(",dreal(one),", ",dimag(one),")."
							!!
							!write(*,*)"[calcConnOnCoarse]: qxl=(",qpts(1, getKindex( qxl, qy )	),	",",&	
							!										qpts(2,	getKindex( qxl, qy )	),"). ",	"oLap=(",dreal(Mxl),", ",dimag(Mxl),")."
							!!
							!write(*,*)"[calcConnOnCoarse]: qxr=(",qpts(1, getKindex( qxr, qy )	),	",",&	
							!										qpts(2,	getKindex( qxr, qy )	),"). ",	"oLap=(",dreal(Mxr),", ",dimag(Mxr),")."
							!!
							!write(*,*)"[calcConnOnCoarse]: qyl=(",qpts(1, getKindex( qx ,qyl )	),	",",&	
							!										qpts(2,	getKindex( qx ,qyl )	),"). ",	"oLap=(",dreal(Myl),", ",dimag(Myl),")."
							!!
							!write(*,*)"[calcConnOnCoarse]: qyr=(",qpts(1, getKindex( qx ,qyr )	),	",",&	
							!										qpts(2,	getKindex( qx ,qyr )	),"). ",	"oLap=(",dreal(Myr),", ",dimag(Myr),")."

						end if
						!
						!write(*,'(a,f15.12,a,f15.12)')"[calcConn]: Mxl=",dreal(Mxl),"+i*",dimag(Mxl)
						!FD SUM OVER NEAREST NEIGHBOURS
						A(1:2,qi,n,m) = A(1:2,qi,n,m) + wbx * bxl(1:2) * dimag( Mxl - one )
						A(1:2,qi,n,m) = A(1:2,qi,n,m) + wbx * bxr(1:2) * dimag( Mxr - one )
						A(1:2,qi,n,m) = A(1:2,qi,n,m) + wby * byl(1:2) * dimag( Myl - one )
						A(1:2,qi,n,m) = A(1:2,qi,n,m) + wby * byr(1:2) * dimag( Myr - one )
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
		!
		!DEBUG
		avg	= avg / real(found,dp)
		write(*,'(a,i6,a,i6,a,f16.12,a,f16.12)')	"[calcConnOnCoarse]: ",found," of ",tot,&
										" checked unk functions had normalization issues;  max delta=",dmax,&
										" avg diff=",avg
		!
		return
	end subroutine


	subroutine calcVeloMat(unk, veloBwf, Velo)
		!calculates matrix elements of the velocity operator
		!	velocity operator is analytically applied to the plane wave basis 
		!	and then weighted by the basCoeff obtained from the solver and stored in veloBwf
		!	
		complex(dp),		intent(in)		:: unk(:,:,:), veloBwf(:,:,:)	!	unk(nR, nK, nWfs) , veloBwf(nR,nK ,2*nWfs)
		complex(dp),		intent(out)		:: Velo(:,:,:,:)   !Velo(3,	nK,	nWfs	, nwFs)		
		integer								:: qi, n,m, ri
		complex(dp),		allocatable		:: fx(:), fy(:)
		!
		allocate(	fx(nR)	)
		allocate(	fy(nR)	)
		!
		do m = 1, nWfs
			do n = 1, nWfs
				do qi = 1, nQ
					!FILL INTEGRATION ARRAY
					do ri = 1, nWfs
						fx(ri)	= dconjg(	myExP( 	dot_product( qpts(:,qi), rpts(:,ri) ))	)	*unk(ri,qi,n)	* veloBwf(ri, qi,	m		)
						fy(ri)	= dconjg(	myExP( 	dot_product( qpts(:,qi), rpts(:,ri) ))	)	*unk(ri,qi,n)	* veloBwf(ri, qi,	nWfs+m	)
					end do
					!INTEGRATE
					Velo(1,qi,n,m)	= nIntegrate(nR, nRx, nRy, dx, dy, fx)
					Velo(2,qi,n,m)	= nIntegrate(nR, nRx, nRy, dx, dy, fy)
					Velo(3,qi,n,m) 	= dcmplx(0.0_dp)
				end do
			end do
		end do
		!	
		!
		return
	end subroutine


	subroutine approxVelo(unk)
		complex(dp),	intent(in)		:: unk(:,:,:)



		return
	end subroutine




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
		integer							:: n, qi, a,b,c
		!
		Fcurv = dcmplx(0.0_dp)
		do n = 1, nWfs
			do qi = 1, nQ
				do c = 1,3
					do b= 1,3
						do a=1,3
							if( myLeviCivita(a,b,c) /= 0) then
								Fcurv(c,qi,n) = Fcurv(c,qi,n) + myLeviCivita(a,b,c) * omega(n,a,b,qi,Velo, En)
							end if
						end do 
					end do
				end do
			end do
		end do
		!
		!
		return
	end subroutine


























!privat
	complex(dp) function unHum(n,m,qi, unk )
		!calculates the matrix element
		!
		!	<u_nq|H(q)|u_mq>
		!on the coarse k point mesh
		integer,		intent(in)		:: n,m, qi
		complex(dp),	intent(in)		:: unk(:,:,:)
		complex(dp),	allocatable		:: f(:)
		complex(dp)						:: phase, ham
		integer							:: ri, at
		!
		allocate(	f(nR)	)
		!
		do ri = 1, nR	
			!SET UP HAMILTONIAN AT CURRENT ri
			ham 	= 0.5_dp * dot_product(qpts(:,qi),qpts(:,qi))
			do at = 1, nAt
				if(insideAt(at,rpts(:,ri)))	then
					ham = ham + atPot(at)
				end if
			end do
			phase	= myExp( dot_product(qpts(:,qi), rpts(:,ri)	)	)
			ham 	= dconjg(phase) * ham * phase
			!FILL INTEGRATION ARRAY
			f(ri)	= dconjg(	unk(ri,qi,n)	) * ham * unk(ri,qi,m)
		end do
		!
		unHum = nIntegrate(nR, nRx,nRy, dx,dy, f)
		!
		!
		return
	end function










	real(dp) function omega(n,a,b,qi,Velo,En)
		!returns curvature tensor element a,b
		!see Wang/Vanderbilt PRB 74, 195118 (2006) eq.(9)
		!
		integer,		intent(in)		:: n, a, b, qi
		complex(dp),	intent(in)		:: Velo(:,:,:,:)!Velo(3,nK,nWfs,nWfs)
		real(dp),		intent(in)		:: En(:,:)!En(nK,nWfs)
		integer							:: m
		!
		omega	= 0.0_dp
		do m = 1, nWfs
			if(m /= n) then
				omega = omega +		dimag(	Velo(a,qi,n,m) * Velo(b,qi,m,n)	) 	/ 	( En(qi,m) - En(qi,n) )**2
			end if 
		end do		
		omega	= -2.0_dp * omega
		!
		!
		return
	end function


	!CONNECTION HELPERS
	complex(dp) function UNKoverlap(n, m, qi, knb, unk)
		!HELPER for calcConn
		!calculates the overlap between unk at qi and at a neigbhouring k point knb
		!	integration only over the first unit cell
		!
		integer,		intent(in)		:: n, m, qi, knb
		complex(dp),	intent(in)		:: unk(:,:,:)  !unk(	nR, nK, nWfs/nG	)
		complex(dp),	allocatable		:: f(:)
		integer							:: xi,yi,ri,rloc, nRx1, nRy1, nR1
		!
		!Set integration range to first unit cell
		!nRx1 	= int(		real(nRx,dp) / real(nSCx,dp)		)
		!nRy1 	= int(		real(nRy,dp) / real(nSCy,dp)		)
		!nR1		= nRx1 * nRy1 
		!allocate(	f(nR1)	)
		!!
		!!fill integration array
		!f 		= dcmplx(0.0_dp)
		!do yi = 1, nRy1
		!	do xi = 1, nRx1
		!		ri		= getRindex(xi,yi)			!overall index, to get correct position from unk
		!		rloc 	= (yi-1) * nRx1 + xi		!for mapping to f array
		!		f(rloc)	= dconjg( unk(ri,qi,n) ) * unk(ri,knb,m)
		!		!write(*,'(a,f10.6,a,f10.6)')	"[overlap] f=",dreal(f(rloc)),"+i*",dimag(f(rloc))
		!	end do
		!end do
		!
		allocate( f(nR)	)
		do ri = 1, nR
			f(ri)	= dconjg( unk(ri,qi,n) ) * unk(ri,knb,m)
		end do



		!integrate
		!UNKoverlap = nIntegrate(nR1, nRx1, nRy1, dx, dy, f	)
		UNKoverlap = sum( f ) / real(nR,dp)
		!write(*,'(a,f10.6,a,f10.6)')"[overlap]=",dreal(overlap),"+i*",dimag(overlap)
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






!	subroutine gaugeConn(Aw, unkW, EnH, Uh, Ah)
!		!	transforms connection from Wannier gauge Aw back to Hamiltonian gauge Ah
!		!	
!		!	Ah 		= Ahbar + i*Dh
!		!	Ahbar 	= (U*) Aw U 
!		!	Dh_(n,m)= (U*) Hw_(n,m) U / (Em - En)
!		!
!		!	see Wang/Vanderbilt PRB 74, 195118 (2006) 
!		real(dp),		intent(in)		:: Aw(:,:,:,:), EnH(:,:)		!Aconn(	3	, nK, nWfs	), En(nK, nWfs)		
!		complex(dp),	intent(in)		:: unkW(:,:,:), Uh(:,:,:)	!unk(	nR	, nK, nWfs	), 	Uh(	nWfs, nWfs,	nK)		
!		real(dp),		intent(out)		:: Ah(:,:,:,:)				!Ah(3, nK, nWfs, nWfs)
!		complex(dp),	allocatable		:: Dh(:,:,:), Atmp(:,:,:)				!Dh(3,		nWfs, nWfs)
!		integer							:: ki, n, m, i
!		!
!		allocate(		Dh(		size(Aw,1),		size(Aw,3), size(Aw,4)	)		)
!		allocate(		Atmp(	size(Aw,1),		size(Aw,3),	size(AW,4)	)		)
!		Ah	= dcmplx(0.0_dp)
!		!
!		!
!		do ki = 1, nK
!			call calcAHbar(ki, Aw, Uh, Atmp)
!			call calcDh(ki, unkW, EnH, Dh)
!
!			Ah(:,ki,:,:)	= Atmp(:,:,:) + i_dp * Dh(:,:,:)	!test
!			!do m = 1, nWfs
!			!	do n = 1, nWfs
!			!		do i = 1, 3
!			!			Ah(i,ki,n,m)	= Ah(i,n,m) + Dh(i,n,m)
!			!		end do
!			!	end do
!			!end do
!		end do
!		!
!		!
!		return
!	end
!
!
!
!	subroutine calcAHbar(ki, Aw, Uh, Ahbar)
!		!rotates A from wannier gauge to hamiltonian gauge	
!		!	Ahbar 	= (U*) Aw U  
!		!
!		integer,		intent(in)		:: ki
!		real(dp),		intent(in)		:: Aw(:,:,:,:)
!		complex(dp),	intent(in)		:: Uh(:,:,:)
!		complex(dp),	intent(out)		:: Ahbar(:,:,:)
!		complex(dp),	allocatable		:: Uconjg(:,:), Atmp(:,:)
!		integer							:: i
!		!
!		allocate(	Uconjg( size(Uh,2), size(Uh,2) 	)		)
!		allocate(	Atmp( 	size(Aw,3), size(Aw,4)	)		)
!		!
!		!
!		Uconjg	= dconjg( transpose(Uh(:,:,ki)	)	)
!		do i = 1, 3
!			Ahbar(i,:,:)	= matmul( 	dcmplx(Aw(i,ki,:,:))	, 	Uh(:,:,ki)				) 
!			Ahbar(i,:,:)	= matmul(	Uconjg(:,:)				,	Ahbar(i,:,:)			)
!			!Atmp(:,:)		= Ahbar(i,:,:)
!			!Ahbar(i,:,:)	= matmul(	Uconjg(:,:)				,	Atmp(:,:)				)
!		end do
!		!
!		return
!	end
!
!
!	subroutine calcDh(ki, unkW, EnH, Dh)
!		!
!		!
!		!
!		integer,		intent(in)		:: ki
!		complex(dp),	intent(in)		:: unkW(:,:,:)		!unk(nR,nK,nWfs)
!		real(dp),		intent(in)		:: EnH(:,:)			!EnH(nK, nWfs)		
!		complex(dp),	intent(out)		:: Dh(:,:,:)		!Dh(	3	, nWfs, nWfs	)	
!		complex(dp),	allocatable		:: Hbar(:,:,:)
!		integer							:: n,m,i
!		!
!		allocate(		Hbar(size(Dh,1), size(Dh,2), size(Dh,3))		)
!		!
!		Dh	= dcmplx(0.0_dp)
!		call calcHbar(ki, unkW, Hbar)
!		do m = 1, size(Dh,3)
!			do n =1, size(Dh,2)
!				if( n /= m) then
!					do i = 1, 3
!						Dh(i,n,m)	= Hbar(i,n,m)	/ (	EnH(ki,m)	-	EnH(ki,n)	)	
!					end do
!				end if
!			end do
!		end do
!		!
!		!
!		return
!	end 
!
!
!
!	!call calcHbar(ki, unkW, Hbar)
!
!
!
!





	













end module berry