module projection
	!module contains all the functions to generate the wannier functions
	!	this includes the calculation of the orthonormalized projected bloch wavefunctions
	!	and the FT of the bwfs to calculate the wannier functions
	use omp_lib
	use mathematics, 	only: 	dp, PI_dp, i_dp, acc, machineP,&
								myExp, nIntegrate, eigSolver,  mySVD, myMatInvSqrt, isUnit, isIdentity, isHermitian
	use sysPara
	use blochWf,		only:	calcBasis
	use output,			only:	writeInterpBands
	implicit none	
	
	private
	public ::	projectUnk		














	contains
!public:

	subroutine projectUnk(ckH, ckW, Uq)	!projectBwf(qi, bWf, loBwf, U(:,:), failCount, smin, smax)
		!does the projection onto Loewdin-orthonormalized Bloch-like states
		!see Marzari, Vanderbilt PRB 56, 12847 (1997) Sec.IV.G.1 for detailed description of the method 
		complex(dp)	,	intent(in)		:: ckH(:,:,:)  ! 	ck(nG,nBands  ,	nQ)		
		complex(dp)	,	intent(out)		:: ckW(:,:,:), Uq(:,:,:)	! Uq(nBands	,nWfs, 	nQ	)	
		real(dp),		allocatable		:: EnP(:,:)	
		complex(dp)	,	allocatable		:: loBwf(:,:), gnr(:,:), A(:,:), Ham(:,:), tmp(:,:)
		
		integer							:: qi, n, m, gi
		!
		allocate(	loBwf(nR,nWfs)		)
		allocate(	gnr(nR,nWfs)		)
		allocate( 	A(nBands,nWfs)		)
		allocate(	Ham(nWfs,nWfs)		)
		allocate(	tmp(nBands,nWfs))
		allocate(	EnP(nWfs,nQ)		)
		!
		if( debugProj ) then
			write(*,*)	"[projectUnk]: debug mode, will test unk normalization after debuging"
		end if
		!
		!

		if( doProjNUM) then
			call genTrialOrb(gnr)
			write(*,*)	"[projectUnk]: generated trial orbitals"
		end if
		!
		!PROJECT
		do qi = 1, nQ
			if( doProj ) then 
				Uq(:,:,qi)	= dcmplx(0.0_dp)
				call calcAmat(qi,ckH(:,:,qi) ,gnr, A) 
				call calcUmat(A, Uq(:,:,qi))
				!ROTATE BLOCH STATES
				!ckW(:,:,qi)	= matmul( ckH(:,:,qi) , Uq(:,:,qi)	 ) !			
				do n = 1, nWfs
					do m = 1, nWfs
						do gi = 1, nGq(qi)
							ckW(gi,n,qi)	= ckW(gi,n,qi) + ckH(gi,m,qi) * Uq(m,n,qi)
						end do
					end do
				end do
			else
				write(*,*)	"[projectUNK]: projection disabled, U= Identity"
				
				Uq(:,:,qi) = dcmplx(0.0_dp)
				do n = 1, size(Uq,1)
					do m = 1, size(Uq,2)
						if(n==m) then
							Uq(n,m,qi)	= dcmplx(1.0_dp)
						end if
					end do
				end do
				ckW(:,:,qi) = ckH(:,:,qi)
			end if
			!
			

			if( debugProj  ) then
				!if( .not. isUnit(ckW(:,:,qi))	)	write(*,*) "[projectUNK]: ckW not a unitary matrix at qi=",qi
				if( .not. isUnit(Uq(:,:,qi)) 	) 	write(*,*) "[projectUNK]: Uq not a unitary matrix at qi=",qi
			end if

		end do
		!

		write(*,*)	"[projectUnk]: done with projections at each k point"	
		!
		return
	end subroutine
		



	subroutine addThopp(qi, U,En, tHopp)
		!	<0|H|R> = (1/sqrt(N)) \sum_q Vq^dagger Eq Vq
		!	where Eq is diagonal matrix (energies in Hamiltonian gauge)
		integer,		intent(in)		:: qi
		complex(dp),	intent(in)		:: U(:,:)	!	U(nBands,nWfs)
		real(dp),		intent(in)		:: En(:,:)	!	En(nBands,nQ)
		complex(dp),	intent(inout)	:: tHopp(:,:,:)	!tHopp(nWfs, nWfs nSC)
		complex(dp)						:: phase
		complex(dp),	allocatable		:: tmp(:,:),EV(:,:), Udag(:,:)
		integer							:: R, m, n
		!
		allocate(	tmp(nBands,nWfs)	)
		allocate(	EV(nWfs,nWfs) 		)
		allocate(	Udag(nWfs,nBands)	)
		!		
		!
		do R = 1, nSC
			!tmp = Eq Vq
			do m = 1, nWfs
				do n = 1, nBands
					tmp(n,m)	= dcmplx(En(n,qi)) * U(n,m)
				end do 
			end do
			! EU = Vq^dagger tmp
			Udag	= dconjg( transpose(U)	)
			EV		= matmul(Udag, tmp)
			!ADD results
			phase			= myExp( - dot_product( qpts(:,qi), Rcell(:,R) ) 		)
			tHopp(:,:,R)	= tHopp(:,:,R) + phase * EV(:,:) /  dcmplx(real(nSC,dp))
		end do
		!
		!tHopp	= tHopp    
		!
		return
	end subroutine


	




!privat:
	subroutine genTrialOrb(gnr)
		complex(dp),	intent(out)	:: gnr(:,:)       !gnr( nR, nWfs)	
		integer						:: n, ri, gammaP
		!
		if(nWfs > nBands) then
			write(*,*)"[genTrialOrb]: warning, nWfs larger then nBands! No propper subspace..."
		end if
		!
		gammaP	= getGammaPoint()
		gnr 	= dcmplx(0.0_dp)
		!
		!SINGLE ATOM
		if( nAt == 1 ) then
			write(*,*)	"[genTrialOrb]: assuming a one atom per unit cell" 
			do n = 1, nWfs
				do ri = 1, nR
					if( 0.0_dp < rpts(1,ri) .and.  rpts(1,ri) < aX .and. 0.0_dp < rpts(2,ri) .and. rpts(2,ri) < aY) then
						if( 	insideAt(1, rpts(:,ri)))	then
							gnr(ri,n)	= potWell(1,n,ri)
						end if
					end if
				end do
			end do
		!TWO ATOMS
		else if( nAt == 2) then
			write(*,*)	"[genTrialOrb]: assuming a two atoms per unit cell" 
			do n = 1, nWfs-1, 2
				do ri = 1, nR
					if( 0.0_dp < rpts(1,ri) .and.  rpts(1,ri) < aX .and. 0.0_dp < rpts(2,ri) .and. rpts(2,ri) < aY) then
						if( 	insideAt(1, rpts(:,ri)))	then
							gnr(ri,n)	= potWell(1,n,ri)
						else if( insideAt(2, rpts(:,ri)) ) then
							gnr(ri,n+1) = potWell(2,n,ri)
						end if
					end if
				end do
			end do
		!DEFAULT		
		else
			write(*,*)	"[genTrialOrb]: to many atoms per unit cell for trial orbitals" 
		end if
		!
		return
	end subroutine



	complex(dp) function potWell(at, n, ri)
		integer, 	intent(in)		:: at, n, ri
		real(dp)					:: x,y, Lx, Ly, xc, yc, kx, ky, A
		!
		x	= rpts(1,ri)
		y	= rpts(2,ri)
		Lx 	= 2.0_dp * atR(1,at)
		Ly 	= 2.0_dp * atR(2,at)
		xc	= atPos(1,at)
		yc	= atPos(2,at)
		kx	= PI_dp/ Lx
		ky	= PI_dp/ Ly
		A	= dsqrt(4.0_dp / (Lx * Ly)	)
		!
		if( n==1 ) then
			potWell = A * dsin( kx*(x -xc + 0.5_dp * Lx) ) * dsin( ky*(y -yc + 0.5_dp * Ly) )
		else if( n==2 ) then
			potWell = A * dcos( kx*(x -xc + 0.5_dp * Lx) ) * dsin( ky*(y -yc + 0.5_dp * Ly) )
		else if( n==3 ) then
			potWell = A * dsin( kx*(x -xc + 0.5_dp * Lx) ) * dcos( ky*(y -yc + 0.5_dp * Ly) )
		end if 


		return
	end function







	subroutine calcAmat(qi,ckH ,gnr, A)
		integer,		intent(in)	:: qi
		complex(dp),	intent(in)	:: ckH(:,:), gnr(:,:) 
		complex(dp),	intent(out)	:: A(:,:) !A(nBands,nWfs)
		!calculates the inner product matrix of bwfs and trial orbital gnr
		if( doProjNUM ) then
			if( qi==1 )		write(*,*)"[calcAmat]: numeric calculation of projection matrix A"
			call calcAmatNUM(qi, ckH, gnr, A)
		else
			if( qi == 1)	write(*,*)"[calcAmat]: analytic calculation of projection matrix A"
			call calcAmatANA(qi, ckH, A)
		end if
		!
		!
		return
	end subroutine


	subroutine calcAmatANA(qi,ckH, A)
		!analytic projection with hard coded integrals
		!	projection onto sin**2, sin cos, cos sin
		integer,		intent(in)	:: qi
		complex(dp),	intent(in)	:: ckH(:,:)
		complex(dp),	intent(out)	:: A(:,:) !A(nBands,nWfs)
		integer						:: n, m
		!
		A	= dcmplx(0.0_dp)
		!reminder: use dconjg on ckH
		!
		!
		!SINGLE ATOM
		if( nAt == 1 ) then
			do n = 1, nWfs
				do m = 1, nBands
					select case( n )
					case(1)
						A(m,n)	= g1Int(qi,m,1, ckH)
					case(2)
						A(m,n)	= g2Int(qi,m,1, ckH)
					case(3)
						A(m,n)	= g3Int(qi,m,1, ckH)
					case default
						write(*,*)"[calcAmatANA]: Warning hit default in single atom switch"
						A(m,n)	= dcmplx(0.0_dp)
					end select
				end do
			end do
		!DUAL ATOMS
		else if( nAt == 2 ) then
			do n = 1, nWfs
				do m = 1, nBands
					select case( n )
					case(1)
						A(m,n)	= g1Int(qi,m,1,ckH)
					case(2)
						A(m,n)	= g1Int(qi,m,2,ckH)
					case(3)
						A(m,n)	= g2Int(qi,m,1,ckH)
					case(4)
						A(m,n)	= g2Int(qi,m,2,ckH)
					case(5)
						A(m,n)	= g3Int(qi,m,1,ckH)
					case(6)
						A(m,n)	= g3Int(qi,m,2,ckH)
					case default
						write(*,*)"[calcAmatANA]: Warning hit default in dual atom switch"
						A(m,n)	= dcmplx(0.0_dp)
					end select
				end do
			end do
		!FALLBACK
		else
			write(*,*)	"[calcAmatANA]: more then two atoms per unit cell. I can only handle two tough"
		end if 
		

		!
		return
	end subroutine


	complex(dp) function g1Int(qi,m, at,ckH)
		integer,		intent(in)	:: qi, m, at
		complex(dp),	intent(in)	:: ckH(:,:)
		complex(dp)					:: num1, num2, denom
		real(dp)					:: kappa, xL, xR, yL, yR, Gx, Gy, xc, yc
		integer						:: gi
		!
		!TRIAL ORBITAL:
		kappa	= PI_dp / ( 2.0_dp * atR(1,at) )
		if( atR(1,at) /= atR(2,at) ) write(*,*)"[g1Int]: warning analytic projection can not handle non cubic wells"
		xc		= atPos(1,at) - atR(1,at)
		yc		= atPos(2,at) - atR(2,at)
		xL 		= atPos(1,at) - atR(1,at) 
		xR		= atPos(1,at) + atR(1,at)
		yL		= atPos(2,at) - atR(2,at)
		yR		= atPos(2,at) + atR(2,at)
		!
		!SUMMATION OVER G:
		g1Int 	= dcmplx(0.0_dp)
		do gi = 1, nGq(qi)
			Gx 		= Gvec(1,gi,qi)
			Gy		= Gvec(2,gi,qi)
			!
			!
			if( 	abs(Gx) < machineP 	.and.	 abs(Gy) < machineP 	) then
				num1	= dcos((xL-xc)*kappa) - dcos((xR-xc)*kappa)
				num2	= dcos((yL-yc)*kappa) - dcos((yR-yc)*kappa)
				denom	= kappa**2
			!
			else if( abs(Gy) < machineP ) then
				num1	= 	myExp(-Gx*(xL+xR)) 	* 	( 			dcos((yL-yc)*kappa) - 				dcos((yR-yc)*kappa)	)
				num2	= 			myExp(Gx*xR)* 	( kappa * 	dcos((xL-xc)*kappa) + i_dp * Gx * 	dsin((xL-xc)*kappa)	)
				num2	= num2  - 	myExp(Gx*xL)* 	( kappa *	dcos((xR-xc)*kappa) + i_dp * Gx *	dsin((xR-xc)*kappa)	)
				denom	= kappa * ( kappa**2 - Gx**2)
			!
			else if( abs(Gx) < machineP) then
				num1	= 	myExp(-Gy*(yL+yR))	*	(			dcos((xL-xc)*kappa) -				dcos((xR-xc)*kappa)	)
				num2	=		-	myExp(Gy*yR)*	( kappa *	dcos((yL-yc)*kappa) + i_dp * Gy *	dsin((yL-yc)*kappa)	)
				num2	= num2  +	myExp(Gy*yL)*	( kappa *	dcos((yR-yc)*kappa) + i_dp * Gy *	dsin((yR-yc)*kappa)	)
				denom	= kappa * ( Gy**2 - kappa**2 )
			!
			else
				num1 	=    	- myExp(-Gx*xL)	* 	( kappa *	dcos((xL-xc)*kappa)	+ 	i_dp * Gx * dsin((xL-xc)*Kappa) )
				num1 	= num1  + myExp(-Gx*xR) * 	( kappa * 	dcos((xR-xc)*kappa)	+	i_dp * Gx * dsin((xR-xc)*Kappa)	)
				num2 	= 		  myExp(Gy*yL) 	* 	( kappa *	dcos((yR-yc)*kappa)	+	i_dp * Gy * dsin((yR-yc)*kappa)	) 
				num2 	= num2  - myExp(Gy*yR) 	* 	( kappa *	dcos((yL-yc)*kappa)	+	i_dp * Gy * dsin((yL-yc)*kappa)	)
				denom	= myExp(Gy*(yL+yR)) * (Gy**2-kappa**2) * (Gx**2-kappa**2)
			!
			end if
			!
			!
			if( abs(denom) < machineP) then
				write(*,'(a,i4,a,f8.4,a,f8.4a,f8.4)') "[g1Int]: warning zero denom at qi=",qi," Gx=",Gx,"Gy=",Gy,"denom=",denom
				denom 	= 1e-5_dp
			end if
			!
			!
			g1Int = g1Int + dconjg(ckH(gi,m)) * num1 * num2 / (dsqrt(vol) * denom)
		end do
		!
		!
		return
	end function


	complex(dp) function g2Int(qi,m, at,ckH)
		integer,		intent(in)	:: qi, m, at
		complex(dp),	intent(in)	:: ckH(:,:)
		complex(dp)					:: num1, num2, denom
		real(dp)					:: kappa, xL, xR, yL, yR, Gx, Gy, xc, yc
		integer						:: gi
		!
		!TRIAL ORBITAL:
		kappa	= PI_dp / (2.0_dp*atR(1,at))
		if( atR(1,at) /= atR(2,at) ) write(*,*)"[g2Int]: warning analytic projection can not handle non cubic wells"
		xc		= atPos(1,at) - atR(1,at)
		yc		= atPos(2,at) - atR(2,at)
		xL 		= atPos(1,at) - atR(1,at) 
		xR		= atPos(1,at) + atR(1,at)
		yL		= atPos(2,at) - 2.0_dp * atR(2,at)
		yR		= atPos(2,at) + 2.0_dp * atR(2,at)
		!
		!SUMMATION OVER G:
		g2Int 	= dcmplx(0.0_dp)
		do gi = 1, nGq(qi)
			!
			Gx 		= Gvec(1,gi,qi)
			Gy		= Gvec(2,gi,qi)
			!
			!
			if( abs(Gx) < machineP .and. abs(Gy) < machineP ) then
				num1	= -  ( 	dcos((xL-xc)*kappa) - dcos((xR-xc)*kappa) )
				num2	= 	   	dsin((yL-yc)*kappa) - dsin((yR-yc)*kappa)
				denom	= 	kappa**2
			!
			else if( abs(Gy) < machineP ) then
				num1	=		  myExp(-Gx*(xL+xR))*	(			dsin((yL-yc)*kappa) - 				dsin((yR-yc)*kappa)	)
				num2	= 		- myExp(-Gx*xR)		*  	( kappa *	dcos((xL-xc)*kappa) +	i_dp * Gx * dsin((xL-xc)*kappa)	)
				num2	= num2 	+ myExp(Gx*xL)		*	( kappa *	dcos((xR-xc)*kappa) +	i_dp * Gx * dsin((xR-xc)*kappa)	)
				denom	= kappa * ( kappa**2 - Gx**2 )
			!
			else if( abs(Gx) < machineP ) then
				num1	= 		  myExp(-Gy*(yL+yR))*	(			dcos((xL-xc)*kappa) -				dcos((xR-xc)*kappa)	)
				num2	=		  myExp(Gy*yR)		*	( kappa *	dsin((yL-yc)*kappa) -	i_dp * Gy * dcos((yL-yc)*kappa)	)
				num2	= num2 	+ myExp(Gy*yL)		*	(-kappa *	dsin((yR-yc)*kappa) +	i_dp * Gy * dcos((yR-yc)*kappa)	)
				denom	= kappa * ( Gy**2 - kappa**2 )
			!
			else  
				num1 	=		- myExp(-Gx*xL)		* 	( kappa*	dcos((xL-xc)*kappa) + 	i_dp * Gx * dsin((xL-xc)*Kappa)	)
				num1 	= num1	+ myExp(-Gx*xR)		* 	( kappa*	dcos((xR-xc)*kappa) +	i_dp * Gx * dsin((xR-xc)*Kappa)	)
				num2 	=		+ myExp(Gy*yR) 		* 	( kappa*	dsin((yL-yc)*kappa) -	i_dp * Gy * dcos((yL-yc)*kappa)	)
				num2 	= num2	+ myExp(Gy*yL) 		* 	(-kappa* 	dsin((yR-yc)*kappa) +	i_dp * Gy *	dcos((yR-yc)*kappa)	) 
				denom	=		+ myExp(Gy*(yL+yR)) * (Gy**2-kappa**2) * (Gx**2-kappa**2)  
			end if
			!
			!
			if( abs(denom) < machineP) then
				write(*,'(a,i4,a,f8.4,a,f8.4a,f8.4)') "[g2Int]: warning zero denom at qi=",qi," Gx=",Gx,"Gy=",Gy,"denom=",denom
				denom 	= 1e-5_dp
			end if
			!
			!
			g2Int = g2Int + dconjg(ckH(gi,m)) * num1 * num2 / (dsqrt(vol) * denom)
		end do
		!
		!
		return
	end function

	complex(dp) function g3Int(qi,m, at,ckH)
		integer,		intent(in)	:: qi, m, at
		complex(dp),	intent(in)	:: ckH(:,:)
		complex(dp)					:: num1, num2, denom
		real(dp)					:: kappa, xL, xR, yL, yR, Gx, Gy, xc, yc
		integer						:: gi
		!
		kappa	= PI_dp / (2.0_dp*atR(1,at))
		if( atR(1,at) /= atR(2,at) ) write(*,*)"[g3Int]: warning analytic projection can not handle non cubic wells"
		xc		= atPos(1,at) - atR(1,at)
		yc		= atPos(2,at) - atR(2,at)
		xL 		= atPos(1,at) - 2.0_dp * atR(1,at) 
		xR		= atPos(1,at) + 2.0_dp * atR(1,at)
		yL		= atPos(2,at) - atR(2,at)
		yR		= atPos(2,at) + atR(2,at)
		g3Int 	= dcmplx(0.0_dp)
		do gi = 1, nGq(qi)
			!
			Gx 		= Gvec(1,gi,qi)
			Gy		= Gvec(2,gi,qi)
			!
			!
			if( abs(Gx) < machineP .and. abs(Gy) < machineP ) then
				num1	= - ( 	dcos((yL-yc)*kappa) - dcos((yR-yc)*kappa) 	)
				num2	= 		dsin((xL-xc)*kappa) - dsin((xR-xc)*kappa)
				denom	= kappa**2
			!
			else if( abs(Gy) < machineP ) then
				num1	= 		  myExp(-Gx*(xL+xR))*	( 			dcos((yL-yc)*kappa) - 				dcos((yR-yc)*kappa) )
				num2	= 		  myExp(Gx*xR)		*	(-kappa* 	dsin((xL-xc)*kappa) + 	i_dp * Gx * dcos((xL-xc)*kappa)	)
				num2	= num2	+ myExp(Gx*xL)		*	( kappa*	dsin((xR-xc)*kappa) - 	i_dp * Gx * dcos((xR-xc)*kappa)	)
				denom	= kappa * ( kappa**2 - Gx**2)
			!
			else if( abs(Gx) < machineP ) then
				num1	= 		  myExp(-Gy*(yL+yR))*	(			dsin((xR-xc)*kappa) -				dsin((xL-xc)*kappa)	)
				num2	=		- myExp(Gy*yR)		*	( kappa*	dcos((yL-yc)*kappa) +	i_dp * Gy * dsin((yL-yc)*kappa)	)
				num2	= num2	+ myExp(Gy*yL)		*	( kappa*	dcos((yR-yc)*kappa) +	i_dp * Gy * dsin((yR-yc)*kappa)	)
				denom	= kappa * ( Gy**2 - kappa**2 )
			!		 		
			else
				!ToDo: revisit
				num1	=  		- i_dp * myExp(Gx*xR) * Gx 		* dcos((xL-xc)*kappa)
				num1	= num1 	+ i_dp * myExp(Gx*xL) * Gx 		* dcos((xR-xc)*kappa)
				num1	= num1 	+ 	    myExp(Gx*xR) * kappa 	* dsin((xL-xc)*kappa)
				num1	= num1 	- 		myExp(Gx*xL) * kappa	* dsin((xR-xc)*kappa)
				!
				num2	=  		- 		myExp(Gy*yR) *  ( kappa	* dcos((yL-yc)*kappa) + i_dp * Gy * dsin((yL-yc)*kappa) )
				num2	= num2 +		myExp(Gy*yL) * 	( kappa * dcos((yR-yc)*kappa) + i_dp * Gy * dsin((yR-yc)*kappa) )
				!
				denom	=  + myExp( Gx*(xL+xR) + Gy*(yL+yR) )	* (Gy**2-kappa**2) * (kappa**2-Gx**2)
				!  
			end if
			!
			!
			if( abs(denom) < machineP) then
				write(*,'(a,i4,a,f8.4,a,f8.4a,f8.4)') "[g3Int]: warning zero denom at qi=",qi," Gx=",Gx,"Gy=",Gy,"denom=",denom
				denom 	= 1e-5_dp
			end if
			!
			!
			g3Int = g3Int + dconjg(ckH(gi,m)) * num1 * num2 / (dsqrt(vol) * denom)
		end do
		!
		!
		return
	end function





	subroutine calcAmatNUM(qi, ckH, gnr, A)
		!numerical caluclation allows for quick change of trial orbitals 
		integer,		intent(in)	:: qi
		complex(dp),	intent(in)	:: ckH(:,:), gnr(:,:) 
		complex(dp),	intent(out)	:: A(:,:) !A(nBands,nWfs)
		complex(dp),	allocatable	:: psi(:), basVec(:)
		integer						:: m,n, xi, gMax
		!
		allocate(	basVec(nG)	)
		allocate(	psi(nBands)	)
		!
		gMax = nGq(qi)
		!
		!INTEGRATION OVER REAL SPACE
		A = dcmplx(0.0_dp)
		do xi = 1, nR
			!GET BASIS VECTOR
			call calcBasis(qi,xi, basVec)
			psi(:)	= matmul( basVec(1:gMax), ckH(1:gMax,1:nBands) ) / dsqrt(vol)
			psi(:)	= dconjg(psi)
			!CALC OVERLAP AT CURRENT GRID POINT
			do n = 1, nWfs
				do m = 1 , nBands 
					A(m,n)	= A(m,n) + psi(m) * gnr(xi,n)
				end do 
			end do
		end do
		!NORMALIZE REAL SPACE INTEGRATION
		A	= A / real(nR,dp)
		!
		!
		return
	end subroutine









	subroutine calcUmat(A, U)
		!U = A S^-0.5 = Z I V
		!	where A= Z d V (singular value decomposition) 
		complex(dp),	intent(in)		:: A(:,:)	!A(nBands,nWfs)
		complex(dp),	intent(out)		:: U(:,:)	!U(nBands,nWfs)
		complex(dp),	allocatable		:: Z(:,:), I(:,:), V(:,:), TMP(:,:)
		real(dp),		allocatable		:: d(:)
		integer							:: m, n, k, lda, ldb, ldc
		complex(dp)						:: alpha, beta
		character*1						:: transa, transb
		!
		allocate(	Z(	nBands, 	nBands	)	)
		allocate(	I(	nBands,		nWfs	)	)
		allocate(	V(	nWfs,		nWfs	)	)
		allocate(	TMP(nBands,		nWfs	)	)
		allocate(	d(	nWfs				)	)
		!
		!SET UP IDENTITY MATRIX
		I	= dcmplx(0.0_dp)
		do m = 1, nBands
			do n = 1, nWfs
				if( n == m) then
					I(m,n)	= dcmplx(1.0_dp)
				end if
			end do
		end do
		!
		!CALC SVD OF A
		call mySVD(A,Z,d,V)
		!
		!CALC U
		tmp=matmul(I,V)
		U=matmul(Z,tmp)

		!calc U with mkl
		!transa	= 'n'
		!transb	= 'n'
		!alpha	= dcmplx(1.0_dp)
		!beta	= dcmplx(0.0_dp)
		!m		= size(I,1)
		!n		= size(V,2)
		!k		= size(I,2)
		!lda		= m
		!ldb		= k
		!ldc		= m
		!call zgemm(transa, transb, m, n, k, alpha, I, lda, V, ldb, beta, TMP, ldc)
		!
		
		!transa	= 'n'
		!transb	= 'n'
		!alpha	= dcmplx(1.0_dp)
		!beta	= dcmplx(0.0_dp)
		!m		= size(Z,1)
		!n		= size(TMP,2)
		!k		= size(Z,2)
		!lda		= m
		!ldb		= k
		!ldc		= m
		!call zgemm(transa, transb, m, n, k, alpha, Z, lda, TMP, ldb, beta, U, ldc)
		!
		!
		return
	end subroutine









	!subroutine calcInvSmat(A, S, smin, smax)
	!	!calculates the sqrt inv. of the overlap matrix S
	!	!
	!	!	S	= A^dagger A
	!	complex(dp),	intent(in)		:: A(:,:)
	!	complex(dp),	intent(out)		:: S(:,:)
	!	real(dp),		intent(inout)	:: smin, smax
	!	complex(dp),	allocatable		:: Sold(:,:), Ssqr(:,:)
	!	integer							:: m,n,k,lda,ldb,ldc, i,j
	!	real(dp)						:: mi, ma
	!	complex(dp)						:: alpha, beta
	!	character*1						:: transa, transb
	!	!
	!	allocate(	Sold(nWfs,nWfs)	)
	!	!
	!	!
	!	!CALCULATE S FROM A 
	!	transa	= 'c'
	!	transb	= 'n'
	!	m		= size(A,1)
	!	n		= size(A,2)
	!	k		= size(A,2)
	!	alpha	= dcmplx(1.0_dp)
	!	beta	= dcmplx(0.0_dp)
	!	lda		= size(A,1)
	!	ldb		= size(A,1)
	!	ldc		= size(A,1)
	!	!call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
	!	 call zgemm(transa, transb, m, n, k, alpha, A, lda, A, ldb, beta, S, ldc)
	!	!
	!	!
	!	!CALCULATE INVERSE SQRT OF S
	!	Sold = S
!
!	!	call myMatInvSqrt(S, mi, ma)
!	!	!
!	!	if( mi < smin ) then
!	!		smin = mi
!	!	end if
!	!	if( ma > smax ) then
!	!		smax = ma
!	!	end if
!	!	!
!	!	!
!	!	!DEBUG
!	!	if(debugHam) then
!	!		! test if S * (S^-0.5)^2 == I
!	!		write(*,*)"[calcInvSmat]: all done, start debugging test"
!	!		allocate(	Ssqr(nWfs,nWfs)	)
!	!		transa	= 'n'
!	!		transb	= 'n'
!	!		alpha	= dcmplx(1.0_dp)
!	!		beta	= dcmplx(0.0_dp)
!	!		call zgemm(transa, transb, m, n, k, alpha, S   , lda, S   , ldb, beta, Ssqr, ldc)
!	!		call zgemm(transa, transb, m, n, k, alpha, Sold, lda, Ssqr, ldb, beta, Sold, ldc)
!	!		if( .not. isIdentity(Sold) ) then
!	!			write(*,*)"[calcInvSmat]: problem with mat inversion, seems to be not inverse square root"
!	!		end if
!	!	end if
!	!	!
!	!	!
!	!	return
!	!end subroutine
!

		


end module projection