module planeWave
	!generates bloch and lattice periodidc functions given a basCoeff matrix
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, acc, machineP,& 
									myExp, myLeviCivita
	use sysPara

	implicit none

	private
	public	::	calcVeloGrad, calcMmat, calcAmatANA, calcConnOnCoarse, calcBasis


	contains







!public
	subroutine calcVeloGrad(ck, v_mat )
		!calculates the velocity operator matrix
		!	Psi_n v Psi_m	= i/hbar Psi_n grad_r Psi_m
		!					= - 1 / hbar sum_G ckn^dag ckm G
		complex(dp),	intent(in)		:: 	ck(:,:,:)
		complex(dp),	intent(out)		::	v_mat(:,:,:,:)
		integer							::	qi, m, n, gi
		!
		v_mat = dcmplx(0.0_dp)
		!
		if(	size(ck,3)/=size(v_mat,4)	) stop	"[calcVeloGrad]: coeff and velo defined on different k meshes, stop now"
			
		!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(qi, m, n, gi)
		do qi = 1, nQ
			do m = 1, nSolve
				do n = 1, nSolve
					!
					!SUM OVER BASIS FUNCTIONS
					do gi = 1 , nGq(qi)
						v_mat(1:2,n,m,qi) = v_mat(1:2,n,m,qi) -  dconjg(ck(gi,n,qi)) *  ck(gi,m,qi) *  Gvec(1:2,gi,qi)
					end do
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		!
		return
	end subroutine


	subroutine calcMmat(qi,knb,gShift, nGq, Gvec, ck, Mmat)
		integer,		intent(in)		:: qi, knb, nGq(:)
		real(dp),		intent(in)		:: gShift(2),  Gvec(:,:,:)
		complex(dp),	intent(in)		:: ck(:,:,:)
		complex(dp),	intent(out)		:: Mmat(:,:)
		integer							:: gi, gj, n, m, cnt
		real(dp)						:: delta(2)
		logical							:: found
		!
		Mmat	= dcmplx(0.0_dp)
		cnt		= 0
		do gi = 1, nGq(qi)
			found 	= .false.
			gj			= 1
			!
			do while( gj<= nGq(knb) .and. (.not. found) ) 
				!find gj, which fullfills the delta condition
				delta(1:2)	=  ( Gvec(1:2,gi,qi)-qpts(1:2,qi) ) 	-  		( Gvec(1:2,gj,knb)-qpts(1:2,knb)-gShift(1:2) )
				!
				if( norm2(delta) < 1e-8_dp )	then
					do n = 1, size(Mmat,2)
						do m = 1, size(Mmat,1)
							Mmat(m,n)	= Mmat(m,n)	+ dconjg(	ck(gi,m,qi)	) * ck(gj,n,knb)
						end do
					end do
					cnt = cnt + 1
					found = .true.
				end if
				!
				gj = gj + 1
			end do
			!if( .not. found  ) write(*,'(a,i5,a,i5,a,f6.2,a,f6.2,a)')	"[calcMmat]: WARNING no matching Gvec found for qi=",qi," q_nn=",knb,&
			!															" gshift=(",gShift(1),",",gShift(2),")."
			!
			!
		end do
		!
		if( cnt /= nGq(qi)	)		write(*,'(a,i8,a,i8)')	"[calcMmat]: WARNING, found ",cnt," neighbouring Gvec, where nGmax(qi)=",nGq(qi)
		!
		!
		return
	end subroutine



	subroutine calcAmatANA(qi,ckH, A_matrix)
		!analytic projection with hard coded integrals
		!	projection onto sin**2, sin cos, cos sin
		integer,		intent(in)	:: 	qi
		complex(dp),	intent(in)	:: 	ckH(:,:)
		complex(dp),	intent(out)	:: 	A_matrix(:,:) !A(nBands,nWfs)
		real(dp)					::	r_state
		integer						:: 	n, m, perAtom, at, state
		!
		A_matrix	= dcmplx(0.0_dp)
		!
		perAtom	= nWfs / nAt
		if( qi == 1 )	write(*,'(a,i2,a,i2,a)')	"[calcAmatANA]: projecting onto ",nWfs," states, distributed over ",nAt," atoms"
		if( qi == 1 )	write(*,'(a,i2,a)')			"[calcAmatANA]: projecting onto ",perAtom," states per atom"
		!
		if( perAtom == 1 ) then
			do n = 1, nWfs
				do m = 1, nBands
					at	= n
					A_matrix(m,n)	= g1Int(qi,m, at, ckH)
				end do
			end do
		else if( perAtom == 3 ) then
			do n = 1, nWfs
				do m = 1, nBands
					!DETERMINE ATOM AND STATE TO PROJECT
					at 		= mod(n, nAt)
					if( at == 0) 	at = nAt
					r_state	= real(n,dp) / real(nAt,dp)
					if(	r_state <= 1) then
						state = 1
					else if( r_state <= 2) then
						state = 2
					else if( r_state <= 3) then
						state = 3
					else 
						state = 0
					write(*,*)	"[calcAmatANA]: WARNING! will set A_matrix component to zero, try to run wannier with use_bloch switch"
					end if
					!
					!DO PROJECTION
					if(	state == 1 ) A_matrix(m,n) 	= g1Int(qi,m, at, ckH)
					if( state == 2 ) A_matrix(m,n)	= g2Int(qi,m, at, ckH)
					if( state == 3 ) A_matrix(m,n)	= g3Int(qi,m, at, ckH)
					if( state == 0 ) A_matrix(m,n)	= dcmplx(0.0_dp)
				end do
			end do
		else
			write(*,*)	"[calcAmatANA]: only 1 or 3 states per atom supported. set A_mat to zero"
			A_matrix = dcmplx(0.0_dp)
		end if
		!
		!
		return
	end subroutine






	subroutine calcConnOnCoarse(ck, nntot, nnlist, nncell, b_k, w_b, A_conn)
		!calculates the Berry connection, based on finite differences
		!nntot, nnlist, nncell list the nearest neighbours (k-space)
		!b_k:	 non zero if qi at boundary of bz and nn accross bz
		!w_b:	 weight of nn 
		complex(dp),	intent(in)			:: 	ck(:,:,:)
		integer,		intent(in)			::	nntot, nnlist(:,:), nncell(:,:,:)
		real(dp),		intent(in)			::	b_k(:,:), w_b(:)
		complex(dp),	intent(out)			::	A_conn(:,:,:,:)
		complex(dp),	allocatable			::	M_matrix(:,:)
		real(dp)							::	gShift(2), gX, gY
		integer								::	qi, nn
		!
		A_conn = dcmplx(0.0_dp)
		allocate(	M_matrix( size(A_conn,2), size(A_conn,3) )			)
		!
		gX = 2.0_dp * PI_dp / aX
		gY = 2.0_dp * PI_dp / aY
		!
		!$OMP PARALLEL DO SCHEDULE(STATIC)	DEFAULT(SHARED) PRIVATE(qi, nn, gShift)
		do qi = 1, nQ
			do nn = 1, nntot
				!GET OVERLAP MATRIX
				gShift(1)	= real(nncell(1,qi,nn),dp) * gX
				gShift(2)	= real(nncell(2,qi,nn),dp) * gY
				call calcMmat(qi, nnlist(qi, nn), gshift, nGq, Gvec, ck, M_matrix)
				!
				!WEIGHT OVERLAPS (Fast Convergence)
				A_conn(1,:,:,qi)	= w_b(nn) * b_k(1,nn) * dimag( log(M_matrix(:,:))	)
				A_conn(2,:,:,qi)	= w_b(nn) * b_k(2,nn) * dimag( log(M_matrix(:,:))	)
				A_conn(3,:,:,qi)	= w_b(nn) * b_k(3,nn) * dimag( log(M_matrix(:,:))	)
				!
				!WEIGHT OVERLAPS (FD)
			end do
		end do
		!$OMP END PARALLEL DO
		!
		!DEBUG
		if( .not. B1condition(b_k, w_b) )	stop '[calcConnCoarse]: B1 condition (2D version) not fullfilled'
		!
		return
	end subroutine




	subroutine calcBasis(qi, ri, basVec)
		!calculates the basis vectors e^i(k+G).r
		!	if |k+G| is larger then the cutoff the basis vector is set to zero
		!	the cutoff enforces a symmetric base at each k point
		integer,	 intent(in)		:: qi, ri
		complex(dp), intent(out)	:: basVec(:)
		integer 				 	:: i 
		!
		basVec	= 0.0_dp
		do i =1, nGq(qi)
			basVec(i) 		= myExp( dot_product( Gvec(1:2,i,qi), rpts(1:2,ri) )		)  !/ dsqrt(vol)
		end do
		!
		return
	end subroutine







!prviat
	logical function B1condition(b_k, w_b)
		! test if
		!		sum_b{w_b * b_a * b_b}	= \delta_ab 
		! is true for all a,b
		real(dp),		intent(in)		:: 	b_k(:,:), w_b(:)
		logical							:: 	xx, xy, yy
		!
		xx	= .false.
		xy	= .false.
		yy	= .false.
		!
		if(		abs(	sum(w_b(:)*b_k(1,:)*b_k(1,:)) - 1.0_dp	)		< 0.1_dp				) 			xx = .true.
		if(		abs(	sum(w_b(:)*b_k(2,:)*b_k(2,:)) - 1.0_dp	)		< 0.1_dp				) 			yy = .true.
		if(		abs(	sum(w_b(:)*b_k(1,:)*b_k(2,:)) 			)		< 0.1_dp				)			xy = .true.
		!
		B1condition = xx .and. yy .and. xy
		!
		return
	end function


		complex(dp) function g1Int(qi,m, at,ckH)
		integer,		intent(in)	:: qi, m, at
		complex(dp),	intent(in)	:: ckH(:,:)
		complex(dp)					:: num1, num2, denom
		real(dp)					:: kappa, xL, xR, yL, yR, Gx, Gy, xc, yc
		integer						:: gi
		!
		!TRIAL ORBITAL:
		kappa	= PI_dp / ( 2.0_dp * atR(1,at) )
		if( atR(1,at) /= atR(2,at) )	stop	"[g1Int]: ERROR analytic projection can not handle non cubic wells"
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
				write(*,'(a,i4,a,f8.4,a,f8.4a,f8.4)') "[g1Int]: WARNING zero denom at qi=",qi," Gx=",Gx,"Gy=",Gy,"denom=",denom
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
		if( atR(1,at) /= atR(2,at) ) 	stop "[g2Int]: WARNING analytic projection can not handle non cubic wells"
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
				write(*,'(a,i4,a,f8.4,a,f8.4a,f8.4)') "[g2Int]: WARNING zero denom at qi=",qi," Gx=",Gx,"Gy=",Gy,"denom=",denom
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
		if( atR(1,at) /= atR(2,at) ) stop	"[g3Int]: WARNING analytic projection can not handle non cubic wells"
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
				write(*,'(a,i4,a,f8.4,a,f8.4a,f8.4)') "[g3Int]: WARNING zero denom at qi=",qi," Gx=",Gx,"Gy=",Gy,"denom=",denom
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






end module planeWave 

























!	subroutine calcConnOnCoarseQUADRATIC(ck, A)
!		!	deprecated, use calcConnOnCoarse
!		!finite difference on lattice periodic unk to calculate the Berry connection A
!		!	A_n(k) 	= <u_n(k)|i \nabla_k|u_n(k)>
!		!		 	= i  <u_n(k)| \sum_b{ w_b * b * [u_n(k+b)-u_n(k)]}
!		!			= i \sum_b{		w_b	 * [  <u_n(k)|u_n(k+b)> -  <u_n(k)|u_n(k)>]		}
!		!
!		! see Mazari, Vanderbilt PRB.56.12847 (1997), Appendix B
!		!
!		complex(dp),	intent(in)		:: ck(:,:,:) 	! ckW(nG, nWfs, nQ)		
!		complex(dp),	intent(out)		:: A(:,:,:,:)			
!		complex(dp),	allocatable		:: Mtmp(:,:)
!		integer							:: n, m, Z, qi, qx, qy, qxl, qxr, qyl, qyr, al, be
!		real(dp)						:: wbx,wby, bxl(2), bxr(2), byl(2), byr(2), delta, &
!											Gxl(2), Gyl(2), Gxr(2), Gyr(2)
!		!
!		if( size(nGq) 		/= nQ ) write(*,*)"[#",myID,";calcConnOnCoarse]: critical WARNING: basis array nGq has wrong size"
!		if(	size(Gvec,3)	/= nQ ) write(*,*)"[#",myID,";calcConnOnCoarse]: critical WARNING: basis array Gvec has wrong size"
!
!		A 		= dcmplx(0.0_dp)
!		Z 		= 4	!amount of nearest neighbours( 2 for 2D cubic unit cell)
!		wbx 	= 2.0_dp / 		( real(Z,dp) * dqx**2 )
!		wby		= wbx
!		write(*,*)	"[calcConnOnCoarse]	weight wb=", wbx
!		!wby 	= 1.0_dp /		( real(Z,dp) * dqy**2 )
!		!b vector two nearest X neighbours:
!		bxl(1) 	= -dqx				
!		bxl(2)	= 0.0_dp
!		bxr(1) 	= +dqx
!		bxr(2)	= 0.0_dp
!		!b vector two nearest Y neighbours:
!		byl(1) 	= 0.0_dp
!		byl(2)	= -dqy
!		byr(1) 	= 0.0_dp
!		byr(2)	= +dqy
!		!
!		!DEBUG WEIGHTS
!		do al = 1, 2
!			do be = 1, 2
!				delta 	= 0.0_dp
!				delta = delta + wbx * bxl(al) * bxl(be)
!				delta = delta + wbx * bxr(al) * bxr(be)
!				delta = delta + wby * byl(al) * byl(be)
!				delta = delta + wby * byr(al) * byr(be)
!				if( al==be .and. abs(delta-1.0_dp) > acc ) then
!					write(*,'(a,i1,a,i1,a,f6.3)') &
!							"[calcConnCoarse]: weights dont fullfill condition for a=",al," b=",be," delta=",delta
!				else if ( al/=be .and. abs(delta) > acc ) then
!					write(*,'(a,i1,a,i1,a,f6.3)') & 
!							"[calcConnCoarse]: weights dont fullfill condition for a=",al," b=",be,"delta=",delta
!				end if
!			end do
!		end do
!		!
!		write(*,'(a,f6.3,a,f6.3)')	"[calcConnOnCoarse]: dqx=",dqx," dqy=",dqy
!		!
!		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,n,qx,qy, qxl, qxr, qyl, qyr, qi, Gxl, Gxr, Gyl, Gyr, Mtmp)
!		allocate(	Mtmp( size(A,2), size(A,3) )		)
!		!$OMP DO COLLAPSE(2)  SCHEDULE(STATIC)
!		do qx = 1, nQx
!			do qy = 1, nQy
!				!GET NEIGHBOURS
!				qxl	= getLeft( qx,nQx)
!				qxr	= getRight(qx,nQx)
!				qyl	= getLeft(  qy,nQy)
!				qyr = getRight( qy,nQy)
!				!
!				!GET GRID POSITION OF NEIGHBOURS
!				qi	= getKindex(qx,qy)
!				qxl	= getKindex(qxl,qy)
!				qxr	= getKindex(qxr,qy)
!				qyl	= getKindex(qx,qyl)
!				qyr	= getKindex(qx,qyr)
!				!introduce test for neighbours
!				!
!				!SHIFT NEIGHBOURS BACK TO FIRST BZ
!				Gxl(:)	= 0.0_dp
!				Gxr(:)	= 0.0_dp
!				Gyl(:)	= 0.0_dp
!				Gyr(:)	= 0.0_dp
!				if( qx == 1 ) 	Gxl(1)	= - 2.0_dp * PI_dp / aX
!				if( qx == nQx)	Gxr(1)	= + 2.0_dp * PI_dp / aX
!				if( qy == 1 ) 	Gyl(2)	= - 2.0_dp * PI_dp / aY
!				if( qy == nQy)	Gyr(2)	= + 2.0_dp * PI_dp / aY
!				!
!				!UNCOMMENT FOR DEBUGGING
!				!write(*,*)"*"
!				!write(*,*)"*"
!				!write(*,*)"*"
!				!write(*,'(a,f6.3,a,f6.3,a)')	"[calcConnOnCoarse]: q_i=(",qpts(1,qi) ,", ",qpts(2,qi) ,")"
!				!write(*,'(a,f6.3,a,f6.3,a)')	"[calcConnOnCoarse]: qxl=(",qpts(1,qxl),", ",qpts(2,qxl),")"
!				!write(*,'(a,f6.3,a,f6.3,a)')	"[calcConnOnCoarse]: qxr=(",qpts(1,qxr),", ",qpts(2,qxr),")"
!				!write(*,'(a,f6.3,a,f6.3,a)')	"[calcConnOnCoarse]: qyl=(",qpts(1,qyl),", ",qpts(2,qyl),")"
!				!write(*,'(a,f6.3,a,f6.3,a)')	"[calcConnOnCoarse]: qyr=(",qpts(1,qyr),", ",qpts(2,qyr),")"
!				!write(*,*)"*"
!				!write(*,'(a,f6.3,a,f6.3,a)')"[calcConnOnCoarse]:  Gxl=",Gxl(1),", ",Gxl(2),")."
!				!write(*,'(a,f6.3,a,f6.3,a)')"[calcConnOnCoarse]:  Gxr=",Gxr(1),", ",Gxr(2),")."
!				!write(*,'(a,f6.3,a,f6.3,a)')"[calcConnOnCoarse]:  Gyl=",Gyl(1),", ",Gyl(2),")."
!				!write(*,'(a,f6.3,a,f6.3,a)')"[calcConnOnCoarse]:  Gyr=",Gyr(1),", ",Gyr(2),")."
!				!
!				!OVLERAPS:
!				!XL
!				!call calcMmat(qi, nnlist(qi,nn), gShift, nGq_glob, Gvec_glob, ck_glob, M_loc(:,:,nn,qi))
!				call calcMmat(qi, qxl, Gxl, nGq, Gvec, ck, Mtmp)
!				do n = 1, size(A,3)
!					do m = 1, size(A,2)
!						A(1:2,m,n,qi)	= A(1:2,m,n,qi) - wbx * bxl(1:2) * dimag( 	log(	Mtmp(m,n) )	 )
!					end do
!				end do
!				!XR
!				call calcMmat(qi, qxr, Gxr, nGq, Gvec, ck, Mtmp)
!				do n = 1, size(A,3)
!					do m = 1, size(A,2)
!						A(1:2,m,n,qi)	= A(1:2,m,n,qi) - wbx * bxr(1:2) * dimag( 	log(	Mtmp(m,n) )	 )
!					end do
!				end do
!				!YL
!				call calcMmat(qi, qyl, Gyl, nGq, Gvec, ck, Mtmp)
!				do n = 1, size(A,3)
!					do m = 1, size(A,2)
!						A(1:2,m,n,qi)	= A(1:2,m,n,qi) - wby * byl(1:2) * dimag( 	log(	Mtmp(m,n) )	 )
!					end do
!				end do
!				!YR
!				call calcMmat(qi, qyr, Gyr, nGq, Gvec, ck, Mtmp)
!				do n = 1, size(A,3)
!					do m = 1, size(A,2)
!						A(1:2,m,n,qi)	= A(1:2,m,n,qi) - wby * byr(1:2) * dimag( 	log(	Mtmp(m,n) )	 )
!					end do
!				end do
!				!
!			end do
!		end do
!		!$OMP END DO
!		!$OMP END PARALLEL
!		!
!		!
!		return
!	end subroutine