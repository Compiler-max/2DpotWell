module projection
	!module contains all the functions to generate the wannier functions
	!	this includes the calculation of the orthonormalized projected bloch wavefunctions
	!	and the FT of the bwfs to calculate the wannier functions
	use omp_lib
	use mathematics, 	only: 	dp, PI_dp, i_dp, acc, machineP, myExp
	use sysPara
	use planeWave,		only:	calcBasis
	implicit none	
	
	private
	public ::	calcAmatANA		














	contains
!public:
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
						write(*,*)	"[calcAmatANA]: error while determining the state to project on. n=",n
						write(*,*)	"[calcAmatANA]: error! will set A_matrix component to zero"
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



!privat:
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


end module projection