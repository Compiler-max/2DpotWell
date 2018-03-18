module ham_PWbasis
	!generates bloch and lattice periodidc functions given a basCoeff matrix
	use omp_lib
	use util_math,	only:	dp, PI_dp,i_dp, acc, machineP,& 
									myExp
	use util_sysPara

	implicit none

	private
	public	::	calcVeloGrad, calcMmat, calcAmatANA, calcBasis, writeUNKs


	contains







!public
	subroutine calcVeloGrad(qi, ck, v_mat )
		!calculates the velocity operator matrix
		!	Psi_n v Psi_m	= i/hbar Psi_n grad_r Psi_m
		!					= - 1 / hbar sum_G ckn^dag ckm G
		integer,		intent(in)		::	qi
		complex(dp),	intent(in)		:: 	ck(:,:)
		complex(dp),	intent(out)		::	v_mat(:,:,:)
		integer							::	m, n, gi
		!
		v_mat 	= dcmplx(0.0_dp)
		!
		!get momentum
		do m = 1, nSolve
			do n = 1, nSolve
				!
				!SUM OVER BASIS FUNCTIONS
				do gi = 1 , nGq(qi)
					v_mat(1:2,n,m) = v_mat(1:2,n,m) -  dconjg(ck(gi,n)) *  ck(gi,m)  *  Gvec(1:2,gi,qi)
				end do
			end do
		end do
		!
		!
		return
	end subroutine




	subroutine calcMmat(qi, q_nn, gShift, nG_qi, nG_nn, Gvec_qi, Gvec_nn, ck_qi, ck_nn, Mmn	)
		integer,		intent(in)		::	qi, q_nn, nG_qi, nG_nn
		real(dp),		intent(in)		::	gShift(2), Gvec_qi(:,:), Gvec_nn(:,:)
		complex(dp),	intent(in)		::	ck_qi(:,:), ck_nn(:,:)
		complex(dp),	intent(out)		::	Mmn(:,:)
		integer							::	gi, gj, cnt, m, n, dG
		real(dp)						::	delta(2)
		logical							::	found
		!
		Mmn	= dcmplx(0.0_dp)
		cnt		= 0
		do gi = 1, nG_qi
			found 	= .false.
			gj			= 1
			!
			do while( gj<= nG_nn .and. (.not. found) ) 
				!find gj, which fullfills the delta condition
				delta(1:2)	=  ( Gvec_qi(1:2,gi)-qpts(1:2,qi) ) 	-  		( Gvec_nn(1:2,gj)-qpts(1:2,q_nn)-gShift(1:2) )
				!
				if( norm2(delta) < 1e-8_dp )	then
					do n = 1, size(Mmn,2)
						do m = 1, size(Mmn,1)
							Mmn(m,n)	= Mmn(m,n)	+ dconjg(	ck_qi(gi,m)	) * ck_nn(gj,n)
						end do
					end do
					cnt = cnt + 1
					found = .true.
				end if
				!
				gj = gj + 1
			end do
			!
		end do
		!
		dG = GmaxGLOBAL-GminGLOBAL
		if( cnt < GminGLOBAL-2*dG	)		write(*,'(a,i3,a,i8,a,i8)')	"[#",myID,"calcMmat]: WARNING, found ",cnt," neighbouring Gvec, where nGmax(qi)=",nG_qi
		!
		!		
		return
	end subroutine


	subroutine calcAmatANA(qi,ckH, A_matrix)
		!analytic projection with hard coded integrals
		!	projection onto sin**2, sin cos, cos sin
		!
		!	Amn = <Psi_m|g_n>
		integer,		intent(in)	:: 	qi
		complex(dp),	intent(in)	:: 	ckH(:,:)
		complex(dp),	intent(out)	:: 	A_matrix(:,:) !A(nBands,nWfs)
		integer						:: 	nWf, m
		!
		A_matrix	= dcmplx(0.0_dp)
		!
		do nWf = 1, nWfs
			if( proj_at(nWf) > nAt ) then
				write(*,'(a,i3,a,i2,a,i4,a)')	"[#",myID,"calcAmatANA]: can not project onto state atom",proj_at(nWf)," (unit cell contains: ",nAt," atoms)"
				write(*,'(a,i3,a,i2,a,i4,a)')	"[#",myID,"calcAmatANA]: projected atom out of bounds, will use unitary matrix..."
				A_matrix(nWf,nWf)	= dcmplx(1.0_dp)
			else
				do m = 1, nBands
					A_matrix(m,nWf)	= infiniteWell(qi, m, ckH, proj_at(nWf), proj_nX(nWf), proj_nY(nWf))
					!if(	proj_stat(nWf) == 1)	A_matrix(m,nWf)	= g1Int(qi,m, proj_at(nWf), ckH)
					!if( proj_stat(nWf) == 2)	A_matrix(m,nWf) = g2Int(qi,m, proj_at(nWf), cKH)
					!if( proj_stat(nWf) == 3)	A_matrix(m,nWf) = g3Int(qi,m, proj_at(nWf), cKH)					
					!if( proj_stat(nWf) > 3 )	then
					!	write(*,'(a,i3,a,i2,a)')	"[#",myID,"calcAmatANA]: can not project onto state #",proj_stat(nWf)," (maximum supported: 3)"
					!	A_matrix(m,nWf) = gDefaultInt(m, nWf)
					!end if
				end do
			end if
		end do
		!
		!
		return
	end subroutine



	subroutine writeUNKs(qi, nG_qi, ck, Gvec)
		integer,						intent(in)		::	qi, nG_qi
		complex(dp),					intent(in)		::	ck(:,:)
		real(dp),						intent(in)		::	Gvec(:,:)
		complex(dp)										::	tmp
		integer											::	ix, loop_b, nbnd, ig, spin
		character(len=1024)								::	unkFORM='(a,i5.5,a,i1)'
		character(len=15)								::	wfnname
		!
		!filename
		spin	= 1
		nbnd	= size(ck,2)
		write(wfnname, unkFORM) 'UNK',qi,'.', spin

		open(unit=205, file=w90_dir//wfnname,	form='formatted', access='stream', action='write', status='replace')
		!write to file
		write(205,*)	nRx, nRy, nRz, qi, nbnd
		do loop_b = 1, nbnd
			do ix = 1, nR


				!SUM G-VEC
				tmp = dcmplx(0.0_dp)
				do ig = 1, nG_qi
					tmp	= tmp & 
						+ ck(ig, loop_b) 	* 	myExp( 		dot_product( Gvec(1:2,ig)-qpts(1:2,qi),	rpts(1:2,ix) )		)  / dcmplx(sqrt(vol))
				end do
				!write result
				write(205,*) dreal(tmp)," ",dimag(tmp)
			end do
		end do
		close(205)
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



































!prviat:
		complex(dp) function infiniteWell(qi,m, ckH, at, nx, ny)
		!	calculates the integral
		!
		!		sin(nx*k*x ) * sin( ny*k*y) 	
		!
		!
		!
		integer,		intent(in)	:: qi, m, at, nx, ny
		complex(dp),	intent(in)	:: ckH(:,:)
		complex(dp)					:: num1, num2, denom
		real(dp)					:: kappaX, kappaY, xL, xR, yL, yR, Gx, Gy, xc, yc, at_rad
		integer						:: gi
		!
		!TRIAL ORBITAL:
		kappaX	= real(nx,dp) * PI_dp / ( 2.0_dp * atR(1,at) )
		kappaY 	= real(ny,dp) * PI_dp / ( 2.0_dp * atR(2,at) )
		if( atR(1,at) /= atR(2,at) )	stop	"[infiniteWell]: ERROR analytic projection can not handle non cubic wells"
		at_rad 	= atR(1,at)
		xc		= atPos(1,at)
		yc		= atPos(2,at)
		xL 		= atPos(1,at) - atR(1,at) 
		xR		= atPos(1,at) + atR(1,at)
		yL		= atPos(2,at) - atR(2,at)
		yR		= atPos(2,at) + atR(2,at)
		!
		!SUMMATION OVER G:
		infiniteWell 	= dcmplx(0.0_dp)
		do gi = 1, nGq(qi)
			Gx 		= Gvec(1,gi,qi)
			Gy		= Gvec(2,gi,qi)
			!
			!
			!
			!
			!CASE 1 (Gx==0 && Gy==0)
			if( 	abs(Gx) < machineP 	.and.	 abs(Gy) < machineP 	) then
				num1	= dcos((xL-xc+at_rad)*kappaX) - dcos((xR-xc+at_rad)*kappaX)
				num2	= dcos((yL-yc+at_rad)*kappaY) - dcos((yR-yc+at_rad)*kappaY)
				denom	= at_rad * kappaX*kappaY
			!----------------------------------------------------------------------------------------------------------------------------------------------
			!
			!	ToDo:
			!CASE 2 (Gx!=0	&& Gy==0)
			else if( abs(Gy) < machineP ) then
				num1	= 	myExp(-Gx*(xL+xR)) 	* 	( 				dcos((yL-yc+at_rad)*kappaY) - 				dcos((yR-yc+at_rad)*kappaY)	)
				num2	= 		-	myExp(Gx*xR)* 	( 	kappaX * 	dcos((xL-xc+at_rad)*kappaX) + i_dp * Gx * 	dsin((xL-xc+at_rad)*kappaX)	)
				num2	= num2  + 	myExp(Gx*xL)* 	( 	kappaX *	dcos((xR-xc+at_rad)*kappaX) + i_dp * Gx *	dsin((xR-xc+at_rad)*kappaX)	)
				denom	= at_rad * kappaY * ( Gx**2 - kappaX**2 )
			!----------------------------------------------------------------------------------------------------------------------------------------------
			!
			!
			!CASE 3 (Gx==0	&& Gy!=0)
			else if( abs(Gx) < machineP) then
				num1	= 	myExp(-Gy*(yL+yR))	*	(				dcos((xL-xc+at_rad)*kappaX) -				dcos((xR-xc+at_rad)*kappaX)	)
				num2	=		-	myExp(Gy*yR)*	( 	kappaY *	dcos((yL-yc+at_rad)*kappaY) + i_dp * Gy *	dsin((yL-yc+at_rad)*kappaY)	)
				num2	= num2  +	myExp(Gy*yL)*	( 	kappaY *	dcos((yR-yc+at_rad)*kappaY) + i_dp * Gy *	dsin((yR-yc+at_rad)*kappaY)	)
				denom	= at_rad * kappaX * ( Gy**2 - kappaY**2 )
			!----------------------------------------------------------------------------------------------------------------------------------------------
			!
			!
			!CASE 4 (Gx!=0	&& Gy!=0)
			else
				num1 	=    	- myExp(-Gx*xL)	* 	( kappaX *	dcos((xL-xc+at_rad)*kappaX)	+ 	i_dp * Gx * dsin((xL-xc+at_rad)*KappaX) )
				num1 	= num1  + myExp(-Gx*xR) * 	( kappaX * 	dcos((xR-xc+at_rad)*kappaX)	+	i_dp * Gx * dsin((xR-xc+at_rad)*KappaX)	)
				num2 	= 		  myExp(Gy*yL) 	* 	( kappaY *	dcos((yR-yc+at_rad)*kappaY)	+	i_dp * Gy * dsin((yR-yc+at_rad)*kappaY)	) 
				num2 	= num2  - myExp(Gy*yR) 	* 	( kappaY *	dcos((yL-yc+at_rad)*kappaY)	+	i_dp * Gy * dsin((yL-yc+at_rad)*kappaY)	)
				denom	= at_rad * myExp(Gy*(yL+yR)) * (Gy**2-kappaY**2) * (Gx**2-kappaX**2)
			!
			end if
			!----------------------------------------------------------------------------------------------------------------------------------------------
			!
			!DEBUG
			if( abs(denom) < machineP) then
				write(*,'(a,i4,a,f8.4,a,f8.4a,f8.4)') "[infiniteWell]: WARNING zero denom at qi=",qi," Gx=",Gx,"Gy=",Gy,"denom=",denom
				denom 	= 1e-8_dp
			end if
			!
			!
			!SUM OVER gi
			infiniteWell = infiniteWell + dconjg(ckH(gi,m)) * num1 * num2 / (dsqrt(vol) * denom)
		end do
		!
		!
		return
	end function





end module ham_PWbasis