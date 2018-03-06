module pol_Berry
	!module contains methods related to the Berry phase theory of modern polarization
	!	i.e. calculation of connection, velocities, curvatures and polarization
	use omp_lib
	use util_sysPara,	only:	Bext, prefactF3, &
								nWfs, nQ, nSolve, vol, &
								qpts, aY, vol,  &
								atPos, atPot, &
								doGaugBack, doNiu, fastConnConv, doVeloNum 
	
	use util_w90Interf,	only:	read_U_matrix, read_M_initial, read_FD_b_vectors, read_wann_centers, &
								readBandVelo
	use util_basisIO,	only:	read_energies, read_velo
	use util_output,	only:	writePolFile, writeVeloHtxt, writeConnTxt, writeEnTXT, readEnTXT

	use pol_Niu,		only:	calcFirstOrdP
	
	implicit none

	private
	public	::					berryMethod

	!cloned from mathematics.f90:
	integer, 		parameter 	:: 	dp 				= kind(0.d0)
	real(dp), 		parameter	:: 	machineP 		= 1e-15_dp
	real(dp),		parameter	::	FDtol			= 1e-8_dp
	real(dp), 		parameter 	::	PI_dp 			= 4 * datan (1.0_dp)
	complex(dp),	parameter 	::	i_dp 			= dcmplx(0.0_dp, 1.0_dp)
	real(dp),		parameter 	::	aUtoEv	 		= 27.211385_dp
	real(dp),		parameter	::	aUtoAngstrm 	= 0.52917721092_dp
	real(dp),		parameter 	::	elemCharge	 	= 1.6021766208 * 1e-19_dp  *1e+6_dp! mu Coulomb

	!read in parameters:
	integer						::	num_wann, num_bands, num_kpts, num_stat, &
									nntot
	integer,		allocatable	::	nnlist(:,:), nncell(:,:,:)
	real(dp),		allocatable	::	b_k(:,:), w_b(:)








	contains




!public
	subroutine berryMethod()
		complex(dp),	allocatable		:: 	U_mat(:,:,:), M_ham(:,:,:,:), M_wann(:,:,:,:), &
											FcurvQ(:,:,:,:),veloQ(:,:,:,:) 		
		real(dp),		allocatable		::	Aconn_H(:,:,:,:), Aconn_W(:,:,:,:), &
											EnQ(:,:), EnClone(:,:), b_k(:,:), w_b(:), &
											w_centers(:,:), berry_W_gauge(:,:),berry_H_gauge(:,:), niu_polF2(:,:), niu_polF3(:,:)
		real(dp)						::	polQuantum, centiMet
		integer							::	n, qi
		!
		num_wann = nWfs		
		num_kpts = nQ
		num_stat = nSolve


		polQuantum 	= elemCharge / ( vol*aUtoAngstrm**2 ) 
		centiMet	= 1e+8_dp
		
		!READ W90 Files
		call read_M_initial(num_bands, num_kpts, nntot, nnlist, nncell, M_ham)
		call read_U_matrix(num_wann, U_mat)
		
		!WARNINGS & INFO
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,*)	"[berryMethod]: detected parameter info:"
		write(*,*)	"	num_kpts  =",num_kpts
		write(*,*)	"	num_bands =",num_bands
		write(*,*)	"	num_wann  =",num_wann
		write(*,*)	"	nntot  	  =",nntot
		if(	num_wann /= num_bands ) 	stop		"[berryMethod]: detected different num_wann(from _u.mat) and num_bands( .mmn file)"
		if(	nntot /= 4 )				write(*,*)	"[berryMethod]: WARNING nntot is not equal 4"
		write(*,*)	"*"

		!fd-scheme allo
		allocate(			w_b(		nntot	)	)
		allocate(			b_k(	3,	nntot	)	)
		!k-space allo
		allocate(			M_wann(		num_wann	,	num_wann	,	nntot		,	num_kpts		)			)		
		allocate(			Aconn_H(		3		, 	num_wann	,	num_wann	,	num_kpts		)			)
		allocate(			Aconn_W(		3		, 	num_wann	,	num_wann	,	num_kpts		)			)
		allocate(			EnQ(				num_stat							,	num_kpts		)			)
		!real-space
		allocate(			w_centers(		3,					num_wann					)			)
		allocate(			berry_W_gauge(	3,					num_wann					)			)
		allocate(			berry_H_gauge(	3,					num_wann					)			)
		allocate(			niu_polF2(		3,					num_wann					)			)
		allocate(			niu_polF3(		3,					num_wann					)			)
		!
		!
		!print unit cell info
		write(*,'(a,f8.3,a)')	"[berryMethod]: unit cell volume=",vol*aUtoAngstrm**2,"[Å^2] "


		!print atoms
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,*)		"[berryMethod]: atom positions:"
		write(*,*)		"	at | centers [Å] | V [eV]"
		do n = 1, size(atPos,2)
				write(*,'(i3,a,f16.8,a,f16.8,a ,f16.8)')	n," | ",atPos(1,n)*aUtoAngstrm,", ",atPos(2,n)*aUtoAngstrm," 	| ",atPot(n)*aUtoEv
		end do
	
		!print w90 centers
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,*)		"[berryMethod]: w90 centers:"
		call read_wann_centers(w_centers)
		write(*,*)		" #wf | 	<r>[Å]	"
		do n = 1, size(w_centers,2)
			write(*,'(i3,a,f16.8,a,f16.8,a,f16.8)')	n," | ",w_centers(1,n),", ",w_centers(2,n),", ",w_centers(3,n)
		end do




		!0th HAM GAUGE
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,*)		"[berryMethod]: start (H) gauge calculation"
		call read_FD_b_vectors(b_k, w_b)
		call calcConnOnCoarse(M_ham, Aconn_H)
		call calcPolViaA(polQuantum, centiMet, Aconn_H, berry_H_gauge)
		!
		!0th WANN GAUGE
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,*)	"*"
		write(*,*)		"[berryMethod]: start (W) gauge calculation"
		call rot_M_matrix(M_ham, U_mat, M_wann)
		call calcConnOnCoarse(M_wann, Aconn_W)
		call calcPolViaA(polQuantum, centiMet, Aconn_W, berry_W_gauge)
		!
		

		!
		!1st ORDER SEMICLASSICS
		if(doNiu) then
			!
			write(*,*)	"*"
			write(*,*)	"*"
			write(*,*)	"[berryMethod]: **************SEMICLASSICS*************************"
		
			allocate(	FcurvQ(	3,	num_wann, num_wann,		num_kpts)		)
			allocate(	veloQ(	3, 	num_stat, num_stat,		num_kpts)		)

			!read energies
			call readEnTXT(EnQ)


			!get velocities
			if(	doVeloNum) then
				if( .not. 	doGaugBack )	call calcVeloBLOUNT(Aconn_W, EnQ, veloQ)	
				if(			doGaugBack )	call calcVeloBLOUNT(Aconn_H, EnQ, veloQ)	
			else
				call read_velo(veloQ)
			end if
			!get curvature (toDo) 
			call calcCurv(FcurvQ)
			!
			!semiclassics
			call calcFirstOrdP(polQuantum, centiMet, Bext, prefactF3, FcurvQ, Aconn_W, veloQ, EnQ, niu_polF2, niu_polF3)
			!
			write(*,*)	"[berryMethod]: done with semiclassics"
		else
			niu_polF2 = 0.0_dp
			niu_polF3 = 0.0_dp
			write(*,*)	"[berryMethod]: semiclassics disabled, first order pol. set to zero"
		end if
		write(*,*)		"*"
		write(*,*)		"*"
		!
		!OUTPUT
		call writePolFile(polQuantum, centiMet, w_centers, berry_H_gauge, berry_W_gauge, niu_polF2, niu_polF3)
		write(*,*)		"[berryMethod]: wrote pol file"
		
		
		!TEST READING THE ENERGY FILE
		allocate(	EnClone(size(EnQ,1),size(EnQ,2)))
		call readEnTXT(EnClone)

		do qi = 1, size(EnQ,2)
			do n = 1, size(EnQ,1)
				if( abs(EnQ(n,qi)-EnClone(n,qi))> 1e-8_dp ) write(*,'(a,i3,a,i5)')	"[berryMethod]: error while reading en txt file at n=",n," qi=",qi
			end do
		end do

		call writeConnTxt( Aconn_W )
		!call writeVeloHtxt( veloQ )		
		!write(*,*)	"[berryMethod]: wrote k-space info files (connection & velocities)"			
		!
		!
		return
	end subroutine










!private

	!w90 subroutine :::   (rave are the centers)
!	csheet = dcmplx(1.0_dp)
!	sheet = 0.0_dp
 !   do nkp = 1, num_kpts
 !      do nn = 1, nntot
 !         do n = 1, num_wann
 !            ! Note that this ln_tmp is defined differently wrt the one in wann_omega
 !            ln_tmp(n,nn,nkp)=wb(nn)*( aimag(log(csheet(n,nn,nkp) * m_matrix(n,n,nn,nkp))) - sheet(n,nn,nkp) )
 !         end do
 !     end do
 !   end do
!
! !   ! recalculate rave
! !   rave = 0.0_dp
! !   do iw = 1, num_wann  
! !      do ind = 1, 3  
! !         do nkp = 1, num_kpts  
! !            do nn = 1, nntot  
! !               rave(ind,iw) = rave(ind,iw) +  bk(ind,nn,nkp) &
! !                    * ln_tmp(iw,nn,nkp)
! !            enddo
! !         enddo
! !      enddo
! !   enddo
 !   rave = -rave/real(num_kpts,dp)




	subroutine calcConnOnCoarse(M_mat, A_conn)
		!calculates the Berry connection, from the overlap matrix
		! see Vanderbilt 1997, eq.(22) & eq.(31)
		!b_k:	 non zero if qi at boundary of bz and nn accross bz
		!w_b:	 weight of nn 
		!
		!fastConnvergence:
		!			\vec{r_n} = 	1/N_k 	\sum_{k,b}	 w_b \vec{b} IMAG[ log M(n,n,nn,qi) ]
		!
		!finiteDifference:
		!			\vec{r_n} = 	1/N_k 	\sum_{k,b}	 w_b \vec{b} i_dp [ M(n,n,nn,qi) - 1.0_dp ]

		complex(dp),	intent(in)			:: 	M_mat(:,:,:,:)
		real(dp),		intent(out)			::	A_conn(:,:,:,:)
		real(dp),		allocatable			::	ln_tmp(:,:)
		real(dp)							::	my_bk(3,nntot), my_wb(nntot), bk_abs, wb_test
		complex(dp)							::	delta
		integer								::	qi, nn, n, m
		!
		if( 		fastConnConv ) write(*,*)	"[calcConnOnCoarse]: use logarithmic formula for connection"
		if( .not.	fastConnConv ) write(*,*)	"[calcConnOnCoarse]: use finite difference formula for connection" 
		if(num_bands /= num_wann ) stop			"[calcConnONCoarse]: WARNING disentanglement not supported"
		!
		A_conn 	= 0.0_dp

		!get FD weights:
		bk_abs	= norm2(	(qpts(1:2,nnlist(1,1)) + real(nncell(1:2,1,1)*2.0_dp * PI_dp * aY / vol) -qpts(1:2,1)) / aUtoAngstrm 	)
		wb_test	= 3.0_dp / ( 6.0_dp * bk_abs**2)
		!
		write(*,*)	"[calcConnOnCoarse]: wb=",wb_test," (ang^2)"
		write(*,*)	"[calcConnOnCoarse]: |b|=",bk_abs," (ang^-1)"		


		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(qi, nn, n,m, delta, ln_tmp, my_bk, my_wb  ) FIRSTPRIVATE(bk_abs, wb_test)
		allocate(		ln_tmp( size(M_mat,1), size(M_mat,2) )			)
		!$OMP DO SCHEDULE(STATIC)	
		do qi = 1, num_kpts
			!SUM OVER NN
			my_bk	= 0.0_dp
			do nn = 1, nntot
				!
				!GET WEIGHTS
				my_bk(1:2,nn)	= 	(	(qpts(1:2, nnlist(qi,nn)) + real(nncell(1:2,qi,nn),dp)* 2.0_dp * PI_dp * aY / vol	) - qpts(1:2,qi) ) / aUtoAngstrm
				if(	 ( bk_abs - norm2(my_bk(:,nn))	)	> 1e-8_dp	) then
					write(*,*)	"*"
					write(*,*)	"my_bk =",my_bk(:,nn)
					write(*,*)	"b_k=",b_k(:,nn)
				end if
				my_wb(nn)		= 3.0_dp / 	(real(nntot+2,dp)*norm2(my_bk(:,nn))**2)

				if( abs(my_wb(nn) - wb_test) > FDtol)	write(*,*)	"[calcConnOnCoarse]: warning, issue with weights detected" 
				!
				!(Fast Convergence)
				if( fastConnConv ) then
					ln_tmp(:,:)			= dimag(	log(	M_mat(:,:,nn,qi)	)		)
					A_conn(1,:,:,qi)	= A_conn(1,:,:,qi)	-	my_wb(nn) * my_bk(1,nn) * ln_tmp(:,:)
					A_conn(2,:,:,qi)	= A_conn(2,:,:,qi)	-	my_wb(nn) * my_bk(2,nn) * ln_tmp(:,:)
					A_conn(3,:,:,qi)	= A_conn(3,:,:,qi)	-	my_wb(nn) * my_bk(3,nn) * ln_tmp(:,:)
				!(Finite Difference)
				else
					do n = 1, size(M_mat,2)
						do m = 1, size(M_mat,1)
							delta = dcmplx(0.0_dp)
							if( n==m ) 	delta = dcmplx(1.0_dp)
							A_conn(1,m,n,qi)		= 	A_conn(1,m,n,qi)	-	my_wb(nn) * my_bk(1,nn) * dreal( i_dp * (M_mat(m,n,nn,qi) - delta)  )					
							A_conn(2,m,n,qi)		= 	A_conn(2,m,n,qi)	-	my_wb(nn) * my_bk(2,nn) * dreal( i_dp * (M_mat(m,n,nn,qi) - delta)  )
							A_conn(3,m,n,qi)		= 	A_conn(3,m,n,qi)	-	my_wb(nn) * my_bk(3,nn) * dreal( i_dp * (M_mat(m,n,nn,qi) - delta)  )	
						end do
					end do				
				end if
			end do
			!
			if( .not. FDcheck(my_wb, my_bk) ) then
				write(*,*)	"[calcConnOnCoarse]: problems with FD scheme at qi=,",qi
				do nn = 1, nntot
					write(*,'(a,i3,a,f10.4,a,f10.4,a,f10.4,a,e14.5)')	"b_k(nn=",nn,")=	(",my_bk(1,nn),", ",my_bk(2,nn),", ",my_bk(3,nn),"); wb=",my_wb(nn)
				end do
				!stop	'[calcConnOnCoarse]: violated B1 condition'
			end if
			!
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		write(*,*)	"[calcConnOnCoarse]: established connection"
		!
		return
	end subroutine


	logical function FDcheck(wb, b_k)
		!look at vanderbilt1997 (B1)
		real(dp),	intent(in)		::	wb(:), b_k(:,:)
		real(dp)					::	B1xy, B1xx, B1yy
		logical						::	xy, xx, yy
		!
		xy	= .false. 
		xx 	= .false.
		yy	= .false.
		!
		B1xy = sum( wb(:)*b_k(1,:)*b_k(2,:)	)
		B1xx = sum( wb(:)*b_k(1,:)*b_k(1,:)	)
		B1yy = sum( wb(:)*b_k(2,:)*b_k(2,:)	)
		!
		if( abs(B1xy-0.0_dp) < FDtol) 		xy = .true.
		if( abs(B1xx-1.0_dp) < FDtol)		xx = .true.
		if( abs(B1yy-1.0_dp) < FDtol)		yy = .true.
		!
		if(	.not. xy ) 	write(*,*)	"[FDcheck]: B1_xy=", B1xy," (expected 0.0)"
		if(	.not. xx ) 	write(*,*)	"[FDcheck]: B1_xx=", B1xx," (expected 1.0)"
		if(	.not. yy ) 	write(*,*)	"[FDcheck]: B1_yy=", B1yy," (expected 1.0)"
		!
		FDcheck	=	(xy .and. xx .and. yy)
		!
		!
	end function







	subroutine calcPolViaA(polQuantum, centiMet, A_mat, centers)
		!calculates the polarization by integrating connection over the brillouin zone
		! r_n 	= <0n|r|0n> 
		!		=V/(2pi)**2 \integrate_BZ <unk|i \nabla_k|unk>
		!		=V/(2pi)**2 \integrate_BZ A(k)
		real(dp),			intent(in)		:: 	polQuantum, centiMet, A_mat(:,:,:,:)			!A(2,	 nWfs, nWfs, nQ	)	
		real(dp),			intent(out)		:: 	centers(:,:)
		integer								::	n
		!
	
		write(*,'(a,e12.4,a)')	"[calcPolViaA]: the pol Quantum is p_quant=",polQuantum,"	[mu C/ Å²]"
		!
		centers = 0.0_dp
		!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(n)
		do n 	= 1, size(A_mat,2)
			!
			!INTEGRATE
			centers(1,n) = sum(A_mat(1,n,n,:)) / real(size(A_mat,4),dp)
			centers(2,n) = sum(A_mat(2,n,n,:)) / real(size(A_mat,4),dp)
			if(size(A_mat,1)==3)	centers(3,n) =  sum(A_mat(3,n,n,:)) / real(size(A_mat,4),dp)
			!
		end do
		!$OMP END PARALLEL DO
		!
		write(*,*)		" #state | 	<r>[Å]			| 	p[	\{mu}C/cm	]"
		do n = 1, size(A_mat,2)
			write(*,'(i3,a,f16.8,a,f16.8,a,f16.8,a,a,e13.4,a,e13.4,a)') n,"  | ",centers(1,n),", ", centers(2,n),",",centers(3,n), "  | ",&
																"(",centers(1,n)*polQuantum*centiMet,", ", centers(2,n)*polQuantum*centiMet, ")."
		end do
		write(*,'(a,e13.4,a,e13.4,a)')	"sum | 				(",sum(centers(1,:))*polQuantum*centiMet,", ",sum(centers(2,:))*polQuantum*centiMet,")."
		!		
		return
	end subroutine





	subroutine rot_M_matrix(M_H, U_mat, M_W)
		!	M^(W)(nn,qi)	= U^*(qi)	M^(H)(nn,qi) U(qnb(nn)) 
		complex(dp),		intent(in)		::	M_H(:,:,:,:), U_mat(:,:,:)
		complex(dp),		intent(out)		::	M_W(:,:,:,:)
		complex(dp),		allocatable		::	U_left(:,:), tmp(:,:)
		integer								:: 	qi, nn, qnn
		!
		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tmp, U_left, nn, qi)
		allocate(		tmp(	size(U_mat,1), size(U_mat,2)	)			)
		allocate(		U_left(	size(U_mat,1), size(U_mat,2)	)			)
		!$OMP DO SCHEDULE(STATIC)
		do qi = 1, size(M_H,4)
			U_left(:,:)			=	dconjg(		 transpose( U_mat(:,:,qi) )			)
			!
			do nn = 1, nntot
				qnn				=	nnlist(qi,nn)
				tmp(:,:)		=	matmul(		U_left(:,:)	,	M_H(:,:,nn,qi)		)
				M_W(:,:,nn,qi)	= 	matmul(		tmp(:,:)	, 	U_mat(:,:,qnn) 		)
			end do
			!
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		return
	end subroutine



	subroutine calcCurv(Fcurv_mat)
		!todo: implement this in case of 3d system
		complex(dp),		intent(out)		:: Fcurv_mat(:,:,:,:)
		!
		write(*,*)	"[calcCurv]: WARNING - calcCurv is not implemented (irrelevant for 2D system)"
		!
		Fcurv_mat = dcmplx(0.0_dp)
		!
		!
		return
	end subroutine



	subroutine calcVeloBLOUNT(A_conn, En_vec ,  v_mat)
		!use Blount 1962 formu.la
		! 
		real(dp),		intent(in)		::	A_conn(:,:,:,:)
		real(dp),		intent(in)		::	En_vec(:,:)
		complex(dp),	intent(out)		::	v_mat(:,:,:,:)
		real(dp),		allocatable		::	v_Band(:,:,:)
		integer							::	n, m, qi
		!
		allocate(	v_Band(	3,	num_wann,	num_kpts)		)
		!
		v_mat	= dcmplx(0.0_dp)
		!
		!DEBUG
		if( size(A_conn,1) /= size(v_mat,1)	)	write(*,*)	"[calcVeloBLOUNT]: A_conn and v_mat have different real space dimensionality"
		if( size(A_conn,2) /= size(v_mat,2) )	stop		"[calcVeloBLOUNT]: A_conn and v_mat have different amount of states covered"
		if( size(A_conn,3) /= size(v_mat,3) )	stop		"[calcVeloBLOUNT]: A_conn and v_mat have different amount of states covered"
		if( size(A_conn,4) /= size(v_mat,4) )	stop		"[calcVeloBLOUNT]: A_conn and v_mat live on different k meshes"
		!
		!GET DIAGONAL ELEMENTS
		call readBandVelo( v_Band ) !=band derivative
		!FILL MATRIX
		do qi = 1, size(A_conn,4)
			do m = 1, size(v_mat,3)
				do n = 1, size(v_mat,2)
					if(n==m)	v_mat(1:3,n,n,qi)	= dcmplx(v_Band(1:3,n,qi) )
					if(n/=m) 	v_mat(1:3,n,m,qi)	= dcmplx(		0.0_dp,		-1.0_dp * (En_vec(m,qi)-En_vec(n,qi)) * A_conn(1:3,n,m,qi)		 )
				end do
			end do
		end do	
		!
		!
		return
	end subroutine




end module pol_Berry

