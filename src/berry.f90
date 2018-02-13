module berry
	!module contains methods related to the Berry phase theory of modern polarization
	!	i.e. calculation of connection, velocities, curvatures and polarization
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP,  myExp, myLeviCivita, aUtoAngstrm, aUtoEv
	use sysPara
	
	use w90Interface,	only:	read_U_matrix, read_M_initial, readBandVelo, read_FD_scheme, read_wann_centers
	use basisIO,		only:	read_energies, read_velo
	use polarization,	only:	calcPolViaA
	use semiClassics,	only:	calcFirstOrdP
	use output,			only:	writePolFile, writeVeloHtxt, writeConnTxt

	implicit none

	private
	public	::					berryMethod



	integer	::					num_wann, num_kpts










	contains




!public
	subroutine berryMethod()
		complex(dp),	allocatable		:: 	U_mat(:,:,:), M_mat(:,:,:,:), M_wann(:,:,:,:), &
											Aconn_H(:,:,:,:), Aconn_W(:,:,:,:), FcurvQ(:,:,:,:),veloQ(:,:,:,:) 		
		real(dp),		allocatable		::	EnQ(:,:), b_k(:,:), w_b(:), &
											w_centers(:,:), berry_W_gauge(:,:),berry_H_gauge(:,:), niu_polF2(:,:), niu_polF3(:,:)
		integer							::	nntot
		integer,		allocatable		:: 	nnlist(:,:), nncell(:,:,:)
		!
		!READ FD SCHEME
		num_wann = nWfs
		num_kpts = nQ
		call read_FD_scheme(nntot, nnlist, nncell, b_k, w_b)
		!
		!k-space
		allocate(			EnQ(		nSolve	 							,	nQ		)			)
		allocate(			U_mat(	nWfs	,	nWfs						,	nQ		)			)
		allocate(			M_mat(	nWfs	, 	nWfs		, 	nntot 		,	nQ		)			)
		allocate(			M_wann(	nWfs	,	nWfs		,	nntot		,	nQ		)			)		
		allocate(			Aconn_H(		3		, 	nWfs	,	nWfs	,	nQ		)			)
		allocate(			Aconn_W(		3		, 	nWfs	,	nWfs	,	nQ		)			)		
		!
		!real-space
		allocate(			w_centers(		3,					nWfs					)			)
		allocate(			berry_W_gauge(	3,					nWfs					)			)
		allocate(			berry_H_gauge(	3,					nWfs					)			)
		allocate(			niu_polF2(		3,					nWfs					)			)
		allocate(			niu_polF3(		3,					nWfs					)			)
		!
		!
		call read_energies(EnQ)
		!read wannier files
		call read_M_initial(M_mat)
		call read_U_matrix(U_mat)
		call read_wann_centers(w_centers)
		!
		!HAM GAUGE
		write(*,*)	"[berryMethod]: start (H) gauge calculation"
		call calcConnOnCoarse(M_mat, nntot, b_k, w_b, Aconn_H)
		call calcPolViaA(Aconn_H, berry_H_gauge)
		!
		!WANN GAUGE
		write(*,*)	"[berryMethod]: start (W) gauge calculation"

		call rot_M_matrix(nntot, nnlist, M_mat, U_mat, M_wann)
		call calcConnOnCoarse(M_wann, nntot, b_k, w_b, Aconn_W)
		call calcPolViaA(Aconn_W, berry_W_gauge)
		!
		!1st ORDER SEMICLASSICS
		if(doNiu) then
			write(*,*)	"[berrryMethod]: now calc first order pol"
			allocate(	FcurvQ(	3,	nWfs, nWfs,		nQ)		)
			allocate(	veloQ(	3, 	nSolve,nSolve,	nQ)		)
			FcurvQ	= dcmplx(0.0_dp)	!does not matter since <FcurvQ,AconnQ> is always zero in 2D
			call calcVelo(U_mat , Aconn_W, EnQ ,  veloQ)
			call calcFirstOrdP(FcurvQ, Aconn_W, veloQ, EnQ, niu_polF2, niu_polF3)
		else
			niu_polF2 = 0.0_dp
			niu_polF3 = 0.0_dp
		end if
		!
		!OUTPUT
		call printNNinfo(nntot, nnlist, w_b)
		call writePolFile(w_centers, berry_H_gauge, berry_W_gauge, niu_polF2, niu_polF3)
		call writeConnTxt( Aconn_W )
		call writeVeloHtxt( veloQ )	!*aUtoEv*aUtoAngstrm )				
		!
		!
		return
	end subroutine










!private
	subroutine calcConnOnCoarse(M_mat, nntot, b_k, w_b, A_conn)
		!calculates the Berry connection, based on finite differences
		!nntot, nnlist, nncell list the nearest neighbours (k-space)
		!b_k:	 non zero if qi at boundary of bz and nn accross bz
		!w_b:	 weight of nn 
		complex(dp),	intent(in)			:: 	M_mat(:,:,:,:)
		integer,		intent(in)			::	nntot
		real(dp),		intent(in)			::	b_k(:,:), w_b(:)
		complex(dp),	intent(out)			::	A_conn(:,:,:,:)
		complex(dp)							::	delta
		integer								::	qi, nn, n, m
		!
		if( 		fastConnConv ) write(*,*)	"[calcConnOnCoarse]: use logarithmic formula for connection"
		if( .not.	fastConnConv ) write(*,*)	"[calcConnOnCoarse]: use finite difference formula for connection" 
		!
		!$OMP PARALLEL DO SCHEDULE(STATIC)	DEFAULT(SHARED) PRIVATE(qi, nn, n,m, delta)
		do qi = 1, nQ
			A_conn(:,:,:,qi) = dcmplx(0.0_dp)
			!SUM OVER NN
			do nn = 1, nntot
				!
				!WEIGHT OVERLAPS (Fast Convergence)
				if( fastConnConv ) then
					A_conn(1,:,:,qi)	= A_conn(1,:,:,qi)	+	w_b(nn) * b_k(1,nn) * dimag( log(M_mat(:,:,nn,qi))	)
					A_conn(2,:,:,qi)	= A_conn(2,:,:,qi)	+	w_b(nn) * b_k(2,nn) * dimag( log(M_mat(:,:,nn,qi))	)
					A_conn(3,:,:,qi)	= A_conn(3,:,:,qi)	+	w_b(nn) * b_k(3,nn) * dimag( log(M_mat(:,:,nn,qi))	)
				!WEIGHT OVERLAPS (Finite Difference)
				else
					do n = 1, nWfs
						do m = 1, nWfs
							delta = dcmplx(0.0_dp)
							if( n==m ) 	delta = dcmplx(1.0_dp)
							A_conn(1,m,n,qi)		= 	A_conn(1,m,n,qi)	+	w_b(nn) * b_k(1,nn) * dimag( M_mat(m,n,nn,qi) - delta  )					
							A_conn(2,m,n,qi)		= 	A_conn(2,m,n,qi)	+	w_b(nn) * b_k(2,nn) * dimag( M_mat(m,n,nn,qi) - delta  )
							A_conn(3,m,n,qi)		= 	A_conn(3,m,n,qi)	+	w_b(nn) * b_k(3,nn) * dimag( M_mat(m,n,nn,qi) - delta  )	
						end do
					end do				
				end if
			end do
		end do
		!$OMP END PARALLEL DO
		!
		return
	end subroutine


	subroutine applyRot(ck, Uq, ckW)
		complex(dp),	intent(in)		::	ck(:,:,:), Uq(:,:,:)
		complex(dp),	intent(out)		::	ckW(:,:,:)
		integer							::	qi, gi, n, m
		!
		ckW	= dcmplx(0.0_dp)
		!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(gi, n, m) ,&
		!$OMP& SCHEDULE(STATIC)
		do qi = 1, nQ
			do gi = 1, nGq(qi)			
				do n = 1, nSolve
					!BOUND STATES
					if( n <= nWfs ) then
						!
						!SUM_M
						do m = 1, num_wann
							ckW(gi,n,qi)	=  ckW(gi,n,qi) + Uq(m,n,qi)   * ck(gi,m,qi)	
						end do
					!CONDUCTING STATES
					else
						ckW(gi,n,qi)	= ck(gi,n,qi)
					end if
					!
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		write(*,*)	"[berry/applyRot]: applied U matrix to basis coefficients"

		!
		return
	end subroutine


	subroutine rot_M_matrix(nntot, nnlist, M_H, U_mat, M_W)
		!	M^(W)(nn,qi)	= U^*(qi)	M^(H)(nn,qi) U(qnb(nn)) 
		integer,			intent(in)		::	nntot, nnlist(:,:)
		complex(dp),		intent(in)		::	M_H(:,:,:,:), U_mat(:,:,:)
		complex(dp),		intent(out)		::	M_W(:,:,:,:)
		complex(dp),		allocatable		::	U_left(:,:), tmp(:,:)
		integer								:: 	qi, nn
		!
		allocate(		tmp(	size(U_mat,1), size(U_mat,2)	)			)
		allocate(		U_left(	size(U_mat,1), size(U_mat,2)	)			)
		!
		do qi = 1, nQ
			do nn = 1, nntot
				tmp(:,:)		=	matmul(		M_H(:,:,nn,qi)	,	U_mat(:,:,nnlist(qi,nn) )			)
				U_left(:,:)		=	dconjg(		 transpose( U_mat(:,:,qi) )			)
				!
				!
				M_W(:,:,nn,qi)	= matmul(	U_left(:,:), tmp(:,:) )
			end do
		end do
		!
		return
	end subroutine



	subroutine calcVelo(U_mat , A_mat, En_vec ,  v_mat)
		complex(dp),	intent(in)		::	U_mat(:,:,:), A_mat(:,:,:,:)
		real(dp),		intent(in)		::	En_vec(:,:)
		complex(dp),	intent(out)		::	v_mat(:,:,:,:)
		!
		!PLANE WAVE GRADIENT
		if( .not. doVeloNum ) then
			write(*,*)	"[beryMethod/calcVelo]: velo via plane wave gradients"
			call read_velo(v_mat)
		!BLOUNT
		else
			write(*,*)	"[beryMethod/calcVelo]: velo via blount formula - WARNING this is deprecated please set doVeloNum = f"
			call calcVeloBLOUNT(U_mat, A_mat, En_vec, v_mat)			
		end if
		!
		return
	end subroutine



	subroutine calcVeloBLOUNT(U_mat , A_mat, En_vec ,  v_mat)
		complex(dp),	intent(in)		::	U_mat(:,:,:), A_mat(:,:,:,:)
		real(dp),		intent(in)		::	En_vec(:,:)
		complex(dp),	intent(out)		::	v_mat(:,:,:,:)
		complex(dp),	allocatable		:: 	Abar(:,:,:), U(:,:), Ucjg(:,:), tmp(:,:)
		real(dp),		allocatable		::	v_Band(:,:,:)
		integer							::	n, m, qi, a

		allocate(			Abar(		3		,	nWfs	,	nWfs				)			)
		allocate(			tmp(					nWfs	,	nWfs				)			)
		allocate(			U(						nWfs	,	nWfs				)			)
		allocate(			Ucjg(					nWfs	,	nWfs				)			)
		allocate(			v_Band(		3		,			nWfs		,	nQ		)			)
		!
		v_mat	= dcmplx(0.0_dp)
		!
		call readBandVelo( v_Band )
		!
		do qi = 1, nQ
			!(H) GAUGE
			if( doGaugBack ) then
				!GET ROTATED QUANTITIES
				U	 = U_mat(:,:,qi)
				Ucjg = transpose( dconjg(U)	)
				do a = 1, 3
					tmp(:,:)	= matmul(	A_mat(a,:,:,qi) 	,	U	)
					Abar(a,:,:)	= matmul(	Ucjg				,	tmp	)
				end do
				!
				!APPLY
				do m = 1, nWfs
					do n = 1, nWfs
						if(n==m)	v_mat(1:3,n,n,qi)	= v_Band(1:3,n,qi)
						if(n/=m) 	v_mat(1:3,n,m,qi)	= - i_dp * dcmplx( En_vec(m,qi)-En_vec(n,qi) ) * Abar(1:3,n,m)
					end do
				end do
			!(W) GAUGE
			else
				do m = 1, nWfs
					do n = 1, nWfs
						if(n==m)	v_mat(1:3,n,n,qi)	= v_Band(1:3,n,qi)
						!if(n/=m) 	v_mat(1:3,n,m,qi)	= - i_dp * dcmplx( En_vec(m,qi)-En_vec(n,qi) ) * A_mat(1:3,n,m,qi)
						if(n/=m) 	v_mat(1:3,n,m,qi)	= - i_dp * dcmplx( En_vec(m,qi)-En_vec(n,qi) ) * A_mat(1:3,n,m,qi)
						!if(n/=m)	v_mat(1:3,n,m,qi)	= - i_dp * dcmplx( En_vec(m,qi)-En_vec(n,qi) ) * Abar(1:3,n,m)
					end do
				end do
			end if
		end do	

		!

		return
	end subroutine


	subroutine printNNinfo(nntot, nnlist, w_b)
		integer, 	intent(in)				::	nntot, nnlist(:,:)
		real(dp),	intent(in)				::	w_b(:)
		integer								::	gammaPt, nn
		real(dp)							::	nnEstimate
		!
		write(*,'(a,i2,a)')	"[berryMethod]: nn info:"
		gammaPt = 1 + int(	nQx*(0.5_dp+0.5_dp*nQy)	)
		write(*,'(a,f6.2,a,f6.2,a)')	"        dqx=",dqx,"; dqy=",dqy,"."
		write(*,'(a,f6.2,a,f6.2,a)')	"        this means for qpt=(",qpts(1,gammaPt),", ",qpts(2,gammaPt),")."
		write(*,*)	" nn  | q_nn(x) | q_nn(y) | w_b "
		write(*,*)	"-------------------------------"
		do nn = 1, nntot
			
			write(*,'(i2,a,f6.2,a,f6.2,a,f6.2)')	nn,"  |  ",	qpts(1,nnlist(gammaPt,nn)),"  |  ",&
																				qpts(2,nnlist(gammaPt,nn)),"  |  ",w_b(nn)
			
			nnEstimate = 2.0_dp / ( 	w_b(nn) * dqx**2  ) !wbx 	= 2.0_dp / 	( real(Z,dp) * dqx**2 ) ; where Z is number of nearest neighbours
			if( int(nnEstimate)/= nntot) write(*,*)	"[berryMethod]: weights suggest ",int(nnEstimate)," nearest neigbours, where ",nntot," are expected"
		end do
		!
		return
	end subroutine

end module berry