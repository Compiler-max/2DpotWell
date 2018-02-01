module berry
	!module contains methods related to the Berry phase theory of modern polarization
	!	i.e. calculation of connection, velocities, curvatures and polarization
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP,  myExp, myLeviCivita, aUtoAngstrm, aUtoEv
	use sysPara
	use w90Interface,	only:	read_U_matrix, readBandVelo, read_FD_scheme, read_wann_centers
	use planeWave,		only:	calcVeloGrad, calcConnOnCoarse
	use polarization,	only:	calcPolViaA
	use semiClassics,	only:	calcFirstOrdP
	use output,			only:	writePolFile, writeVeloHtxt, writeConnTxt

	implicit none

	private
	public	::	berryMethod



	integer								::	num_wann, num_kpts










	contains




!public
	subroutine berryMethod(ck, EnQ)
		complex(dp),	intent(in)		::	ck(:,:,:)
		real(dp),		intent(in)		::	EnQ(:,:)
		real(dp)						:: 	pBerry(3), pNiuF2(3), pNiuF3(3), pPei(3)
		real(dp),		allocatable		:: 	R_real(:,:)
		complex(dp),	allocatable		:: 	U_matrix(:,:,:),	ck_wann(:,:,:), AconnQ(:,:,:,:), FcurvQ(:,:,:,:),veloQ(:,:,:,:) 		
		real(dp)						::	pWann(3)
		real(dp),		allocatable		::	v_Band(:,:,:), krel(:,:), b_k(:,:), w_b(:), w_centers(:,:), b_centers(:,:), niu_centers(:,:)
		integer							::	nntot
		integer,		allocatable		:: 	nnlist(:,:), nncell(:,:,:)
		!					
		!COARSE
		allocate(			ck_wann(				GmaxGLOBAL	, 	nSolve	,  	nQ		)			)
		allocate(			AconnQ(			3		, 	nWfs	,	nWfs	,	nQ		)			)
		allocate(			FcurvQ(			3		,	nWfs	,	nWfs	,	nQ		)			)
		allocate(			veloQ(			3		, 	nSolve	,	nSolve	,	nQ		)			)
		allocate(			R_real(			3		,							nSC		)			)
		allocate(			v_Band(			3		,			nSolve		,	nQ		)			)
		allocate(			krel(			3,									nQ		)			)
		allocate(			U_matrix(					nWfs	,	nWfs	,	nQ		)			)
		allocate(			w_centers(		3,					nWfs					)			)
		allocate(			b_centers(		3,					nWfs					)			)
		allocate(			niu_centers(	3,					nWfs					)			)
		!
		!READ IN QUANTITIES
		num_wann = nWfs
		num_kpts = nQ
		call read_FD_scheme(nntot, nnlist, nncell, b_k, w_b)
		call read_wann_centers(w_centers)
		!
		!ROTATE
		call read_U_matrix(R_real, U_matrix, krel)
		call applyRot(ck, U_matrix, ck_wann)
		!
		!!CONNECTION (via K space)
		!call calcConnOnCoarse(ck_wann, nntot, nnlist, nncell, b_k, w_b, AconnQ)

		call calcPolViaA(AconnQ,pBerry)
		!write(*,*)"[berryMethod]: coarse rotated pol =(",pBerry(1),", ",pBerry(2),", ", pBerry(3),")."
		!!
		!!1st ORDER SEMICLASSICS
		if(doNiu) then
			write(*,*)	"[berrryMethod]: now calc first order pol"
			FcurvQ	= dcmplx(0.0_dp)	!does not matter since <FcurvQ,AconnQ> is always zero in 2D
			call calcVelo(ck , U_matrix , AconnQ, EnQ ,  veloQ)
			!
			call calcFirstOrdP(FcurvQ, AconnQ, veloQ, EnQ, pNiuF2, pNiuF3)
			write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[berryMethod]: pNiuF2=(",pNiuF2(1),", ",pNiuF2(2),", ",pNiuF2(3),")."
			write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[berryMethod]: pNiuF3=(",pNiuF3(1),", ",pNiuF3(2),", ",pNiuF3(3),")."
		else
			niu_centers = 0.0_dp
			pNiuF2 = 0.0_dp
			pNiuF3 = 0.0_dp
		end if
		!!
	
		!!
		!!OUTPUT
		call writePolFile(pWann, pBerry, pNiuF2, pNiuF3, pPei )
		call writeConnTxt( AconnQ )
		call writeVeloHtxt( veloQ)	!*aUtoEv*aUtoAngstrm )				
		!!
		!write(*,*)	"[berrryMethod]: all done"
		!
		!
		!
		return
	end subroutine










!private
	subroutine applyRot(ck, Uq, ckW)
		complex(dp),	intent(in)		::	ck(:,:,:), Uq(:,:,:)
		complex(dp),	intent(out)		::	ckw(:,:,:)
		integer							::	qi, gi, n, m
		!
		ckW	= dcmplx(0.0_dp)
		if(	useRot ) then
			do qi = 1, nQ
				do gi = 1, nGq(qi)			! u^(H) = u^(W) U -> u^(W) = u^(H) U^dagger
					do n = 1, num_wann
						!SUM OVER m
						do m = 1, num_wann
							ckW(gi,n,qi)	=  ckW(gi,n,qi) + Uq(m,n,qi)   * ck(gi,m,qi)	
						end do		
					end do	
					do n = num_wann+1, nSolve
							ckW(gi,n,qi)	= ck(gi,n,qi)
					end do
				end do
			end do
			write(*,*)	"[berry/applyRot]: applied U matrix to basis coefficients"
		else 
			ckW	= ck
			write(*,'(a,a)')	"[berry/applyRot]: rotations disabled.",&
									" Will use initial electronic structure coeff"
		end if
		!
		!
		return
	end subroutine



	subroutine calcVelo(ck, U_mat , A_mat, En_vec ,  v_mat)
		complex(dp),	intent(in)		::	ck(:,:,:), U_mat(:,:,:), A_mat(:,:,:,:)
		real(dp),		intent(in)		::	En_vec(:,:)
		complex(dp),	intent(out)		::	v_mat(:,:,:,:)
		complex(dp),	allocatable		:: 	Abar(:,:,:), U(:,:), Ucjg(:,:), tmp(:,:)
		real(dp),		allocatable		::	v_Band(:,:,:)
		integer							::	n, m, qi, a
		!
		if( .not. doVeloNum ) then
			!PLANE WAVE GRADIENT
			write(*,*)	"[beryMethod/calcVelo]: velo via plane wave gradients"
			call calcVeloGrad( ck, v_mat)
			!BLOUNT
		else
			write(*,*)	"[beryMethod/calcVelo]: velo via blount formula - warning this is deprecated please set doVeloNum = f"
			allocate(			Abar(		3		,	nWfs	,	nWfs				)			)
			allocate(			tmp(					nWfs	,	nWfs				)			)
			allocate(			U(						nWfs	,	nWfs				)			)
			allocate(			Ucjg(					nWfs	,	nWfs				)			)
			allocate(			v_Band(		3		,			nWfs		,	nQ		)			)
			!
			call readBandVelo( v_Band )
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
		end if
		!
		!
		return
	end subroutine



end module berry