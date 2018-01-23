module wannInterp
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP, myExp, myLeviCivita, &
								eigSolverFULL, rotMat, myCommutat, isHermitian, isUnit
	use sysPara,		only:	kpts, pw90GaugeB

	!use 
	implicit none

	private
	public ::					DoWannInterpol

	integer					::	num_wann


	contains

!public
	subroutine DoWannInterpol(ki, rHopp, tHopp, R_real, EnH, U_int, AconnH, FcurvH, veloH)
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: rHopp(:,:,:,:), tHopp(:,:,:)
		real(dp),		intent(in)		:: R_real(:,:)
		complex(dp),	intent(out)		:: U_int(:,:)
		real(dp),		intent(out)		:: EnH(:)
		complex(dp),	intent(out)		:: AconnH(:,:,:), FcurvH(:,:,:), veloH(:,:,:)
		complex(dp),	allocatable		:: HW(:,:), HaW(:,:,:), AW(:,:,:), FWtens(:,:,:,:), FWmat(:,:,:)
		integer							:: a, b, c
		!
		num_wann	= size(tHopp,1) 
		allocate(	HW(				num_wann, 	num_wann			)		)
		allocate(	HaW(	3	,	num_wann, 	num_wann			)		)
		allocate(	AW(		3	,	num_wann, 	num_wann			)		)
		allocate(	FWtens(	3,3	,	num_wann,	num_wann			)		)
		allocate(	FWmat( 	3	,	num_wann,	num_wann			)		)
		!
		EnH	= 0.0_dp
		AconnH	= dcmplx(0.0_dp)
		FcurvH	= dcmplx(0.0_dp)
		veloH	= dcmplx(0.0_dp)
		!
		!
		!GET K SPACE QUANTITIES		(EnH, U_int, HW, HaW, AW, FWtens)
		call wannInterpolator(ki, tHopp, rHopp, R_real, EnH, U_int(:,:), HW, HaW, AW, FWtens)
		!
		!WORK IN K SPACE
		if( pw90GaugeB ) then
			if(ki == 1) write(*,*)	"[DoWannInterpol]: start gauging back" 	
			call gaugeBack(ki, Hw, HaW, AW, FWtens, EnH, U_int(:,:), AconnH(:,:,:), FcurvH(:,:,:), veloH)	
		else
			if(ki ==1)	write(*,*)	"[DoWannInterpol]: Gauge trafo DISABLED	"
			!CONNECTION
			AconnH(1:3,:,:) 		= AW(1:3,:,:)
			!VELOCITIES
			call calcVeloNew(ki, EnH, U_int(:,:), HaW, AW, veloH)
			!
			!CURVATURE TO MATRIX
			do c = 1, 3
				do b = 1, 3
					do a = 1,3
						FcurvH(c,:,:)	= myLeviCivita(a,b,c) * FWtens(a,b,:,:)
					end do
				end do
			end do
		
		end if
		!
		!
		return
	end subroutine
























!privat
	subroutine wannInterpolator(ki, H_tb,r_tb, R_real, En_vec, U_mat, H_mat, Ha_mat, A_mat,Om_tens)
		integer,		intent(in)		::	ki
		complex(dp),	intent(in)		::	H_tb(:,:,:), r_tb(:,:,:,:)
		real(dp),		intent(in)		:: 	R_real(:,:)
		real(dp),		intent(out)		::	En_vec(:)
		complex(dp),	intent(out)		::	U_mat(:,:), H_mat(:,:), Ha_mat(:,:,:), A_mat(:,:,:), Om_tens(:,:,:,:)
		integer							:: R, a, b
		complex(dp)						:: phase
		!
		H_mat	= dcmplx(0.0_dp)
		Ha_mat	= dcmplx(0.0_dp)
		A_mat	= dcmplx(0.0_dp)
		Om_tens = dcmplx(0.0_dp)
		!
		!SET UP K SPACE MATRICES
		do R = 1, size(R_real,2)
			phase				= myExp( 	dot_product(kpts(1:2,ki),R_real(1:2,R))		) !/ dcmplx(real(nSC,dp))
			!
			H_mat(:,:)			= H_mat(:,:)	 	+ 			phase 								* H_tb(:,:,R)
			do a = 1, 3
				Ha_mat(a,:,:)	= Ha_mat(a,:,:) 	+ 			phase * i_dp * dcmplx(R_real(a,R))	* H_tb(:,:,R)
				A_mat(a,:,:)	= A_mat(a,:,:)		+ 			phase								* r_tb(a,:,:,R)
				!
				do b = 1, 3
					Om_tens(a,b,:,:)	= Om_tens(a,b,:,:)	+  	phase * i_dp * dcmplx(R_real(a,R)) 	* r_tb(b,:,:,R)
					Om_tens(a,b,:,:)	= Om_tens(a,b,:,:)	-  	phase * i_dp * dcmplx(R_real(b,R)) 	* r_tb(a,:,:,R)
				end do
			end do
		end do
		!
		!ENERGY INTERPOLATION
		U_mat(:,:)	= H_mat(:,:)
		if( .not. isHermitian(U_mat)	)		 	write(*,*)	"[wannInterpolator]: warning Ham is not hermitian"
		call eigSolverFULL(U_mat(:,:),	En_vec(:))
		U_mat	= transpose( dconjg(U_mat))
		if( .not. isUnit(U_mat) ) 					write(*,*)	"[wannInterpolator]: eigen solver gives non unitary U matrix"
		!
		!
		return
	end subroutine	

	subroutine calcVeloNew(ki, En_vec, U_int, Ha_mat, A_mat, v_mat)
		integer,		intent(in)		::	ki
		real(dp),		intent(in)		::	En_vec(:)
		complex(dp),	intent(in)		::	U_int(:,:), Ha_mat(:,:,:), A_mat(:,:,:)
		complex(dp),	intent(out)		::	v_mat(:,:,:)
		complex(dp),	allocatable		:: 	Hbar(:,:,:), Abar(:,:,:), Ucjg(:,:), tmp(:,:), vec(:)
		integer							::	m, n, i
		!
		!
		allocate(		Hbar(	3,	num_wann 	,	num_wann		)	)		
		allocate(		Abar(	3,	num_wann 	,	num_wann		)	)
		allocate(		Ucjg(		num_wann	,	num_wann		)	)
		allocate(		tmp(		num_wann	,	num_wann		)	)
		allocate(		vec(		num_wann					)	)
		Ucjg			= dconjg(	transpose(U_int)	)
		!
		do i = 1, 3
			!ROTATE TO HAM GAUGE
			if( pw90GaugeB ) then
				if(i==1 .and. ki==1 )	write(*,*)	"[wannInterp/calcVeloNew]: (H) gauge velocities used, gauge back activated"

				tmp			= matmul(	Ha_mat(i,:,:)	, U_int			)	
				Hbar(i,:,:)	= matmul(	Ucjg			, tmp		)	
				!
				tmp			= matmul(	A_mat(i,:,:)		, U_int	)	
				Abar(i,:,:)	= matmul(	Ucjg				, tmp	)
				!APPLY ROTATION
				do m = 1, num_wann
					do n = 1, num_wann
						if( n==m )	v_mat(i,n,n) = Hbar(i,n,n)
						if( n/=m )	v_mat(i,n,m) = - i_dp * dcmplx( En_vec(m) - En_vec(n) ) * Abar(i,n,m) 
						!DEBUG
						if( n/=m .and. abs(Hbar(i,n,m)) > 0.1_dp ) then
							write(*,'(a,i1,a,i3,a,i3,a,f8.4,a,f8.4,a,f8.4)')"[calcVeloNeW]: found off diag band deriv i=",i,&
									" n=",n," m=",m, "Hbar_nm=",dreal(Hbar(i,n,m)), "+i*",dimag(Hbar(i,n,m))," abs=",abs(Hbar(i,n,m))
						end if
						!
					end do
				end do
			!NO GAUGE BACK
			else
				if(i==1 .and. ki==1 )	write(*,*)	"[wannInterp/calcVeloNew]: (W) gauge velocities used"
				do m = 1, num_wann
					do n = 1, num_wann
						if(n==m)	v_mat(i,n,m)	= Ha_mat(i,n,m)
						!if(n/=m)	v_mat(i,n,m) = - i_dp * dcmplx( En_vec(m) - En_vec(n) ) * A_mat(i,n,m) 
						if(n/=m)	v_mat(i,n,m) =	- i_dp *dcmplx( En_vec(m) - En_vec(n) ) * A_mat(i,n,m) 
					end do	
				end do
			end if
		end do
		!
		!
		return
	end subroutine




	subroutine gaugeBack(ki, Hw, HaW, AW, FWtens, EnH, U, AconnH, FcurvH, veloH)
		!transform from wannier gauge back to hamiltonian gauge
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: Hw(:,:)
		complex(dp),	intent(inout)	:: HaW(:,:,:), AW(:,:,:), FWtens(:,:,:,:)
		real(dp),		intent(out)		:: EnH(:)
		complex(dp),	intent(out)		:: U(:,:), AconnH(:,:,:), FcurvH(:,:,:), veloH(:,:,:)
		complex(dp),	allocatable		:: DH(:,:,:)
		!
		allocate(	DH(2,num_wann,num_wann)	)
		!
		!COPY
		U	= HW
		!GET U MAT & ENERGIES
		call eigSolverFULL(U, EnH(:))
		U = dconjg( transpose(U))
		!ROTATE WITH u
		call calcBarMat(U, HaW, AW, FWtens)
		!CONNECTION
		call calcA(EnH(:), AW, HaW, AconnH, DH)
		!VELOCITIES
		call calcVeloNew(ki, EnH, U, HaW, AW, veloH)
		!CURVATURE
		call calcCurv(FcurvH)!returns 0.0!!!!!
		!
		!		
		return
	end subroutine


















!calcRmat Helpers:






!gaugeBack HELPERS
	subroutine calcBarMat(U, HaW, AW, FW)
		!Helper for gaugeBack
		!	 for any quantity O: \bar(O) = U^dag O U
		complex(dp),	intent(in)		:: U(:,:)
		complex(dp),	intent(inout)	:: HaW(:,:,:), Aw(:,:,:), FW(:,:,:,:)
		complex(dp),	allocatable		:: Uc(:,:)
		integer							:: a,b
		!
		allocate(	Uc( size(U,2), size(U,1) )		)
		!	
		Uc	= dconjg( transpose(U)	)
		!
		do a = 1, 2
			HaW(a,:,:) 	= matmul(	HaW(a,:,:)	, 	U				)
			HaW(a,:,:)	= matmul(	Uc			, 	HaW(a,:,:)		)
			!
			AW(a,:,:) 	= matmul(	AW(a,:,:) 	,	U				)
			AW(a,:,:) 	= matmul(	Uc			,	AW(a,:,:)		)
			!
			do b = 1, 2
				FW(a,b,:,:) 	= matmul(	FW(a,b,:,:) 	,	U				)
				FW(a,b,:,:) 	= matmul(	Uc				,	FW(a,b,:,:)		)
			end do
		end do
		!
		!
		return
	end subroutine


	subroutine calcA(EnH, AW, HaW, AconnH, DH)
		! Helper for gaugeBack
		!	A^(H) = \bar{A}^(H) + i D^(H)
		real(dp),		intent(in)		:: EnH(:)
		complex(dp),	intent(in)		:: Aw(:,:,:), HaW(:,:,:)
		complex(dp),	intent(out)		:: AconnH(:,:,:), DH(:,:,:)
		integer							:: m, n 
		!
		AconnH	= dcmplx(0.0_dp)
		DH		= dcmplx(0.0_dp)
		!SET UP D MATRIX
		do m = 1, num_wann
			do n = 1, num_wann
				if( n /= m) then
					DH(1:2,n,m)	= HaW(1:2,n,m) / ( EnH(m) - EnH(n) + machineP	)
				else
					DH(:,n,m)	= dcmplx(0.0_dp)
				end if
			end do
		end do
		!
		!CALC CONNECTION
		AconnH(1:2,:,:) = AW(1:2,:,:) + i_dp * DH(1:2,:,:)
		!
		!
		return
	end subroutine


	subroutine calcVelo(EnH, AW, HaW, veloH)
		!Helper for gaugeBack
		!	v_nm	= \bar{Ha}_nm - i_dp * (Em-En) * \bar{A}_nm
		real(dp),		intent(in)		:: EnH(:)
		complex(dp),	intent(in)		:: AW(:,:,:), HaW(:,:,:)
		complex(dp),	intent(out)		:: veloH(:,:,:)
		integer							:: n, m
		!
		veloH	= dcmplx(0.0_dp)
		!
		do m = 1, num_wann
			do n = 1, num_wann
				if( n==m ) veloH(1:2,n,m)	= HaW(1:2,n,m)
				if( n/=m ) veloH(1:2,n,m)	= - i_dp * dcmplx(	EnH(m) - EnH(n)		) * AW(1:2,n,m)
			end do
		end do
		!
		return
	end subroutine


	subroutine calcCurv(FcurvH)
		complex(dp),	intent(out)		:: FcurvH(:,:,:)
		!
		!ToDO
		FcurvH	= dcmplx(0.0_dp)
		!
		return
	end subroutine

	









end module wannInterp