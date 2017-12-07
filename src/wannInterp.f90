module wannInterp
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP, myExp, myLeviCivita, nIntegrate, &
								eigSolverFULL, rotMat, myCommutat, isHermitian, isUnit
	use sysPara

	!use 
	implicit none

	private
	public ::					DoWannInterpol

	contains

!public
	subroutine DoWannInterpol(ckW, rHopp, tHopp, R_real, EnH, U_int, AconnH, FcurvH, veloH)
		complex(dp),	intent(in)		:: ckW(:,:,:), rHopp(:,:,:,:), tHopp(:,:,:)
		real(dp),		intent(in)		:: R_real(:,:)
		complex(dp),	intent(out)		:: U_int(:,:,:)
		real(dp),		intent(out)		:: EnH(:,:)
		complex(dp),	intent(out)		:: AconnH(:,:,:,:), FcurvH(:,:,:,:), veloH(:,:,:,:)
		complex(dp),	allocatable		:: HW(:,:), HaW(:,:,:), AW(:,:,:), FWtens(:,:,:,:), FWmat(:,:,:)
		integer							:: ki, a, b, c
		!
		allocate(	HW(				nWfs, 	nWfs			)		)
		allocate(	HaW(	3	,	nWfs, 	nWfs			)		)
		allocate(	AW(		3	,	nWfs, 	nWfs			)		)
		allocate(	FWtens(	3,3	,	nWfs,	nWfs			)		)
		allocate(	FWmat( 	3	,	nWfs,	nWfs			)		)
		!
		EnH	= 0.0_dp
		AconnH	= dcmplx(0.0_dp)
		FcurvH	= dcmplx(0.0_dp)
		veloH	= dcmplx(0.0_dp)
		!
		!
		do ki = 1, nK
			!call interpolateMat(ki, tHopp, rHopp, HW, HaW, AW, FWtens)
			call wannInterpolator(ki, tHopp, rHopp, R_real, EnH, U_int(:,:,ki), HW, HaW, AW, FWtens)
			if( doGaugBack ) then
				if(ki == 1) write(*,*)	"[DoWannInterpol]: start gauging back" 	
				call gaugeBack(Hw, HaW, AW, FWtens, EnH(:,ki), U_int(:,:,ki), AconnH(:,:,:,ki), FcurvH(:,:,:,ki), veloH(:,:,:,ki))	
			else
				if(ki ==1)	write(*,*)	"[DoWannInterpol]: Gauge trafo DISABLED	"
				!CONNECTION
				AconnH(1:3,:,:,ki) 		= AW(1:3,:,:)
				!VELOCITIES
				!call calcVeloNOIntP(ki, ckW, U, HaW, EnH, AconnH, veloH)
				call calcVeloNew(ki, EnH, U_int(:,:,ki), ckW, HaW, AW, veloH)
				!
				!CURVATURE TO MATRIX
				do c = 1, 3
					do b = 1, 3
						do a = 1,3
							FcurvH(c,:,:,ki)	= myLeviCivita(a,b,c) * FWtens(a,b,:,:)
						end do
					end do
				end do
			
			end if
		end do	
		!
		write(*,*)	"[DoWannInterpol]: calculated interpolated energy, connection, curvature, velocity "
		!
		return
	end subroutine




	subroutine wannInterpolator(ki, H_tb,r_tb, R_real, En_vec, U_mat, H_mat, Ha_mat, A_mat,Om_tens)
		integer,		intent(in)		::	ki
		complex(dp),	intent(in)		::	H_tb(:,:,:), r_tb(:,:,:,:)
		real(dp),		intent(in)		:: 	R_real(:,:)
		real(dp),		intent(out)		::	En_vec(:,:)
		real(dp),		allocatable		:: 	dummy(:)
		complex(dp),	intent(out)		::	U_mat(:,:), H_mat(:,:), Ha_mat(:,:,:), A_mat(:,:,:), Om_tens(:,:,:,:)
		integer							:: R, a, b, c, n, m
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
		call eigSolverFULL(U_mat(:,:),	En_vec(:,ki))
		!U_mat	= transpose( dconjg(U_mat))
		if( .not. isUnit(U_mat) ) 					write(*,*)	"[wannInterpolator]: eigen solver gives non unitary U matrix"
		!
		!
		return
	end subroutine








	subroutine calcVeloNew(ki, En_vec, U_int, ckW, Ha_mat, A_mat, v_mat)
		integer,		intent(in)		::	ki
		real(dp),		intent(in)		::	En_vec(:,:)
		complex(dp),	intent(in)		::	U_int(:,:), ckW(:,:,:), Ha_mat(:,:,:), A_mat(:,:,:)
		complex(dp),	intent(out)		::	v_mat(:,:,:,:)
		complex(dp),	allocatable		:: 	Hbar(:,:,:), Abar(:,:,:), Ucjg(:,:), tmp(:,:), vec(:)
		integer							::	m, n, i, gi
		!
		!
		allocate(		Hbar(	3,	nWfs 	,	nWfs		)	)		
		allocate(		Abar(	3,	nWfs 	,	nWfs		)	)
		allocate(		Ucjg(		nWfs	,	nWfs		)	)
		allocate(		tmp(		nWfs	,	nWfs		)	)
		allocate(		vec(		nWfs					)	)
		Ucjg			= dconjg(	transpose(U_int)	)
		v_mat			= dcmplx(0.0_dp)
		!
		if(	doVeloNUM ) then
			!INTERPOLATE VELOCITIES
			if(ki==1)	write(*,*)"[calcVeloNew]: velocities are calculated via TB approach"
			do i = 1, 3
				!v_mat(i,:,:,ki)	= Ha_mat(i,:,:)
				!
				!NO GAUGE BACK
				!do m = 1, nWfs
				!	do n = 1, nWfs
				!		if(n==m)	v_mat(i,n,m,ki)	= Ha_mat(i,n,m)
				!		if(n/=m)	v_mat(i,n,m,ki) = - i_dp * dcmplx( En_vec(m,ki) - En_vec(n,ki) ) * A_mat(i,n,m) 
				!	end do
				!end do


				!ROTATE TO HAM GAUGE
				tmp			= matmul(	Ha_mat(i,:,:)	, U_int			)	
				Hbar(i,:,:)	= matmul(	Ucjg			, tmp		)	
				!
				tmp			= matmul(	A_mat(i,:,:)		, U_int	)	
				Abar(i,:,:)	= matmul(	Ucjg				, tmp	)
				!APPLY ROTATION
				do m = 1, nWfs
					do n = 1, nWfs
						if( n==m )	v_mat(i,n,n,ki) = Hbar(i,n,n)
						if( n/=m )	v_mat(i,n,m,ki) = - i_dp * dcmplx( En_vec(m,ki) - En_vec(n,ki) ) * Abar(i,n,m) 
						!DEBUG
						if( n/=m .and. abs(Hbar(i,n,m)) > 0.1_dp ) then
							write(*,'(a,i1,a,i3,a,i3,a,f8.4,a,f8.4,a,f8.4)')"[calcVeloNeW]: found off diag band deriv i=",i,&
									" n=",n," m=",m, "Hbar_nm=",dreal(Hbar(i,n,m)), "+i*",dimag(Hbar(i,n,m))," abs=",abs(Hbar(i,n,m))
						end if
						!
					end do
				end do
				!
			end do
		else	
			!GRAD OF BASIS FUNCTIONS
			if(ki==1)	write(*,*)"[calcVeloNew]:velocities are calculated analytically, with the plane wave coefficients"
			if( nK /= nQ) write(*,*)"[calcVeloNew]: warning analytic approach does not support different k mesh spacing"
			!do gi = 1, nGq(ki)
			!	!ROTATE BACK TO (H) GAUGE
			!	do n = 1, nWfs
			!		do m = 1, nWfs
			!			tmp(m,n)	= dconjg(ckW(gi,m,ki)) * ckW(gi,n,ki)
			!		end do
			!	end do
			!	tmp		= matmul( tmp	,	Ucjg 	)
			!	tmp		= matmul( U_int	,	tmp	)
			!	!SUM OVER BASIS FUNCTIONS
			!	do m = 1, nWfs
			!		do n = 1, nWfs
			!			v_mat(1:2,n,m,ki)	= v_mat(1:2,n,m,ki) + i_dp * Gvec(1:2,gi,ki) * tmp(n,m)
			!		end do
			!	end do
			!end do

			!OLD VERSION (no gauge back)
			do m = 1, nWfs
				do n = 1, nWfs
					do gi = 1 , nGq(ki)
						v_mat(1:2,n,m,ki) = v_mat(1:2,n,m,ki) +  dconjg(ckW(gi,n,ki)) *  ckW(gi,m,ki) * i_dp * Gvec(1:2,gi,ki)
					end do
				end do
			end do
		end if	
		!
		return
	end subroutine

























	subroutine calcVeloNOIntP(ki, ckW, U, HaW, EnH, AconnH, veloH)
		!deprecated, use calcVeloNew
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: ckW(:,:,:), U(:,:), HaW(:,:,:), AconnH(:,:,:,:)
		real(dp),		intent(in)		:: EnH(:,:) 
		complex(dp),	intent(out)		:: veloH(:,:,:,:)
		complex(dp),	allocatable		:: Hbar(:,:,:), Abar(:,:,:)
		integer							:: n,m, gi, i

		allocate(	Hbar( size(HaW,1),size(HaW,2),size(HaW,3) 			)			)
		allocate(	Abar( size(AconnH,1),size(AconnH,2),size(AconnH,3) 	)			)

		!TB approach
		if( doVeloNUM ) then
			if(ki==1)	write(*,*)"[calcVeloNOIntP]: velocities are calculated via TB approach"
			
			do i = 1, 3
				Hbar(i,:,:)	= matmul(	Haw(i,:,:),	U	)
				Hbar(i,:,:)	= matmul( dconjg(transpose(U)),	Hbar(i,:,:)	)
				Abar(i,:,:)	= matmul(	AconnH(i,:,:,ki),	U	)
				Abar(i,:,:)	= matmul( dconjg(transpose(U)),	AconnH(i,:,:,ki)	)
			end do

			do m = 1, nWfs
				do n = 1, nWfs
					if( n==m ) veloH(1:2,n,n,ki)	= Hbar(1:2,n,n)
					if( n/=m ) veloH(1:2,n,m,ki)	= - i_dp * ( EnH(m,ki)-EnH(n,ki) ) * AconnH(1:2,n,m,ki) 
				end do
			end do
		!ANALYTIC APPROACH
		else	
			if(ki==1)	write(*,*)"[calcVeloNOIntP]: velocities are calculated analytically"
			if( nK /= nQ) write(*,*)"[calcVeloNOIntP]: warning analytic approach does not support different k mesh spacing"
			do m = 1, nWfs
				do n = 1, nWfs
					do gi = 1 , nGq(ki)
						veloH(1:2,n,m,ki) = veloH(1:2,n,m,ki) &
											+  dconjg(ckW(gi,n,ki)) *  ckW(gi,m,ki) * i_dp * Gvec(1:2,gi,ki)
					end do
				end do
			end do
		end if

		return
	end subroutine








!privat


	subroutine interpolateMat(ki, tHopp, rHopp, HW, HaW, AW, FW)
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: tHopp(:,:,:), rHopp(:,:,:,:)
		complex(dp),	intent(out)		:: HW(:,:), HaW(:,:,:), AW(:,:,:), FW(:,:,:,:)
		integer							:: R, a, b
		complex(dp)						:: phase
		!
		HW	= dcmplx(0.0_dp)
		HaW	= dcmplx(0.0_dp)
		AW  = dcmplx(0.0_dp)
		FW	= dcmplx(0.0_dp)
		!
		do R = 1, nSC
			phase			= myExp(	dot_product(kpts(:,ki),Rcell(:,R))	)  ! / dsqrt(real(nSC,dp) )
			!HAM
			HW(:,:)		= HW(:,:) 		+ phase		 					* tHopp(:,:,R)
			!HAM DERIVATIVE
			HaW(1,:,:)	= HaW(1,:,:)	+ phase * i_dp *  Rcell(1,R) 	* tHopp(:,:,R) 
			HaW(2,:,:)	= HaW(2,:,:)	+ phase * i_dp *  Rcell(2,R) 	* tHopp(:,:,R) 
			!CONNECTION
			AW(1,:,:)	= AW(1,:,:) 	+ phase 						* rHopp(1,:,:,R) 
			AW(2,:,:)	= AW(2,:,:) 	+ phase 						* rHopp(2,:,:,R)
			!CURVATURE
			do a = 1, 2
				do b = 1, 2
					FW(a,b,:,:) = FW(a,b,:,:) + phase * i_dp * Rcell(a,R) * rHopp(b,:,:,R)
					FW(a,b,:,:) = FW(a,b,:,:) - phase * i_dp * Rcell(b,R) * rHopp(a,:,:,R) 	
				end do
			end do
		end do
		!
		!
		return
	end subroutine



	subroutine gaugeBack(Hw, HaW, AW, FWtens, EnH, U, AconnH, FcurvH, veloH)
		!transform from wannier gauge back to hamiltonian gauge
		complex(dp),	intent(in)		:: Hw(:,:)
		complex(dp),	intent(inout)	:: HaW(:,:,:), AW(:,:,:), FWtens(:,:,:,:)
		real(dp),		intent(out)		:: EnH(:)
		complex(dp),	intent(out)		:: U(:,:), AconnH(:,:,:), FcurvH(:,:,:), veloH(:,:,:)
		complex(dp),	allocatable		:: DH(:,:,:)
		integer							:: ki
		!
		allocate(	DH(2,nWfs,nWfs)	)





		!1. CALC 
		U	= HW
		!GET U MAT & ENERGIES
		call eigSolverFULL(U, EnH)
		U = dconjg( transpose(U))
		!ROTATE WITH u
		call calcBarMat(U, HaW, AW, FWtens)
		!CONNECTION
		call calcA(EnH, AW, HaW, AconnH, DH)
		!VELOCITIES
		call calcVelo(EnH, AW, HaW, veloH)
		!CURVATURE
		call calcCurv(FWtens, DH, AW, FcurvH)

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
		do m = 1, nWfs
			do n = 1, nWfs
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
		do m = 1, nWfs
			do n = 1, nWfs
				if( n==m ) veloH(1:2,n,m)	= HaW(1:2,n,m)
				if( n/=m ) veloH(1:2,n,m)	= - i_dp * dcmplx(	EnH(m) - EnH(n)		) * AW(1:2,n,m)
			end do
		end do
		!
		return
	end subroutine


	subroutine calcCurv(FW, DH, AW, FcurvH)
		complex(dp),	intent(in)		:: FW(:,:,:,:), DH(:,:,:), AW(:,:,:)
		complex(dp),	intent(out)		:: FcurvH(:,:,:)
		complex(dp),	allocatable		:: FH(:,:,:,:)
		!
		!ToDO
		FcurvH	= dcmplx(0.0_dp)

		return
	end subroutine


!!CONVERT OMEGA TENSOR TO VECTOR
!		do m = 1, nWfs
!			do n = 1, nWfs
!				do ki = 1, nK
!					!EVAL CROSS PRODUCT
!					do c = 1,3
!						do b= 1,3
!							do a=1,3
!								if( myLeviCivita(a,b,c) /= 0) then
!									Fh(c,ki,n,m) = Fh(c,ki,n,m) + myLeviCivita(a,b,c) * FhTens(a,b,ki,n,m)
!								end if
!							end do 
!						end do
!					end do
!					!
!				end do
!			end do
!		end do
!
	









end module wannInterp