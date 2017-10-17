module wannInterp
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP, myExp, myLeviCivita, nIntegrate, eigSolver, rotMat, myCommutat
	use sysPara
	use output,			only: 	writeEnH

	!use 
	implicit none

	private
	public ::					DoWannInterpol

	contains

!public
	subroutine DoWannInterpol(rHopp, tHopp, EnH, AconnH, FcurvH, veloH)
		complex(dp),	intent(in)		:: rHopp(:,:,:,:), tHopp(:,:,:)
		real(dp),		intent(out)		:: EnH(:,:)
		complex(dp),	intent(out)		:: AconnH(:,:,:,:), FcurvH(:,:,:,:)
		complex(dp),	intent(out)		:: veloH(:,:,:,:)
		complex(dp),	allocatable		:: U(:,:), HW(:,:), HaW(:,:,:), AW(:,:,:), FW(:,:,:,:)
		integer							:: ki
		!
		
		allocate(	U(				nWfs, 	nWfs			)		)
		allocate(	HW(				nWfs, 	nWfs			)		)
		allocate(	HaW(	2	,	nWfs, 	nWfs			)		)
		allocate(	AW(		2	,	nWfs, 	nWfs			)		)
		allocate(	Fw(		2,2	,	nWfs,	nWfs			)		)
	
		
		do ki = 1, nK
			call interpolateMat(ki, tHopp, rHopp, HW, HaW, AW, FW)
			if( doGaugBack ) then
				if(ki == 1) write(*,*)	"[DoGaugeTrafo]: start gauging back" 	
				call gaugeBack(Hw, HaW, AW, FW, EnH(:,ki), U, AconnH(:,:,:,ki), FcurvH(:,:,:,ki), veloH(:,:,:,ki))	
			else
				if(ki ==1)	write(*,*)	"[DoGaugeTrafo]: Gauge trafo DISABLED	"
				!AconnH(1:2,:,:,ki)	= AW(1:2,:,:)
				call eigSolver(HW,EnH(:,ki))
				AconnH(1:2,:,:,ki) 		= AW(1:2,:,:)
				veloH(1:2,:,:,ki) 		= HaW(1:2,:,:) 
				FcurvH(1:2,:,:,ki)		= dcmplx(0.0_dp)
				!call calcCurv(FW, DH, AW, FcurvH)
			end if
		end do	

		write(*,*)	"[DoGaugeTrafo]: calculated (H) gauge energy, connection, curvature, velocity"
		!
		if( writeBin )	call writeEnH(EnH)
		!
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
			phase			= myExp(	dot_product(kpts(:,ki),Rcell(:,R))	)   / dsqrt(real(nQ,dp) )
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



	subroutine gaugeBack(Hw, HaW, AW, FW, EnH, U, AconnH, FcurvH, veloH)
		!transform from wannier gauge back to hamiltonian gauge
		complex(dp),	intent(in)		:: Hw(:,:)
		complex(dp),	intent(inout)	:: HaW(:,:,:), AW(:,:,:), FW(:,:,:,:)
		real(dp),		intent(out)		:: EnH(:)
		complex(dp),	intent(out)		:: U(:,:), AconnH(:,:,:), FcurvH(:,:,:), veloH(:,:,:)
		complex(dp),	allocatable		:: DH(:,:,:)
		integer							:: ki
		!
		allocate(	DH(2,nWfs,nWfs)	)


		!1. CALC 
		U	= HW
		!GET U MAT & ENERGIES
		call eigSolver(U, EnH)
		U = dconjg( transpose(U))
		!ROTATE WITH u
		call calcBarMat(U, HaW, AW, FW)
		!CONNECTION
		call calcA(EnH, AW, HaW, AconnH, DH)
		!VELOCITIES
		call calcVelo(EnH, AW, HaW, veloH)
		!CURVATURE
		call calcCurv(FW, DH, AW, FcurvH)

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
				veloH(1:2,n,m)	= HaW(1:2,n,m) - i_dp * dcmplx(	EnH(m) - EnH(n)		) * AW(1:2,n,m)
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