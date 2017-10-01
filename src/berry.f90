module berry
	!module contains methods related to the Berry phase theory of modern polarization
	!	i.e. calculation of connection, velocities, curvatures and polarization
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, myLeviCivita, nIntegrate
	use sysPara
	use gaugeTrafo,		only:	DoGaugeTrafo
	use polarization,	only:	calcPolViaA
	use semiClassics,	only:	calcFirstOrdP
	use output,			only:	writeConnCurv
	implicit none

	private
	public	::	berryMethod
	contains




!public
	subroutine berryMethod(unkW, tHopp, pBerry, pNiu)
		complex(dp),	intent(in)		:: unkW(:,:,:), tHopp(:,:,:)
		real(dp),		intent(out)		:: pBerry(2), pNiu(3)
		real(dp),		allocatable		:: EnH(:,:)
		complex(dp),	allocatable		:: AconnH(:,:,:,:), FcurvH(:,:,:,:), veloH(:,:,:,:)
		!
		allocate(			EnH(		nWfs	, nK							)			)
		allocate(			AconnH(		3		, nWfs,nWfs			,	nK		)			)
		allocate(			FcurvH(		3		, nWfs, nWfs		,	nK		)			)
		allocate(			veloH(		3		, nWfs,nWfs			,	nK		)			)
		!
		write(*,*)	"[berrryMethod]: hello from Berry"
		!GET CONNECTION, CURVATURE & VELOCITY
		call DoGaugeTrafo(unkW, tHopp, EnH, AconnH, FcurvH, veloH)

		write(*,*)	"[berrryMethod]: interpolation done"
		!INTEGRATE CONNECTION
		call calcPolViaA(AconnH,pBerry)
		write(*,*)	"[berrryMethod]: calculated zero order pol"
		
		!ToDo: SWITCH FOR USING NIU
		if(doNiu) then
			write(*,*)	"[berrryMethod]: now calc first order pol"
			call calcFirstOrdP(FcurvH, AconnH, veloH, EnH, pNiu)
			write(*,*)	"[berrryMethod]: all done"
		end if


		!OUTPUT
		call writeConnCurv(AconnH, FcurvH)
		!
		!
		return
	end subroutine







end module berry






	!subroutine berryMethod(unk, En, pBerry, pNiu)
	!	complex(dp),	intent(in)		:: unk(:,:,:)
	!	real(dp),		intent(in)		:: En(:,:)
	!	real(dp),		intent(out)		:: pBerry(2), pNiu(3)
	!	complex(dp),	allocatable		:: velo(:,:,:,:)
	!	real(dp),		allocatable		:: Aw(:,:,:,:), FW(:,:,:)
	!	!
	!	allocate(			Aw(			3		, nWfs,nWfs			,	nQ		)			)
	!	allocate(			velo(		3		, nWfs,nWfs			,	nQ		)			)
	!	allocate(			FW(			3		, nWfs				,	nQ		)			)
!
!	!	call calcConnOnCoarse(unk, Aw)
!	!	call calcPolViaA(Aw,pBerry)
!	!	
!	!	!ToDo: SWITCH FOR USING NIU
!	!	if(doNiu) then
!	!		!Get zero order velocities and curvatures
!	!		!call calcVeloMat(unk, VeloBwf, velo)
!	!		call FDvelocities(unk, velo)
!	!		call calcCurv(En, velo, Fw)
!	!		!use them for calc of first order pol 
!	!		call calcFirstOrdP(Fw, Aw, velo, En, pNiu)
!	!	end if
!
!
!	!	!OUTPUT
!	!	call writeConnCurv(Aw, Fw)
!	!	!
!	!	!
!	!	return
	!end subroutine









	!subroutine calcVeloMat(unk, veloBwf, Velo)
	!	!calculates matrix elements of the velocity operator
	!	!	velocity operator is analytically applied to the plane wave basis 
	!	!	and then weighted by the basCoeff obtained from the solver and stored in veloBwf
	!	!	
	!	complex(dp),		intent(in)		:: unk(:,:,:), veloBwf(:,:,:,:)	!	unk(nR, nK, nWfs) , veloBwf(2,nR,nG,nQ)
	!	complex(dp),		intent(out)		:: Velo(:,:,:,:)   !Velo(3,	nK,	nWfs	, nwFs)		
	!	integer								:: qi, n,m, ri
	!	complex(dp),		allocatable		:: fx(:), fy(:)
	!	complex(dp)							:: BWFphase
	!	!
	!	!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(m, n, qi, ri, BWFphase, fx, fy)
	!	allocate(	fx(nR)	)
	!	allocate(	fy(nR)	)
	!	!
	!	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!	do m = 1, nWfs
	!		do n = 1, nWfs
	!			do qi = 1, nQ
	!				!FILL INTEGRATION ARRAY
	!				do ri = 1, nWfs
	!					BWFphase= myExP( 	dot_product( qpts(:,qi), rpts(:,ri) ))
	!					fx(ri)	= dconjg(	BWFphase * unk(ri,n,qi)	)	* veloBwf(1,ri, m, qi )
	!					fy(ri)	= dconjg(	BWFphase * unk(ri,n,qi)	)	* veloBwf(2,ri, m, qi )
	!				end do
	!				!INTEGRATE
	!				Velo(1,n,m,qi)	= nIntegrate(nR, nRx, nRy, dx, dy, fx)
	!				Velo(2,n,m,qi)	= nIntegrate(nR, nRx, nRy, dx, dy, fy)
	!				Velo(3,n,m,qi) 	= dcmplx(0.0_dp)
	!			end do
	!		end do
	!	end do
	!	!$OMP END DO
	!	!$OMP END PARALLEL	
	!	!
	!	return
	!end subroutine









