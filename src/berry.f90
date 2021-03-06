module berry
	!module contains methods related to the Berry phase theory of modern polarization
	!	i.e. calculation of connection, velocities, curvatures and polarization
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, myLeviCivita, nIntegrate
	use sysPara
	use effTB,			only:	TBviaKspace
	use wannInterp,		only:	DoWannInterpol
	use polarization,	only:	calcPolViaA
	use semiClassics,	only:	calcFirstOrdP
	use peierls,		only:	peierlsMethod
	use output,			only:	writeConnCurv
	implicit none

	private
	public	::	berryMethod
	contains




!public
	subroutine berryMethod(ckW, EnQ, Uq, pBerry, pNiuF2, pNiuF3, pPei)
		complex(dp),	intent(in)		:: ckW(:,:,:), Uq(:,:,:)
		real(dp),		intent(in)		:: EnQ(:,:)
		real(dp),		intent(out)		:: pBerry(2), pNiuF2(3), pNiuF3(3), pPei(3)
		real(dp),		allocatable		:: EnK(:,:)
		complex(dp),	allocatable		:: AconnK(:,:,:,:), FcurvK(:,:,:,:), veloK(:,:,:,:) , tHopp(:,:,:), rHopp(:,:,:,:) 
		!
		allocate(			tHopp(					nWfs	, 	nWfs	,	nSc		)			)
		allocate(			rHopp(		2		,	nWfs	, 	nWfs	, 	nSC		)			)			
		allocate(			EnK(					nWfs	, 				nK		)			)
		allocate(			AconnK(		3		, 	nWfs	,	nWfs	,	nK		)			)
		allocate(			FcurvK(		3		, 	nWfs	, 	nWfs	,	nK		)			)
		allocate(			veloK(		3		, 	nWfs	,	nWfs	,	nK		)			)
					
		!
		write(*,*)	"[berrryMethod]: hello from Berry"
		

		!SET UP EFFECTIVE TIGHT BINDING MODELL
		call TBviaKspace(ckW, EnQ, Uq, tHopp, rHopp)
		write(*,*)	"[berryMethod]: set up effective tight binding model (k-Space method)"

		!INTERPOLATE CONN,CURV, VELO
		call DoWannInterpol(ckW, rHopp, tHopp, EnK, AconnK, FcurvK, veloK)
		write(*,*)	"[berrryMethod]: interpolation done"

		
		
		!INTEGRATE CONNECTION
		call calcPolViaA(AconnK,pBerry)
		write(*,'(a,f12.8,a,f12.8,a)')	"[berrryMethod]: calculated zero order pol=(",pBerry(1),", ",pBerry(2),")."
		



		!1st ORDER SEMICLASSICS
		if(doNiu) then
			write(*,*)	"[berrryMethod]: now calc first order pol"
			call calcFirstOrdP(FcurvK, AconnK, veloK, EnK, pNiuF2, pNiuF3)
			write(*,'(a,f8.4,a,f8.4,a,f8.4,a)')	"[berryMethod]: pNiuF2=(",pNiuF2(1),", ",pNiuF2(2),", ",pNiuF2(3),")."
			write(*,'(a,f8.4,a,f8.4,a,f8.4,a)')	"[berryMethod]: pNiuF3=(",pNiuF3(1),", ",pNiuF3(2),", ",pNiuF3(3),")."
		end if


		!PEIERLS SUBSTITUTION
		if(doPei) then
			write(*,*)	"[berrryMethod]: now calc first order pol via peierls sub."
			call peierlsMethod(ckW, tHopp, pPei)
			write(*,'(a,f8.4,a,f8.4,a,f8.4,a)')	"[berryMethod]: pPei=(",pPei(1),", ",pPei(2),", ",pPei(3),")."
		end if

		

		!OUTPUT
		!call writeConnCurv(AconnK, FcurvK)
		write(*,*)	"[berrryMethod]: all done"
		!
		!
		!
		return
	end subroutine



end module berry