module peierls
	use mathematics,	only:	dp, PI_dp, i_dp, myExp, crossP, nIntegrate, eigSolver
	use sysPara
	use	wannier,		only:	calcHopping
	use gaugeTrafo,		only:	DoWannInterpol
	use	polarization,	only:	calcPolViaA
	use output,			only:	writePeierls
	implicit none
	

	private
	public ::	peierlsMethod


	contains







!public:
	subroutine	peierlsMethod(wnf, pPei)
		complex(dp),	intent(in)		:: wnf(:,:,:)	! tHopp(nWfs,nWfs,nSC)
		real(dp),		intent(out)		:: pPei(3)
		complex(dp),	allocatable		:: tHopp(:,:,:), rHopp(:,:,:,:), AconnP(:,:,:,:), FcurvP(:,:,:,:), veloP(:,:,:,:)
		real(dp),		allocatable		:: EnP(:,:)
		integer							:: R, ki
		complex(dp)						:: phase
		!
		allocate(			tHopp(		nWfs	, 	nWfs	,	nSc				)			)
		allocate(			rHopp(	2	,	nWfs, 	nWfs, 	nSC					)			)		
		allocate(			EnP(nWfs,nK)			)
		allocate(			AconnP(3,nWfs,nWfs,nK)	)
		allocate(			FcurvP(3,nWfs,nWfs,nK)	)
		allocate(			veloP(3,nWfs,nWfs,nK)	)
		!
		AconnP 	= dcmplx(0.0_dp) 
		FcurvP 	= dcmplx(0.0_dp)
		veloP	= dcmplx(0.0_dp)



		write(*,*)	"[peierlsMethod]: start with peierls sub"

		!GET TIGHT BINDING MAT ELEMENTS
		call calcHopping(wnf, tHopp, rHopp)
		!
		!DO PEIERLS SUBSTITUTION
		do R = 1, nSC
			tHopp(:,:,R)	= tHopp(:,:,R) * shift(R0,R)
		end do
		
		!INTERPOLATE
		call DoWannInterpol(rHopp, tHopp, EnP, AconnP, FcurvP, veloP)


		!CALC POL
		call calcPolViaA(AconnP, pPei)

		write(*,*)	"[peierlsMethod]: writing done, by.."
		deallocate(	EnP		)
		deallocate(	AconnP	)
		deallocate(	FcurvP	)
		deallocate(	veloP	)
		!
		!
		return
	end subroutine










!privat:
	real(dp) function shift(ri, rj)
		!integrates the vector potential A(r) analytically
		!	A(r)	= -0.5 cross_p[r,B]
		!	return 	= Integrate^\vec{rMax]_\vec{rMin} 		\vec{A(r')}.\vec{dr'}
		!			= -0.5 cross_p[rMin,B].(rMax - rMin)
		integer,		intent(in)		:: ri, rj
		real(dp)						:: rMin(2), rMax(2), rU(3), rL(3), integrateA
		rMin		= Rcell(:,ri)
		rMax		= Rcell(:,rj)
		rU(1:2)		= rMax(1:2)
		rU(3)		= 0.0_dp
		rL(1:2)		= rMin(1:2)
		rL(3)		= 0.0_dp
		!
		integrateA	= -0.5_dp * dot_product( crossP(rL,Bext)	, (rU-rL)	)
		shift		= myExp( integrateA)
		!
		return
	end function


	

	subroutine genUnk(ki, Ham, unkP)
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: Ham(:,:)
		complex(dp),	intent(out)		:: unkP(:,:)
		complex(dp),	allocatable		:: basVec(:)
		real(dp)						:: kVec(2)
		complex(dp)						:: phase
		integer							:: i, xi
		!
		allocate(basVec(nWfs))
		unkP	= dcmplx(0.0_dp)
		!
		do xi = 1, nR
			!CALC BASIS VECTOR
			basVec	= dcmplx(0.0_dp)
			do i =1, nWfs
				kVec(:) = kpts(:,ki) + Gvec(:,i)
				if( norm2(kVec) < Gcut ) then
					basVec(i) 		=  myExp( dot_product( kVec, rpts(:,xi) )		)
				else
					basVec(i) 		= dcmplx( 0.0_dp )
				end if
			end do
			!
			!GENERATE UNKs
			phase			= myExp( -1.0_dp * dot_product( kpts(:,ki), rpts(:,xi) )		)
			unkP(xi,:)		= phase * dsqrt(real(nSC,dp)) * matmul(basVec,Ham) 
		end do
		!
		!
		return
	end subroutine



end module peierls