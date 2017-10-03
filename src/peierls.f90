module peierls
	use mathematics,	only:	dp, PI_dp, i_dp, myExp, crossP, nIntegrate, eigSolver
	use sysPara
	use gaugeTrafo,		only:	DoGaugeTrafo
	use	polarization,	only:	calcPolViaA
	use output,			only:	writePeierls
	implicit none
	

	private
	public ::	peierlsMethod


	contains







!public:
	subroutine	peierlsMethod(tHopp, pPei)
		complex(dp),	intent(inout)	:: tHopp(:,:,:)	! tHopp(nWfs,nWfs,nSC)
		real(dp),		intent(out)		:: pPei(3)
		complex(dp),	allocatable		:: Ham(:,:), unkP(:,:,:), AconnP(:,:,:,:), FcurvP(:,:,:,:), veloP(:,:,:,:)
		real(dp),		allocatable		:: EnP(:,:)
		integer							:: R, ki
		complex(dp)						:: phase
		!
		allocate( 	Ham(nWfs,nWfs)			)
		allocate(	EnP(nWfs,nK)			)
		allocate(	unkP(nR, nWfs, nK)		)
		allocate(	AconnP(2,nWfs,nWfs,nK)	)
		allocate(	FcurvP(2,nWfs,nWfs,nK)	)
		allocate(	veloP(2,nWfs,nWfs,nK)	)
		!
		write(*,*)	"[peierlsMethod]: start with peierls sub"

		!
		!DO PEIERLS SUBSTITUTION
		do R = 1, nSC
			tHopp(:,:,R)	= tHopp(:,:,R) * shift(R0,R)
		end do
		write(*,*)	"[peierlsMethod]: done with sub, solve Ham now"
		!SET UP K SPACE HAMILTONIAN & SOLVE
		do ki = 1, nK
			Ham				= dcmplx(0.0_dp)
			do R = 1, nSC
				phase		= myExp( 	dot_product ( kpts(:,ki), Rcell(:,R) )		)
				Ham(:,:)	= Ham(:,:)	+ phase * tHopp(:,:,R)
			end do
			!SOLVE:
			call eigSolver(	Ham, EnP(:,ki)	) 
			call genUnk(ki, Ham, unkP(:,:,ki))
		end do
		!
		!write(*,*)	"[peierlsMethod]: solved Ham, calc connection, etc."
		!!CALC CONNECTION & POL
		!call DoGaugeTrafo(unkP, tHopp, EnP, AconnP, FcurvP, veloP)
		!call calcPolViaA(AconnP, pPei(1:2))
		!write(*,*)	"[peierlsMethod]: calculated pol"
		!!
		!call writePeierls(unkP, AconnP, FcurvP)
		write(*,*)	"[peierlsMethod]: writing done, by.."
		deallocate( Ham		)
		deallocate(	EnP		)
		deallocate(	unkP	)
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