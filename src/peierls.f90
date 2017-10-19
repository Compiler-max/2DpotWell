module peierls
	use mathematics,	only:	dp, PI_dp, i_dp, myExp, crossP, nIntegrate, eigSolver
	use sysPara
	use	wannier,		only:	calcHopping
	use effTB,			only:	calcConnOnCoarse
	use wannInterp,		only:	DoWannInterpol
	use	polarization,	only:	calcPolViaA
	use output,			only:	writePeierls
	implicit none
	

	private
	public ::	peierlsMethod


	contains







!public:
	subroutine	peierlsMethod(tHopp, pPei)
		complex(dp),	intent(in)		:: tHopp(:,:,:)	! tHopp(nWfs,nWfs,nSC)
		real(dp),		intent(out)		:: pPei(3)
		complex(dp),	allocatable		:: Hp(:,:), unkP(:,:,:), tshift(:,:,:), AconnP(:,:,:,:)
		real(dp),		allocatable		:: EnP(:,:)
		integer							:: R, ki
		complex(dp)						:: phase
		real(dp)						:: shft
		!
		allocate(			Hp(			nWfs	,	nWfs				)			)
		allocate(			unkP(		nR		,	nWfs	,	nK		)			)
		allocate(			tshift(		nWfs	, 	nWfs	,	nSc		)			)
		allocate(			EnP(					nWfs	,	nK		)			)
		allocate(			AconnP(3,	nWfs	,	nWfs	,	nK		)			)
		!
		pPei	= 0.0_dp
		write(*,*)	"[peierlsMethod]: start with peierls sub"

		
		!
		!DO PEIERLS SUBSTITUTION
		do R = 1, nSC
			shft			= shift(R0,R)
			tshift(:,:,R)	= tHopp(:,:,R) * shft
			write(*,'(a,i3,a,f10.4)')	"[peierlsMethod]: R=",R," shift=",shft
		end do
		write(*,*)	"[peierlsMethod]: substiution of hopping parameters done"

		!GET CONNECTION
		AconnP	= dcmplx(0.0_dp)
		do ki = 1, nK
			Hp	= dcmplx(0.0_dp)
			!FT to k space
			do R = 1, nSC
				phase	= myExp(	dot_product(kpts(:,ki),Rcell(:,R))	)  / dsqrt(real(nQ,dp) )
				!Hp(:,:)	= Hp(:,:) + phase * )
			end do
			!SOLVE ELECTRONIC STRUCTURE
			call eigSolver(Hp,EnP(:,ki))
			call genUnk(ki, Hp, unkP(:,:,ki))
		end do
		!GENERATE CONNECTION
		call calcConnOnCoarse(unkP, AconnP)

		if( nK /= nQ ) then
			write(*,*)	"[peierlsMethod]: WARNING, coarse & mesh do not have same grid spacing... "
			write(*,*)	"[peierlsMethod]: ... the FD implementation of Berry conn. is wrong in that case!!! "
			write(*,*)	"[peierlsMethod]: ... will set pPei to zero "
		else
			write(*,*)	"[peierlsMethod]: calculated Berry connection."
		end if

		!CALC POL
		call calcPolViaA(AconnP, pPei)

		if( nK /= nQ ) pPei = 0.0_dp

		write(*,*)	"[peierlsMethod]: calculated polarization, by.."
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
		shift		= myExp( integrateA )
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