module peierls
	use mathematics,	only:	dp, PI_dp, i_dp, myExp, crossP, nIntegrate, eigSolver, isUnit
	use sysPara
	use projection,		only:	projectUnk
	use effTB,			only:	calcConnOnCoarse
	use wannInterp,		only:	DoWannInterpol
	use	polarization,	only:	calcPolViaA
	use output,			only:	writePeierls
	implicit none
	

	private
	public ::	peierlsMethod


	contains







!public:
	subroutine	peierlsMethod(ck, tHopp, pPei)
		complex(dp),	intent(in)		:: ck(:,:,:), tHopp(:,:,:)	! tHopp(nWfs,nWfs,nSC)
		real(dp),		intent(out)		:: pPei(3)
		complex(dp),	allocatable		:: Hp(:,:), ckP(:,:,:),Up(:,:,:), tshift(:,:,:), AconnP(:,:,:,:)
		real(dp),		allocatable		:: EnP(:,:)
		integer							:: R, ki, ri
		complex(dp)						:: phase
		real(dp)						:: shft
		!
		allocate(			Hp(			nWfs	,	nWfs				)			)
		allocate(			ckP(		nWfs	,	nWfs	,	nK		)			)
		allocate(			tshift(		nWfs	, 	nWfs	,	nSc		)			)
		allocate(			EnP(					nWfs	,	nK		)			)
		allocate(			Up(			nWfs	,	nWfs	, 	nK		)			)
		allocate(			AconnP(3,	nWfs	,	nWfs	,	nK		)			)
		!
		pPei	= 0.0_dp
		write(*,*)	"[peierlsMethod]: start with peierls sub"

		
		!
		!DO PEIERLS SUBSTITUTION
		do R = 1, nSC
			shft			= 1.0_dp!shift(R0,R)
			tshift(:,:,R)	= tHopp(:,:,R) * shft
			write(*,'(a,i3,a,f10.4)')	"[peierlsMethod]: R=",R," shift=",shft
		end do
		write(*,*)	"[peierlsMethod]: substiution of hopping parameters done"

		
		!ELECTRONIC STRUCTURE
		do ki = 1, nK
			!SET UP HAMILTONIAN
			Hp	= dcmplx(0.0_dp)
			do R = 1, nSC
				phase	= myExp( dot_product(kpts(:,ki),Rcell(:,R))	) / dsqrt(real(nSC,dp))
				Hp(:,:)	= Hp(:,:) + phase * tshift(:,:,R)
			end do
			!SOLVE HAM	
			call eigSolver(Hp(:,:),EnP(:,ki))
			!
			ckP(:,:,ki)	= HP(:,:)
			if( .not. isUnit(ckP(:,:,ki))	) write(*,*) "[peierlsMethod]: ckP not unitary at ki=",ki
		end do


	


		!GENERATE CONNECTION
		AconnP	= dcmplx(0.0_dp)
		call calcConnOnCoarse(ckP, AconnP) 


		!CALC POL
		call calcPolViaA(AconnP, pPei)

		!WRITE UNKs & ENERGIES
		if( writeBin ) call writePeierls(ckP, EnP)


		!DEBUG
		if( nK /= nQ ) then
			write(*,*)	"[peierlsMethod]: WARNING, coarse & mesh do not have same grid spacing... "
			write(*,*)	"[peierlsMethod]: ... the FD implementation of Berry conn. is wrong in that case!!! "
			write(*,*)	"[peierlsMethod]: ... will set pPei to zero "
			pPei = 0.0_dp
		else
			write(*,*)	"[peierlsMethod]: calculated Berry connection."
		end if
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


	

	!subroutine genUnk(ki, Ham, unkP)
	!	integer,		intent(in)		:: ki
	!	complex(dp),	intent(in)		:: Ham(:,:)
	!	complex(dp),	intent(out)		:: unkP(:,:)
	!	complex(dp),	allocatable		:: basVec(:)
	!	real(dp)						:: kVec(2)
	!	complex(dp)						:: phase
	!	integer							:: i, xi
	!	!
	!	allocate(basVec(nWfs))
	!	unkP	= dcmplx(0.0_dp)
	!	!
	!	do xi = 1, nR
	!		!CALC BASIS VECTOR
	!		basVec	= dcmplx(0.0_dp)
	!		do i =1, nWfs
	!			kVec(:) = kpts(:,ki) + Gvec(:,i)
	!			if( norm2(kVec) < Gcut ) then
	!				basVec(i) 		=  myExp( dot_product( kVec, rpts(:,xi) )		)
	!			else
	!				basVec(i) 		= dcmplx( 0.0_dp )
	!			end if
	!		end do
	!		!
	!		!GENERATE UNKs
	!		phase			= myExp( -1.0_dp * dot_product( kpts(:,ki), rpts(:,xi) )		)
	!		unkP(xi,:)		= phase * dsqrt(real(nSC,dp)) * matmul(basVec,Ham) 
	!	end do
	!	!
	!	!
	!	return
	!end subroutine



end module peierls