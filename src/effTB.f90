module effTB
	!module is deprecated and should be removed in future
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP, myExp, isHermitian
	use sysPara
	use w90Interface,	only:	readFDscheme
	use planeWave,		only:	calcConnOnCoarse

	implicit none

	private
	public	::					TBviaKspace
	contains	


	subroutine TBviaKspace(ckQ, EnQ, Uq, tHopp, rHopp)
		complex(dp),		intent(in)		:: 	ckQ(:,:,:), Uq(:,:,:)
		real(dp),			intent(in)		:: 	EnQ(:,:)
		complex(dp),		intent(out)		:: 	tHopp(:,:,:), rHopp(:,:,:,:)
		complex(dp),		allocatable		:: 	Atmp(:,:,:,:)
		complex(dp),		allocatable		:: 	Htmp(:,:,:)
		real(dp),			allocatable		::	b_k(:,:), w_b(:)
		integer,			allocatable		:: 	nnlist(:,:), nncell(:,:,:)
		complex(dp)							:: 	phase
		integer								:: 	nntot, R, qi
		!
		allocate(	Atmp(	2,		nWfs,	nWfs,	nQ	)			)
		allocate(	Htmp(			nWfs,	nWfs,	nQ	)			)
		

		
		!SET UP K SPACE QUANTITIES
		call readFDscheme(nntot, nnlist, nncell, b_k, w_b)
		call calcConnOnCoarse(ckQ, nntot, nnlist, nncell, b_k, w_b, Atmp)
		call calcHtmp(EnQ, Uq, Htmp)
		!FT TO REAL SPACE
		rHopp	= dcmplx(0.0_dp)
		tHopp	= dcmplx(0.0_dp)
		do R = 1, nSC
			do qi = 1, nQ
				phase			= myExp( -1.0_dp * dot_product(qpts(1:2,qi),Rcell(1:2,R))		)  /dcmplx(size(tHopp,3))!/ dsqrt(real(nSC,dp) )
				!
				!RHOPP
				rHopp(1,:,:,R)	= rHopp(1,:,:,R) 	+ phase * Atmp(1,:,:,qi) 
				rHopp(2,:,:,R)	= rHopp(2,:,:,R) 	+ phase * Atmp(2,:,:,qi) 
				!
				!THOPP
				tHopp(:,:,R)	= tHopp(:,:,R)		+ phase * Htmp(:,:,qi)
			end do
		end do
		!
		!
		return
	end subroutine

	
!private

	subroutine calcHtmp(EnQ, Uq, Htmp)
		!Slater Koster Interpolation of tight binding paramters
		!	Htmp = Uq^dag Eq Uq (wann review)
		real(dp),		intent(in)		:: EnQ(:,:)
		complex(dp),	intent(in)		:: Uq(:,:,:)
		complex(dp),	intent(out)		:: Htmp(:,:,:)
		complex(dp),	allocatable		:: Udag(:,:), U(:,:), Ediag(:,:), tmp(:,:)
		integer							:: qi, n
		!
		allocate(	U(		nWfs, nwfs	)		)
		allocate(	tmp(	nWfs, nWfs	)		)
		allocate(	Udag(	nWfs, nwfs	)		)
		allocate(	Ediag(	nWfs, nWfs	)		)
		Htmp	= dcmplx(0.0_dp)
		!
		!
		!SET UP
		do qi = 1, nQ
			U(:,:)	= Uq(:,:,qi)
			Udag	= dconjg( transpose( U ) )
			Ediag	= dcmplx(0.0_dp)
			do n = 1, nWfs
				Ediag(n,n)	= EnQ(n,qi)
			end do
			tmp				= matmul(	Ediag	, 	U		)
			Htmp(:,:,qi)	= matmul( 	Udag	, 	tmp		)
			!
			if( .not. isHermitian(Htmp(:,:,qi)) ) 	then
				write(*,'(a,i3)')	"[TBviaKspace]: generated H matrix is not hermitian, at qi=",qi
			end if
		end do
		!
		!
		return
	end subroutine





end module
