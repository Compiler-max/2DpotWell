module blochWf
	!generates bloch and lattice periodic functions given a basCoeff matrix
	use mathematics,	only:	dp, PI_dp,i_dp, myExp, myLeviCivita, eigSolver, nIntegrate
	use sysPara,		only: 	readInp, getKindex, getRindex, &
									dim, aX, aY,vol, nAt, atR, atPos, atPot,&
									nG, nG0, Gcut, nK, nKx, nKy, nWfs, nSC, nSCx, nSCy, nR, nRx, nRy, dx, dy, dkx, dky, &
									Gvec, atPos, atR, kpts, rpts, gaugeSwitch

	implicit none

	private
	public	::	genBlochWf, BwFisLattSym, calcVeloBwf, genUnk


	contains







!public
	subroutine genBlochWf(ki,basCoeff, bWf)
		!generates the bloch wavefunctions, with  the basCoeff from eigSolver
		integer		, intent(in)	:: ki
		complex(dp)	, intent(in)	:: basCoeff(:,:)
		complex(dp)	, intent(out)	:: bWf(:,:)	!bWf(nRpts,nG)			
		complex(dp)	, allocatable	:: basVec(:)
		integer 				 	:: xi,n
		allocate(	basVec(nG)	)
		!
		do xi = 1, nR
				call calcBasVec(ki,xi, basVec)
				bWf(xi,:) = matmul(	 basVec , basCoeff	)  /  dsqrt(vol)
		end do
		!
		return 
	end


	subroutine calcVeloBwf(ki,basCoeff, velobWf)
		!generates the bloch wavefunctions, with  the basCoeff from eigSolver
		integer		, intent(in)	:: ki
		complex(dp)	, intent(in)	:: basCoeff(:,:)
		complex(dp)	, intent(out)	:: velobWf(:,:,:)	!VeloBwf(	nR		,	nK		, 2*nWfs)			
		complex(dp)	, allocatable	:: basVec(:)
		integer 				 	:: xi,n
		allocate(	basVec(nG)	)
		!
		do xi = 1, nR
				call calcVeloBasVec(ki,xi, basVec)
				velobWf(xi,ki,:) = matmul(	 basVec , basCoeff	)  /  dsqrt(vol)
		end do
		!
		return 
	end


	logical function BwFisLattSym(bWf)
		!ToDo
		!checks if bwf(k) = bwf(k+G)
		complex(dp),	intent(in)		:: bWf(:,:,:) !nR, nK , nG or nWfs
		integer							:: k00, k10, k01, k11, n ! edge point indices

		BwFisLattSym = .true.
		k00 = getKindex(	1	, 1		)
		k10	= getKindex(	nKx	, 1		)
		k01 = getKindex(	1	, nKy	)
		k11 = getKindex(	nKx , nKy	)
		write(*,'(a,i3,a,i3,a,i3,a,i3)')"[isLattSym]: k00 =",k00,", k10=",k10,", k01=",k01,", k11=",k11 


		do n = 1, size(bwf,3) ! loop states

		end do

		return
	end


	subroutine genUnk(ki, bWf, unk)
		! generates the lattice periodic part from given bloch wave functions
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: bWf(:,:) !lobWf(	nR, nWfs)
		complex(dp),	intent(out)		:: unk(:,:)   !unk(	nR, nWfs)
		integer							:: xi, n
		complex(dp)						:: phase
		!
		do n = 1, nWfs
			do xi = 1, nR
				phase = myExp( -1.0_dp 	*	 dot_product( kpts(:,ki) , rpts(:,xi)	) 			)
				unk(xi,n) = phase * bWf(xi,n)
			end do
		end do
		!
		return
	end














!privat
	subroutine calcBasVec(ki, ri, basVec)
		!calculates the basis vectors e^i(k+G).r
		!	if |k+G| is larger then the cutoff the basis vector is set to zero
		!	the cutoff enforces a symmetric base at each k point
		integer,	 intent(in)  :: ki, ri
		complex(dp), intent(out) :: basVec(:)
		real(dp)				 :: tmp(dim)
		integer 				 ::	i 
		!

		do i =1, nG
			tmp(:) = kpts(:,ki) + Gvec(:,i)
			!
			if( norm2(tmp) < Gcut ) then
				basVec(i) = myExp( 		dot_product( tmp, rpts(:,ri) )			)
			else
				basVec(i) = dcmplx( 0.0_dp )
			end if
		end do
		!
		return
	end



	!VELOCITY HELPERS
	subroutine calcVeloBasVec(ki,ri,basVec)
		!the velocity basi
		integer,		intent(in)		:: ki, ri
		complex(dp),	intent(out)		:: basVec(:)
		real(dp)				 :: tmp(2)
		integer 				 ::	i 
		!
		do i =1, nG
			tmp(:) = kpts(:,ki) + Gvec(:,i)
			!
			if( norm2(tmp) < Gcut ) then
				!X COMPONENT
				basVec(i) 		= i_dp * (	kpts(1,ki) + Gvec(1,i)	) * myExp( 		dot_product( tmp, rpts(:,ri) )			)
				!Y COMPONENT
				basVec(i+nG)	= i_dp * (	kpts(2,ki) + Gvec(2,i)	) * myExp(		dot_product( tmp, rpts(:,ri) )			)
			else
				basVec(i) = dcmplx( 0.0_dp )
			end if
		end do
		!		
		!
		return
	end



end module blochWf













