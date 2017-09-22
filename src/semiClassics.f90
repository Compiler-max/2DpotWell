module semiClassics
	!this module uses a semiclassic approach to calculate the first ordrer correction
	!	to the polariztion induced by a perturbive magnetic field
	! 	see Niu PRL 112, 166601 (2014)
	use mathematics,	only:	dp, PI_dp, i_dp, myExp, myLeviCivita
	use sysPara,		only:	Bext

	implicit none



	private
	public ::			calcFirstOrdP



	contains



!TODO CHECK INDEXING OF VELO ACONN FCURV AFTER THEY HAVE BEEN CALCULATED


!public


	subroutine	calcFirstOrdP(Fcurv, Aconn, Velo, En, p1)
		!calculates the first order polarization p1 according to
		!	P'= -int_dk [0.5 (Curv.Velo)*B_ext + a']
		real(dp),		intent(in)		::	Fcurv(:,:,:), Aconn(:,:,:,:)	!Fcurv(3,nWfs, nQ)
		complex(dp),	intent(in)		:: 	Velo(:,:,:,:)		!	 Velo(3, nWfs,nWfs, nQ)	
		real(dp),		intent(in)		::	En(:,:)				!	En(			nWfs, nQ)						
		real(dp),		intent(out)		:: 	p1(3)
		complex(dp), 	allocatable		::	f(:,:)
		real(dp)						::	pn(3)
		real(dp)						:: 	Fmat(3,3)
		complex(dp)						:: 	densCorr(3)
		integer							:: 	n, ki, nSize, kSize
		!
		nSize	= size(Velo,3)
		kSize	= size(Velo,4)
		allocate(	f(3,nSize )		)
		p1 = 0.0_dp
		if(		kSize /= size(En,2)		) then
			write(*,*)"[calcFirstOrdP]: WARNING Energy and connection live on different k meshes!"
		end if
		!
		!
		do n = 1, nSize
			f 	= dcmplx(0.0_dp)
			pn	= 0.0_dp
			!FILL INTEGRATION ARRAY
			do ki = 1, kSize
				!PHASE SPACE DENSITY CORRECTION
				densCorr	= 0.5_dp * dot_product(		Fcurv(:,ki,n), Aconn(:,ki,n,n) 	)		* Bext
				f(:,ki)		= f(:,ki) + densCorr
				!POSITIONAL SHIFT
				call calcFmat(n,ki,Velo,En, Fmat)
				f(:,ki)	= f(:,ki) + matmul(Fmat, Bext) 
			end do
			!INTEGRATE
			do ki = 1, kSize
				pn = pn + f(:,ki)
			end do
			!SUM OVER n
			p1 = p1 + pn 
		end do
		!
		!
		!NORMALIZE
		p1 = p1 / kSize !	?!
		!
		return
	end subroutine


!privat
	subroutine	calcFmat(n0,ki, Velo ,En, Fmat)
		!calculates the linear response F matrix for magnetic field perturbation
		!F is derived in the semiclassical wavepacket approach (again see Niu PRL 112, 166601 (2014))
		integer,		intent(in)		:: n0, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  !V(3,nWfs,nWfs,nK)
		real(dp),		intent(in)		:: En(:,:)			!En(nK nWfs)
		real(dp),		intent(out)		:: Fmat(:,:)
		integer							:: i
		!
		Fmat = 0.0_dp		
		call addF2(n0, ki, Velo, En, Fmat)
		call addF3(n0, ki, Velo, En, Fmat)
		!
		return
	end subroutine



	!subroutine	addF1(ki, Fmat)
	!	integer,		intent(in)		:: ki
	!	complex(dp),	intent(inout)	:: Fmat(:,:)
	!	
	!	return
	!end


	subroutine	addF2(n0,ki, Velo ,En, Fmat)
		!
		!	F^(2)_ij = + Re \sum_{n/=0,m/=0} \eps_{j,k,l} * (V^k_nm V^l_m0 V^i_mn) / ( (E0-En)**2 (E0-Em) )
		!
		integer,		intent(in)		:: n0, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  !V(3,nK,nWfs,nWfs)
		real(dp),		intent(in)		:: En(:,:)			!	En(	nK, nWfs)	
		real(dp),		intent(out)		:: Fmat(:,:)
		complex(dp)						:: Vtmp
		real(dp)						:: eDiff
		integer							:: i, j, k, l, n,m, nSize, kSize
		!
		nSize	= size(Velo,3)
		kSize	= size(Velo,2)
		!loop spacial indices
		do j = 1, 3
			do i = 1, 3
				do k = 1, 3
					do l = 1, 3
						!loop bands
						do n = 1, nSize
							do m = 1, nSize
								if( n/=n0 .and. m/=n0) then
									!VELOCITIES
									Vtmp		= Velo(k,ki,n,m) * Velo(l,ki,m,n0) * Velo(i,ki,n0,n) 
									!ENERGIES
									eDiff		= ( 	En(ki,n0) - En(ki,n)	 )**2 	* 	 ( 	En(ki,n0) - En(ki,m)	)
									!MATRIX
									Fmat(i,j) 	= Fmat(i,j) +  myLeviCivita(j,k,l) * dreal(		Vtmp / dcmplx(eDiff)	)
								end if
							end do
						end do
						!
						!
					end do
				end do
			end do
		end do
		!
		!
		return
	end subroutine




	subroutine	addF3(n0,ki, Velo ,En, Fmat)	
		!
		!	F^(2)_ij = +- Re \sum_{n/=0} \eps_{j,k,l}  * (v^k_0 V^l_n0 V^i_0n) / ( (E0-En)**3  )
		!
		integer,		intent(in)		:: n0, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  	!Velo(		3		,	nK		,nWfs, nwFs)
		real(dp),		intent(in)		:: En(:,:)			!En(	nK	,	nWfs)
		real(dp),		intent(out)		:: Fmat(:,:)
		complex(dp)						:: Vtmp
		real(dp)						:: eDiff
		integer							:: i, j, k, l, n, nSize
		!
		nSize 	=	size(Velo,3)
		!loop spacial indices
		do j = 1, 3
			do i = 1, 3
				do k = 1, 3
					do l = 1,3
						!loop bands
						do n = 1, nSize
							if( n/=n0 ) then
								!VELOCITIES
								Vtmp		= Velo(k,ki,n0,n0) * Velo(l,ki,n,n0) * Velo(i,ki,n0,n) 
								!ENERGIES
								eDiff		= ( 	En(ki,n0) - En(ki,n)	 )**3 	
								!MATRIX
								Fmat(i,j) 	= Fmat(i,j) -  myLeviCivita(j,k,l) * dreal(		Vtmp / dcmplx(eDiff)	)
							end if
						end do
						!
						!
					end do
				end do
			end do
		end do
		!
		!
		return
	end subroutine





end module semiClassics