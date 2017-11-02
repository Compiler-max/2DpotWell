module semiClassics
	!this module uses a semiclassic approach to calculate the first ordrer correction
	!	to the polariztion induced by a perturbive magnetic field
	! 	see Niu PRL 112, 166601 (2014)
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, myLeviCivita
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
		complex(dp),	intent(in)		::	Fcurv(:,:,:,:), Aconn(:,:,:,:)	!Fcurv(3,nWfs, nQ)
		complex(dp),	intent(inout)	:: 	Velo(:,:,:,:)		!	 Velo(3, nWfs,nWfs, nQ)	
		real(dp),		intent(in)		::	En(:,:)				!	En(			nWfs, nQ)						
		real(dp),		intent(out)		:: 	p1(3)
		complex(dp), 	allocatable		::	f(:,:)
		real(dp)						::	pn(3)
		real(dp)						:: 	Fmat(3,3)
		real(dp)						:: 	densCorr(3)
		integer							:: 	n, ki, nSize, kSize
		!
		nSize	= size(Velo,3)
		kSize	= size(Velo,4)
		allocate(	f(3,kSize )		)
		if(		kSize /= size(En,2)		) then
			write(*,*)"[calcFirstOrdP]: WARNING Energy and velocities live on different k meshes!"
		end if
		!
		!
		write(*,*)"[calcFirstOrdP]: start calculating P' via semiclassic approach"
		p1 = 0.0_dp
		do n = 1, nSize
			f 	= dcmplx(0.0_dp)
			!FILL INTEGRATION ARRAY
			do ki = 1, kSize
				!PHASE SPACE DENSITY CORRECTION
				densCorr	= 0.5_dp * dot_product(		dreal(Fcurv(:,n,n,ki)), dreal(Aconn(:,n,n,ki) )	)		* Bext
				!f(:,ki)		= f(:,ki) + densCorr
				if( norm2(densCorr) > acc ) then
					write(*,*)	"[calcFirstOrdP]: warning the densCorr is none zero, norm2(densCorr)",norm2(densCorr)
				end if
				!POSITIONAL SHIFT
				Fmat	= 0.0_dp
				call calcFmat(n,ki,Velo,En, Fmat)
				f(:,ki)	= f(:,ki) + matmul(Fmat, Bext) 
			end do
			!INTEGRATE over k-space
			pn	= 0.0_dp
			do ki = 1, kSize
				pn = pn + f(:,ki)  / kSize
			end do
			!SUM OVER n
			p1 = p1 + pn
		end do
		!
		!
		!NORMALIZE
		p1 = p1  !	?!
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
		!
		Fmat = 0.0_dp		
		call addF2(n0, ki, Velo, En, Fmat)
		call addF3(n0, ki, Velo, En, Fmat)
		!
		return
	end subroutine




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
									Vtmp		= Velo(k,n,m,ki) * Velo(l,m,n0,ki) * Velo(i,n0,n,ki) 
									!ENERGIES
									eDiff		= ( 	En(n0,ki) - En(n,ki)	 )**2 	* 	 ( 	En(n0,ki) - En(m,ki)	)
									!MATRIX
									Fmat(i,j) 	= Fmat(i,j) +  myLeviCivita(j,k,l) * dreal(	Vtmp ) / eDiff	
									if(abs(dimag(Vtmp)) > acc ) write(*,*)	"[addF2]: non vanishing imag part detected:",dimag(Vtmp)
									if( eDiff < machineP ) write(*,*) "[addF2]: warning vanishing ediff"
									!write(*,'(a,e10.3,a,e10.3)')"[addF2]: |Vtmp|=",abs(Vtmp), "eDiff=",eDiff
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
								Vtmp		= Velo(k,n0,n0,ki) * Velo(l,n,n0,ki) * Velo(i,n0,n,ki) 
								!ENERGIES
								eDiff		= ( 	En(n0,ki) - En(n,ki)	 )**3 	
								!MATRIX
								Fmat(i,j) 	= Fmat(i,j) -  myLeviCivita(j,k,l) * dreal(		Vtmp / dcmplx(eDiff)	)
								if(abs(dimag(Vtmp)) > acc ) write(*,*)	"[addF3]: non vanishing imag part detected",dimag(Vtmp)
								if( eDiff < machineP ) write(*,*) "[addF3]: warning vanishing ediff"
								!write(*,'(a,e10.3,a,e10.3)')"[addF3]: |Vtmp|=",abs(Vtmp), "eDiff=",eDiff
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