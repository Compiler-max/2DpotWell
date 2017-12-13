module semiClassics
	!this module uses a semiclassic approach to calculate the first ordrer correction
	!	to the polariztion induced by a perturbive magnetic field
	! 	see Niu PRL 112, 166601 (2014)
	use mathematics,	only:	dp, PI_dp, i_dp, machineP, acc, myExp, myLeviCivita
	use sysPara,		only:	Bext, prefactF3, atPos, nAt

	implicit none



	private
	public ::			calcFirstOrdP, calcFmat



	contains



!TODO CHECK INDEXING OF VELO ACONN FCURV AFTER THEY HAVE BEEN CALCULATED


!public


	subroutine	calcFirstOrdP(Fcurv, Aconn, Velo, En, pF2, pF3)
		!calculates the first order polarization p1 according to
		!	P'= -int_dk [0.5 (Curv.Velo)*B_ext + a']
		complex(dp),	intent(in)		::	Fcurv(:,:,:,:), Aconn(:,:,:,:), Velo(:,:,:,:)	!Fcurv(3,nWfs, nQ)
		real(dp),		intent(in)		::	En(:,:)				!	En(			nWfs, nQ)						
		real(dp),		intent(out)		:: 	pF2(3), pF3(3)
		complex(dp)						::	pnF2(3), pnF3(3)
		complex(dp)						:: 	F2(3,3), F3(3,3)
		real(dp)						:: 	densCorr(3)
		integer							:: 	n, ki, nSize, kSize
		!
		nSize	= size(Velo,3)
		kSize	= size(Velo,4)
		!
		!
		if(		kSize /= size(En,2)		) then
			write(*,*)"[calcFirstOrdP]: WARNING Energy and velocities live on different k meshes!"
		end if
		!
		!
		write(*,*)"[calcFirstOrdP]: start calculating P' via semiclassic approach"
		pF2 = 0.0_dp
		pF3 = 0.0_dp
		do n = 1, nSize
			pnF2	= dcmplx(0.0_dp)
			pnF3	= dcmplx(0.0_dp)
			!
			do ki = 1, kSize		
				!PHASE SPACE DENSITY CORRECTION
				densCorr	= 0.5_dp * dot_product(		dreal(Fcurv(:,n,n,ki)), dreal(Aconn(:,n,n,ki) )	)		* Bext
				!f(:,ki)		= f(:,ki) + densCorr
				if( norm2(densCorr) > acc ) then
					write(*,*)	"[calcFirstOrdP]: warning the densCorr is none zero, norm2(densCorr)=",norm2(densCorr)
				end if
				!POSITIONAL SHIFT
				F2	= dcmplx(0.0_dp)
				F3	= dcmplx(0.0_dp)
				call addF2(n,ki,Velo,En, F2)
				call addF3(n,ki,Velo,En, F3)
				!Integrate
				pnF2		= pnF2		+ matmul(F2, Bext) / dcmplx(kSize)
				pnF3		= pnF3		+ matmul(F3, Bext) / dcmplx(kSize)
			end do
		
			!!!SINGLE ATOM
			!if( nAt == 1 ) then
			!	pnF2(:)	= pnF2(:) - atPos(:,1)		!calc center w.r.t. atom center
			!	pnF3(:)	= pnF3(:) - atPos(:,1)
			!!DOUBLE ATOM
			!else if( nAt == 2 ) then
			!	if( mod(n,2)== 0 ) then
			!		pnF2(:)	= pnF2(:) - atPos(:,2)
			!		pnF3(:)	= pnF3(:) - atPos(:,2)
			!	else	
			!		pnF2(:)	= pnF2(:) - atPos(:,1)
			!		pnF3(:)	= pnF3(:) - atPos(:,1)
			!	end if
			!end if
			!write to standard out
			write(*,'(a,i5,a,e12.5,a,e12.5,a,e12.5,a)')	"[calcFirstOrdP]: pNiuF2(n=", n,") =(" ,dreal(pnF2(1)),&
																		", ",dreal(pnF2(2)),", ", dreal(pnF2(3)),")."
			write(*,'(a,i5,a,e12.5,a,e12.5,a,e12.5,a)')	"[calcFirstOrdP]: pNiuF3(n=", n,") =(" ,dreal(pnF3(1)),&
																		", ",dreal(pnF3(2)),", ", dreal(pnF3(3)),")."
			!SUM OVER n
			pF2 = pF2 + dreal(pnF2)
			pF3 = pF3 + dreal(pnF3)
			!DEBUG
			if( norm2(dimag(pnF2)) > acc	) write(*,*)"[calcFirstOrdP]: F2 found complex pol. contribution from band n=",n 
			if( norm2(dimag(pnF3)) > acc	) write(*,*)"[calcFirstOrdP]: F3 found complex pol. contribution from band n=",n 
		end do
		!
		!
		return
	end subroutine



	subroutine	calcFmat(nZero,ki, Velo ,En, Fmat)
		!calculates the linear response F matrix for magnetic field perturbation
		!F is derived in the semiclassical wavepacket approach (again see Niu PRL 112, 166601 (2014))
		integer,		intent(in)		:: nZero, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  !V(3,nWfs,nWfs,nK)
		real(dp),		intent(in)		:: En(:,:)			!En(nK nWfs)
		complex(dp),	intent(out)		:: Fmat(:,:)
		!
		Fmat = dcmplx(0.0_dp)		
		call addF2(nZero, ki, Velo, En, Fmat)
		call addF3(nZero, ki, Velo, En, Fmat)
		!
		return
	end subroutine



!privat
	subroutine	addF2(nZero,ki, Velo ,En, Fmat)
		!
		!	F^(2)_ij = + Re \sum_{n/=0,m/=0} \eps_{j,k,l} * (V^k_nm V^l_m0 V^i_mn) / ( (E0-En)**2 (E0-Em) )
		!
		integer,		intent(in)		:: nZero, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  !V(3,nK,nWfs,nWfs)
		real(dp),		intent(in)		:: En(:,:)			!	En(	nK, nWfs)	
		complex(dp),	intent(out)		:: Fmat(:,:)
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
								if( n/=nZero .and. m/=nZero) then
									!VELOCITIES
									Vtmp		= Velo(k,n,m,ki) * Velo(l,m,nZero,ki) * Velo(i,nZero,n,ki) 
									!ENERGIES
									eDiff		= ( 	En(nZero,ki) - En(n,ki)	 )**2 	* 	 ( 	En(nZero,ki) - En(m,ki)	)
									!MATRIX
									Fmat(i,j) 	= Fmat(i,j) +  myLeviCivita(j,k,l) * 	Vtmp  / dcmplx(eDiff)	
									!if(abs(dimag(Vtmp)) > acc ) write(*,*)	"[addF2]: non vanishing imag part detected:",dimag(Vtmp)
									if( abs(eDiff) < machineP ) write(*,*) "[addF2]: warning vanishing n0=",nZero,"n=",n," m=",m
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




	subroutine	addF3(nZero,ki, Velo ,En, Fmat)	
		!
		!	F^(2)_ij = +- Re \sum_{n/=0} \eps_{j,k,l}  * (v^k_0 V^l_nZero V^i_0n) / ( (E0-En)**3  )
		!
		integer,		intent(in)		:: nZero, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  	!Velo(		3		,	nK		,nWfs, nwFs)
		real(dp),		intent(in)		:: En(:,:)			!En(	nK	,	nWfs)
		complex(dp),	intent(out)		:: Fmat(:,:)
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
							if( n/=nZero ) then
								!VELOCITIES
								Vtmp		= Velo(k,nZero,nZero,ki) * Velo(l,n,nZero,ki) * Velo(i,nZero,n,ki) 
								!ENERGIES
								eDiff		= ( 	En(nZero,ki) - En(n,ki)	 )**3 	
								!MATRIX
								Fmat(i,j) 	= Fmat(i,j) + prefactF3 * myLeviCivita(j,k,l) *	Vtmp / dcmplx(eDiff)
								!if(abs(dimag(Vtmp)) > acc ) write(*,*)	"[addF3]: non vanishing imag part detected",dimag(Vtmp)
								if( abs(eDiff) < machineP ) write(*,*) "[addF3]: warning vanishing n0=",nZero,"n=",n
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