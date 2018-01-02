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
		real(dp)						::	pnF2(3), pnF3(3)
		real(dp)						:: 	F2(3,3), F3(3,3)
		real(dp)						:: 	densCorr(3)
		integer							:: 	n, ki, nSize, kSize
		!
		nSize	= size(Aconn,3)
		kSize	= size(Velo,4)
		!
		write(*,*)	"read ub energies at q=1"
		do n = 1, size(En,1)
			write(*,*)	En(n,2)
		end do
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
			pnF2	= 0.0_dp
			pnF3	= 0.0_dp
			!
			do ki = 1, kSize		
				!PHASE SPACE DENSITY CORRECTION
				densCorr	= 0.5_dp * dot_product(		dreal(Fcurv(:,n,n,ki)), dreal(Aconn(:,n,n,ki) )	)		* Bext
				!f(:,ki)		= f(:,ki) + densCorr
				if( norm2(densCorr) > acc ) then
					write(*,*)	"[calcFirstOrdP]: warning the densCorr is none zero, norm2(densCorr)=",norm2(densCorr)
				end if
				!POSITIONAL SHIFT
				F2	= 0.0_dp
				F3	= 0.0_dp
				call addF2(n,ki,Velo,En, F2)
				call addF3(n,ki,Velo,En, F3)
				!Integrate
				pnF2		= pnF2		+ matmul(F2, Bext) / real(kSize,dp)
				pnF3		= pnF3		+ matmul(F3, Bext) / real(kSize,dp)
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
			write(*,'(a,i5,a,e12.5,a,e12.5,a,e12.5,a)')	"[calcFirstOrdP]: pNiuF2(n=", n,") =(" ,pnF2(1),&
																		", ",pnF2(2),", ", pnF2(3),")."
			write(*,'(a,i5,a,e12.5,a,e12.5,a,e12.5,a)')	"[calcFirstOrdP]: pNiuF3(n=", n,") =(" , pnF3(1),&
																		", ",pnF3(2),", ", pnF3(3),")."
			!SUM OVER n
			pF2 = pF2 + pnF2
			pF3 = pF3 + pnF3
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
		real(dp),		intent(out)		:: Fmat(:,:)
		!
		Fmat = 0.0_dp		
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
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  
		real(dp),		intent(in)		:: En(:,:)			!
		real(dp),		intent(out)		:: Fmat(:,:)
		complex(dp)						:: Vtmp
		real(dp)						:: eDiff, eDiff1, eDiff2
		integer							:: i, j, k, l, n,m, nSize
		!
		nSize	= size(Velo,3)
		
		!loop bands
		do n = 1, nSize
			do m = 1, nSize
				if( n/=nZero .and. m/=nZero) then
					!loop matrix indices
					do j = 1, 3
						do i = 1, 3
							!loop levi civita
							do k = 1, 3
								do l = 1, 3				
									!VELOCITIES
									Vtmp		= Velo(k,n,m,ki) * Velo(l,m,nZero,ki) * Velo(i,nZero,n,ki) 
									!ENERGIES
									eDiff1		= 	( 	En(nZero,ki) - En(n,ki)		)**2 
									eDiff2		=  	( 	En(nZero,ki) - En(m,ki)		)
									eDiff		= 	eDiff1 * eDiff2	
									!MATRIX
									Fmat(i,j) 	= Fmat(i,j) +   myLeviCivita(j,k,l) *  dreal( Vtmp )  / eDiff	
									!if(abs(dimag(Vtmp)) > acc ) write(*,*)	"[addF2]: non vanishing imag part detected:",dimag(Vtmp)
									if( abs(eDiff) < machineP )  then
										write(*,*)	"[addF2]: warning for k point = ",ki
										write(*,'(a,i3,a,i3,a,i3,a,e14.6)') "[addF2]: warning vanishing n0=",nZero,"n=",n," m=",m," eDiff=",eDiff
										write(*,'(a,i3,a,i3,a,e14.6)')	"[addF2]: ( E(",nZero,")-E(",n,") )**2=", eDiff1
										write(*,'(a,i3,a,i3,a,e14.6)')	"[addF2]: ( E(",nZero,")-E(",m,") )   =", eDiff2
										write(*,*)	"[addF2]: E(nZero=",nZero,")=",En(nZero,ki)
										write(*,*)	"[addF2]: E(n=",n,")=",En(n,ki)
										write(*,*)	"[addF2]: E(m=",m,")=",En(m,ki)
									end if
									!write(*,'(a,e10.3,a,e10.3)')"[addF2]: |Vtmp|=",abs(Vtmp), "eDiff=",eDiff				
								end do
							end do
							!
						end do
					end do
					!
				end if
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
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  
		real(dp),		intent(in)		:: En(:,:)			
		real(dp),		intent(out)		:: Fmat(:,:)
		complex(dp)						:: Vtmp
		real(dp)						:: eDiff
		integer							:: i, j, k, l, n, nSize
		!
		nSize 	=	size(Velo,3)
		!loop bands
		do n = 1, nSize
			if( n/=nZero ) then
				!loop matrix indices
				do j = 1, 3
					do i = 1, 3
						!loop levi civita
						do k = 1, 3
							do l = 1,3				
								!VELOCITIES
								Vtmp		= Velo(k,nZero,nZero,ki) * Velo(l,n,nZero,ki) * Velo(i,nZero,n,ki) 
								!ENERGIES
								eDiff		= ( 	En(nZero,ki) - En(n,ki)	 )**3 	
								!MATRIX
								Fmat(i,j) 	= Fmat(i,j) + prefactF3 * myLeviCivita(j,k,l) *	 dreal( Vtmp ) / eDiff
								!if(abs(dimag(Vtmp)) > acc ) write(*,*)	"[addF3]: non vanishing imag part detected",dimag(Vtmp)
								if( abs(eDiff) < machineP ) write(*,*) "[addF3]: warning vanishing n0=",nZero,"n=",n," eDiff=",eDiff
								!write(*,'(a,e10.3,a,e10.3)')"[addF3]: |Vtmp|=",abs(Vtmp), "eDiff=",eDiff
							end do								!
						end do
						!
					end do
				end do
				!
			end if
		end do
		!
		!
		return
	end subroutine





end module semiClassics