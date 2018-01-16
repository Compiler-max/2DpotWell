module semiClassics
	!this module uses a semiclassic approach to calculate the first ordrer correction
	!	to the polariztion induced by a perturbive magnetic field
	! 	see Niu PRL 112, 166601 (2014)
	use mathematics,	only:	dp, PI_dp, i_dp, machineP, acc, myExp, myLeviCivita
	use sysPara,		only:	Bext, prefactF3, atPos, nAt, nWfs
	use output,			only:	printMat
	implicit none



	private
	public ::			calcFirstOrdP, calcFmat



	contains



!TODO CHECK INDEXING OF VELO ACONN FCURV AFTER THEY HAVE BEEN CALCULATED


!public
	subroutine	calcFirstOrdP(Fcurv, Aconn, Velo, En, pF2, pF3)
		!calculates the first order polarization p1 according to
		!	P'= -int_dk [0.5 (Curv.Velo)*B_ext + a']
		complex(dp),	intent(in)		::	Fcurv(:,:,:,:), Aconn(:,:,:,:), Velo(:,:,:,:)	
		real(dp),		intent(in)		::	En(:,:)			
		real(dp),		intent(out)		:: 	pF2(3), pF3(3)
		!real(dp)						::	pnF2(3), pnF3(3)
		real(dp)						:: 	F2(3,3), F3(3,3), F2k(3,3), F3k(3,3)
		real(dp)						:: 	densCorr(3)
		integer							:: 	n, ki, kSize, i, j
		!
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
		write(*,*)"[calcFirstOrdP]: will use ",size(Velo,3)," states"

		!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED)  &
		!$OMP PRIVATE(n, ki, kSize, i, j, densCorr, F2, F2k, F3, F3k) &
		!$OMP REDUCTION(+: pF2, pF3)
		do n = 1, nWfs
			!
			F2 = 0.0_dp
			F3 = 0.0_dp
			!
			!GET RESPONSE MATRIX
			do ki = 1, kSize		
				!
				!PHASE SPACE DENSITY CORRECTION
				densCorr	= 0.5_dp * dot_product(		dreal(Fcurv(:,n,n,ki)), dreal(Aconn(:,n,n,ki) )	)		* Bext
				!f(:,ki)		= f(:,ki) + densCorr
				if( norm2(densCorr) > acc ) then
					write(*,*)	"[calcFirstOrdP]: warning the densCorr is none zero, norm2(densCorr)=",norm2(densCorr)
				end if
				!POSITIONAL SHIFT
				call getF2(n,ki,Velo,En, F2k)
				call getF3(n,ki,Velo,En, F3k)
				!sum over K
				F2 = F2 + F2k
				F3 = F3 + F3k
			end do
			!
			!NORMALIZE
			F2 = F2 / real(kSize,dp)
			F3 = F3  / real(kSize,dp)
			!
			write(*,*)	"Bext = ",Bext(1)," ",Bext(2)," ",Bext(3)
			!
			!APPLY FIELD 
			!do i = 1, 3
			!	do j = 1, 3
			!		pF2(i)	= pF2(i) + F2(i,j) * Bext(j)
			!		pF3(i)	= pF3(i) + F3(i,j) * Bext(j)
			!	end do
			!end do
			pF2 = matmul(F2,Bext) 
			pF3 = matmul(F3,Bext) 
			!
			!
			!write to standard out
			write(*,'(a,i5,a,e12.5,a,e12.5,a,e12.5,a)')	"[calcFirstOrdP]: pNiuF2(n=", n, ") =(" ,pF2(1), ", ", pF2(2), ", ", pF2(3),")."
			write(*,'(a,i5,a,e12.5,a,e12.5,a,e12.5,a)')	"[calcFirstOrdP]: pNiuF3(n=", n, ") =(" ,pF3(1), ", ", pF3(2), ", ", pF3(3),")."
		end do
		!$OMP END PARALLEL DO
		!
		write(*,'(a,e12.5,a,e12.5,a,e12.5,a)')	"[calcFirstOrdP]: pNiuF2 =(" ,pF2(1), ", ", pF2(2), ", ", pF2(3),")."
		write(*,'(a,e12.5,a,e12.5,a,e12.5,a)')	"[calcFirstOrdP]: pNiuF3 =(" ,pF3(1), ", ", pF3(2), ", ", pF3(3),")."
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
		real(dp)						:: F2(3,3), F3(3,3)
		!
		Fmat	=	0.0_dp		
		call getF2(nZero, ki, Velo, En, F2)
		call getF3(nZero, ki, Velo, En, F3)
		!
		!
		Fmat	=	F2 + F3
		!
		return
	end subroutine












!privat
	subroutine	getF2(nZero,ki, Velo ,En, F2)
		!
		!	F^(2)_ij = + Re \sum_{n/=0,m/=0} \eps_{j,k,l} * (V^k_nm V^l_m0 V^i_mn) / ( (E0-En)**2 (E0-Em) )
		!
		integer,		intent(in)		:: nZero, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  
		real(dp),		intent(in)		:: En(:,:)			!
		real(dp),		intent(out)		:: F2(3s,3)
		complex(dp)						:: Vtmp
		real(dp)						:: eDiff, eDiff1, eDiff2
		integer							:: i, j, k, l, n,m, nSize
		!
		nSize	=	size(Velo,3)
		F2		=	0.0_dp
		!loop bands
		do n = 1, nSize
			do m = 1, nSize
				if( n/=nZero .and. m/=nZero) then
					!
					!loop matrix indices
					do j = 1, 3
						do i = 1, 3
							!
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
									F2(i,j) 	= F2(i,j) +   real(myLeviCivita(j,k,l),dp) *  dreal( Vtmp )  / eDiff	
									!DEGENERATE WARNING
									if( abs(eDiff) < machineP )  then
										write(*,*)	"[addF2]: warning for k point = ",ki
										write(*,'(a,i3,a,i3,a,i3,a,e14.6)') "[addF2]: warning degenerate bands n0=",nZero,"n=",n," m=",m," eDiff=",eDiff
										write(*,'(a,i3,a,i3,a,e14.6)')	"[addF2]: ( E(",nZero,")-E(",n,") )**2=", eDiff1
										write(*,'(a,i3,a,i3,a,e14.6)')	"[addF2]: ( E(",nZero,")-E(",m,") )   =", eDiff2
										write(*,*)	"[addF2]: E(nZero=",nZero,")=",En(nZero,ki)
										write(*,*)	"[addF2]: E(n=",n,")=",En(n,ki)
										write(*,*)	"[addF2]: E(m=",m,")=",En(m,ki)
									end if	
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




	subroutine	getF3(nZero,ki, Velo ,En, F3)	
		!
		!	F^(2)_ij = +- Re \sum_{n/=0} \eps_{j,k,l}  * (v^k_0 V^l_nZero V^i_0n) / ( (E0-En)**3  )
		!
		integer,		intent(in)		:: nZero, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  
		real(dp),		intent(in)		:: En(:,:)			
		real(dp),		intent(out)		:: F3(3,3)
		complex(dp)						:: Vtmp
		real(dp)						:: eDiff
		integer							:: i, j, k, l, n, nSize
		!
		nSize 	=	size(Velo,3)
		F3		=	0.0_dp
		!loop bands
		do n = 1, nSize
			if( n/=nZero ) then
				!
				!loop matrix indices
				do j = 1, 3
					do i = 1, 3
						!
						!loop levi civita
						do k = 1, 3
							do l = 1,3				
								!VELOCITIES
								Vtmp		= Velo(k,nZero,nZero,ki) * Velo(l,n,nZero,ki) * Velo(i,nZero,n,ki) 
								!ENERGIES
								eDiff		= ( 	En(nZero,ki) - En(n,ki)	 )**3 	
								!MATRIX
								F3(i,j) 	= F3(i,j) + real(prefactF3 * myLeviCivita(j,k,l),dp) *	 dreal( Vtmp ) / eDiff
								!DEGENERATE WARNING
								if( abs(eDiff) < machineP ) write(*,*) "[addF3]: warning degenerate bands n0=",nZero,"n=",n," eDiff=",eDiff
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