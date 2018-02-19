module semiClassics
	!this module uses a semiclassic approach to calculate the first ordrer correction
	!	to the polariztion induced by a perturbive magnetic field
	! 	see Niu PRL 112, 166601 (2014)
	use omp_lib
	implicit none



	private
	public ::			calcFirstOrdP, calcFmat

	integer, 		parameter 	:: 	dp 				= kind(0.d0)
	real(dp), 		parameter	:: 	machineP 		= 1e-15_dp
	real(dp)					:: 	acc				= 1e-7_dp
	real(dp),		parameter	::	aUtoAngstrm 	= 0.52917721092_dp
	real(dp),		parameter 	::	elemCharge	 	= 1.6021766208 * 1e-19_dp  *1e+6_dp! mu Coulomb


	contains



!TODO CHECK INDEXING OF VELO ACONN FCURV AFTER THEY HAVE BEEN CALCULATED


!public
	subroutine	calcFirstOrdP(cell_vol, Bext, prefactF3, Fcurv, Aconn, Velo, En, pol_F2, pol_F3)
		!calculates the first order polarization p1 according to
		!	P'= -int_dk [0.5 (Curv.Velo)*B_ext + a']
		real(dp),		intent(in)		::	cell_vol, Bext(3), prefactF3, Aconn(:,:,:,:), En(:,:)		
		complex(dp),	intent(in)		::	Fcurv(:,:,:,:), Velo(:,:,:,:)			
		real(dp),		intent(out)		::  pol_F2(:,:), pol_F3(:,:)
		!real(dp)						::	pnF2(3), pnF3(3)
		real(dp)						:: 	F2(3,3), F3(3,3), F2k(3,3), F3k(3,3), sumF2(3), sumF3(3)
		real(dp)						:: 	densCorr(3), polQuantum
		integer							:: 	n, ki, kSize, ind
		!
		kSize		= size(Velo,4)
		polQuantum	= elemCharge / ( cell_vol * aUtoAngstrm **2)
		!
		!
		write(*,*)"[calcFirstOrdP]: start calculating P' via semiclassic approach"
		write(*,*)"[calcFirstOrdP]: will use ",size(Velo,3)," states"

		pol_F2 = 0.0_dp
		pol_F3 = 0.0_dp
		
		!$OMP PARALLEL DEFAULT(SHARED)  &
		!$OMP PRIVATE(n, ki, densCorr, F2, F2k, F3, F3k)
		!$OMP DO SCHEDULE(STATIC)
		do n = 1, size(pol_F2,2)
			F2 = 0.0_dp
			F3 = 0.0_dp
			!
			!GET RESPONSE MATRIX
			do ki = 1, kSize		
				!
				!PHASE SPACE DENSITY CORRECTION
				densCorr(1:3)	= 0.5_dp * dot_product(		dreal(Fcurv(1:3,n,n,ki)), Aconn(1:3,n,n,ki)	)		* Bext
				if( norm2(densCorr) > acc ) write(*,*)	"[calcFirstOrdP]: WARNING the densCorr is none zero, norm2(densCorr)=",norm2(densCorr)
				!
				!POSITIONAL SHIFT
				call getF2(n,ki,Velo,En, F2k)
				call getF3(prefactF3, n,ki,Velo,En, F3k)
				!sum over K
				F2 = F2 + F2k
				F3 = F3 + F3k
			end do
			!
			!NORMALIZE
			F2 = F2 / real(kSize,dp)
			F3 = F3  / real(kSize,dp)
		
			!APPLY MATRIX 
			pol_F2(:,n) = matmul(F2,Bext) 
			pol_F3(:,n) = matmul(F3,Bext) 
			!
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		!

		do ind = 1, 3
			sumF2(ind) 	= sum( 	mod(pol_F2(ind,:)*aUtoAngstrm, 	polQuantum)			)
			sumF3(ind)	= sum(	mod(pol_F3(ind,:)*aUtoAngstrm,	polQuantum)			)
		end do
		!
		!PRINT F2
		write(*,*)															"[calcFirstOrdP]: F2 matrix contribution:"
		write(*,*)															" #state | 		<r>[Å]			| 		p[mu C / Å]"
		do n = 1, size(pol_F2,2)	
			write(*,'(i3,a,e13.4,a,e13.4,a,e13.4,a,a,e13.4,a,e13.4,a)')		n," | ", pol_F2(1,n)*aUtoAngstrm,", ",pol_F2(2,n)*aUtoAngstrm, ", ", pol_F2(3,n)*aUtoAngstrm,&
																			" | ", " (",mod(pol_F2(1,n)*aUtoAngstrm,polQuantum) ,", ",mod(pol_F2(2,n)*aUtoAngstrm,polQuantum),")"
		end do
		write(*,'(a,e13.4,a,e13.4,a,e13.4,a)')								"sum | 						|	(", sumF2(1),", ",sumF2(2), ", ", sumF2(3),")."
		!
		!PRINT F3
		write(*,*)															"[calcFirstOrdP]: F3 matrix contribution:"
		write(*,*)															" #state | 		<r>[Å]			| 		p[mu C / Å]"
		do n = 1, size(pol_F3,2)	
			write(*,'(i3,a,e13.4,a,e13.4,a,e13.4,a,a,e13.4,a,e13.4,a)')		n," | ", pol_F3(1,n)*aUtoAngstrm,", ",pol_F3(2,n)*aUtoAngstrm, ", ", pol_F3(3,n)*aUtoAngstrm,&
																			" | ", " (",mod(pol_F3(1,n)*aUtoAngstrm,polQuantum) ,", ",mod(pol_F3(2,n)*aUtoAngstrm,polQuantum),")"
		end do
		write(*,'(a,e13.4,a,e13.4,a,e13.4,a)')								"sum | 						|	(", sumF3(1),", ",sumF3(2), ", ", sumF3(3),")."
		!
		!PRINT TOT
		write(*,*)															"[calcFirstOrdP] total first order pol:"
		write(*,'(a,e13.4,a,e13.4,a,e13.4,a)')								"p'= (",sumF2(1)+sumF3(1),", ",sumF2(2)+sumF3(2),", ",sumF2(3)+sumF3(3),") [mu C/ Å] "

		!
		!DEBUG
		if(	kSize /= size(En,2)	)	 write(*,*)"[calcFirstOrdP]: WARNING Energy and velocities live on different k meshes!"
		!
		!
		return
	end subroutine



	subroutine	calcFmat(prefactF3, nZero,ki, Velo ,En, Fmat)
		!calculates the linear response F matrix for magnetic field perturbation
		!F is derived in the semiclassical wavepacket approach (again see Niu PRL 112, 166601 (2014))
		integer,		intent(in)		:: nZero, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  !V(3,nWfs,nWfs,nK)
		real(dp),		intent(in)		:: prefactF3, En(:,:)			!En(nK nWfs)
		real(dp),		intent(out)		:: Fmat(3,3)
		real(dp)						:: F2(3,3), F3(3,3)
		!
		Fmat	=	0.0_dp		
		call getF2(nZero, ki, Velo, En, F2)
		call getF3(prefactF3, nZero, ki, Velo, En, F3)
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
		real(dp),		intent(out)		:: F2(3,3)
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
									!F2(i,j) 	= F2(i,j) +   myLeviCivita(j,k,l) *  dreal( Vtmp )  / eDiff	
									!DEGENERATE WARNING
									if( abs(eDiff) < machineP )  then
										write(*,*)	"[addF2]: WARNING for k point = ",ki
										write(*,'(a,i3,a,i3,a,i3,a,e14.6)') "[addF2]: WARNING degenerate bands n0=",nZero,"n=",n," m=",m," eDiff=",eDiff
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




	subroutine	getF3(prefactF3, nZero,ki, Velo ,En, F3)	
		!
		!	F^(2)_ij = +- Re \sum_{n/=0} \eps_{j,k,l}  * (v^k_0 V^l_nZero V^i_0n) / ( (E0-En)**3  )
		!
		integer,		intent(in)		:: nZero, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)  
		real(dp),		intent(in)		:: prefactF3, En(:,:)			
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
								F3(i,j) 	= F3(i,j) + real(prefactF3,dp) * real(myLeviCivita(j,k,l),dp) *	 dreal( Vtmp ) / eDiff
								!F3(i,j) 	= F3(i,j) + prefactF3 * myLeviCivita(j,k,l) *	 dreal( Vtmp ) / eDiff
								!DEGENERATE WARNING
								if( abs(eDiff) < machineP ) write(*,*) "[addF3]: WARNING degenerate bands n0=",nZero,"n=",n," eDiff=",eDiff
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



	!helper
	integer function myLeviCivita(i,j,k)
		!Hard coded Levi Civita tensor
		integer,		intent(in)		:: i,j,k
		logical							:: even, odd
		!
		!
		even	= (i==1 .and. j==2 .and. k==3) .or. (i==2 .and. j==3 .and. k==1) .or. (i==3 .and. j==1 .and. k==2)
		odd		= (i==3 .and. j==2 .and. k==1) .or. (i==1 .and. j==3 .and. k==2) .or. (i==2 .and. j==1 .and. k==3)
		!
		if(even) then
			myLeviCivita	=  1
		else if (odd) then
			myLeviCivita	= -1
		else
			myLeviCivita	=  0
		end if
		!
		!DEBUGGING
		if(even .and. odd) then
			write(*,*)"[myLeviCivita]: myLeviCivita detected even and odd, contact the programer he fucked up"
		end if
		!
		return
	end function



end module semiClassics