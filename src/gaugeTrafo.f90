module gaugeTrafo
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, myLeviCivita, nIntegrate, eigSolver, rotMat, myCommutat
	use sysPara
	use berry,			only:	calcWaveMat
	use wannier,		only:	genUnkW, calcWannMat
	!use 
	implicit none

	private
	public ::	calcConnCurv, testIfReal

	contains

!public
	subroutine calcConnCurv(unk, wnF, Ah, Fh, Vh, EnH)
		!uses the Wannier functions to generate the four matrix elements Hw, Hwa, Aw, Fw
		complex(dp),	intent(in)		:: unk(:,:,:), wnF(:,:,:)
		complex(dp),	intent(out)		:: Ah(:,:,:,:), Fh(:,:,:,:), Vh(:,:,:,:)
		real(dp),		intent(out)		:: EnH(:,:)
		complex(dp),	allocatable		:: Hw(:,:,:), Hwa(:,:,:,:), Aw(:,:,:,:), Fw(:,:,:,:,:), Fhtens(:,:,:,:,:)
		integer							:: n,m,ki,a,b,c
		real							:: T0,T1,T, Pt0, Pt1, Pt
		!
		allocate(	Hw(		nK, nWfs, nWfs)			)
		allocate(	Hwa(	3, nK, nWfs, nWfs)		)
		allocate(	Aw(		3, nK, nWfs, nWfs)		)
		allocate(	Fw(		3, 3, nK, nWfs, nWfs)	)
		allocate(FhTens(	3, 3, nK, nWfs, nWfs)	)
		!

		!GET DESIRED MATRICES VIA K SPACE (calcWaveMat) OR VIA R SPACE (calcWannMat)
		call cpu_time(T0)
		select case(connSwitch)
			case(0)
				write(*,*)"[calcConnCurv]: via K space"
				call calcWaveMat(unk, Hw, Hwa, Aw, Fw)
			case(1)
				write(*,*)"[calcConnCurv]: via R space"
				call calcWannMat(unk, Hw, Hwa, Aw, Fw)
			case default
				write(*,*)"[calcConnCurv]: via K space"
				call calcWaveMat(unk, Hw, Hwa, Aw, Fw)
		end select
		call cpu_time(T1) 
		T	= T1 -T0
		write(*,'(a,f12.8,a)')"[calcConnCurv]: the 4 matrices on coarse mesh where calculated in ",T," seconds"
		!
		!
		call cpu_time(Pt0)

		!select case(gaugeBack)
		!	case(0)
		!		call copyMat(Hw, Hwa, Aw, Fw, Ah, FhTens, Vh, EnH)
		!	case(1)
!
!		!	case(default)
!
		!end select
		call shiftToHamGauge(Hw, Hwa, Aw, Fw, Ah, FhTens, Vh, EnH)
		call cpu_time(Pt1)
		Pt	= Pt1-Pt0
		write(*,'(a,f12.8,a)')"[calcConnCurv]: the conversion to Ham gauge took ",Pt," seconds"
		!
		!CONVERT OMEGA TENSOR TO VECTOR
		do m = 1, nWfs
			do n = 1, nWfs
				do ki = 1, nK
					!EVAL CROSS PRODUCT
					do c = 1,3
						do b= 1,3
							do a=1,3
								if( myLeviCivita(a,b,c) /= 0) then
									Fh(c,ki,n,m) = Fh(c,ki,n,m) + myLeviCivita(a,b,c) * FhTens(a,b,ki,n,m)
								end if
							end do 
						end do
					end do
					!
				end do
			end do
		end do
		!
		!
		return
	end subroutine


	subroutine copyMat(Hw, Hwa, Aw, Fw, Ah, FhTens, Vh, EnH)
		complex(dp),	intent(in)		:: Hw(:,:,:), Hwa(:,:,:,:), Aw(:,:,:,:), Fw(:,:,:,:,:)
		complex(dp),	intent(out)		:: Ah(:,:,:,:), Fhtens(:,:,:,:,:), Vh(:,:,:,:), EnH(:,:)

		!todo
		return
	end


	subroutine testIfReal(Ah, Fh)
		complex(dp),	intent(in)		:: Ah(:,:,:,:), Fh(:,:,:,:)
		integer							:: a, ki, n, m
		integer							:: cntA, cntF, tot
		real(dp)						:: avgA, avgF
		!
		cntA	= 0
		cntF	= 0
		tot 	= 0
		avgA	= 0.0_dp
		avgF 	= 0.0_dp
		!
		do m = 1, nWfs
			do n = 1, nWfs
				do ki = 1, nK
					do a = 1, 3
						!check connection
						avgA	= avgA + abs(dimag(Ah(a,ki,n,m)))
						if( abs(dimag(Ah(a,ki,n,m))) > acc	) then
							cntA = cntA +1
						end if 

						!check curvature
						avgF	= avgF + abs(dimag(Fh(a,ki,n,m)))
						if( abs(dimag(Fh(a,ki,n,m))) > acc	) then
							cntF = cntF + 1
						end if
						tot = tot + 1
					end do
				end do
			end do
		end do

		avgA	= avgA / tot
		avgF	= avgF / tot
		write(*,'(a,i8,a,i8,a,f15.12,a,f15.12)')	"[testIfReal]: ",cntA, " of ",tot," = ",cntA/real(tot,dp),&
														"% checked conn. vals. are not strictly real, average imag.=",avgA
		write(*,'(a,i8,a,i8,a,f15.12,a,f15.12)')	"[testIfReal]: ",cntF, " of ",tot," = ",cntF/real(tot,dp),&
														"% checked curv. vals. are not strictly real, average imag.=",avgF												
		return
	end	subroutine







!privat
	subroutine shiftToHamGauge(Hw, Hwa, Aw, Fw, Ah, FhTens, Vh, EnH)
		!Shift given matrices to Ham Gauge
		!	see Wang/Vanderbilt PRB 74, 195118 (2006)
		!
		!	ToDo: calc Fmat
		complex(dp),	intent(in)		:: Hw(:,:,:)		!Hw(nKw, nWfs, nWfs)
		complex(dp),	intent(in)		:: Hwa(:,:,:,:)		!Hwa(3,nKw, nWfs, nWfs)
		complex(dp),	intent(in)		:: Aw(:,:,:,:)		!Aw(3,nKw, nWfs, nWfs)
		complex(dp),	intent(in)		:: Fw(:,:,:,:,:)	!Fw(3,3,nKw,nWfs,nWfs)
		complex(dp),	intent(out)		:: Ah(:,:,:,:), FhTens(:,:,:,:,:), Vh(:,:,:,:)
		real(dp),		intent(out)		:: EnH(:,:)
		complex(dp),	allocatable		:: Dh(:,:,:), U(:,:), Ahbar(:,:,:)
		integer							:: ki, i
		!
		allocate( 	Dh(		3, nWfs, nWfs	)		)
		allocate(	U(		nWfs, nWfs		)		)
		allocate(	Ahbar(	3, nWfs, nWfs	)		)
		EnH 	= 0.0_dp
		FhTens	= dcmplx(0.0_dp)
		!
		do ki = 1, nK
			call HamSolve(ki, Hw, 	U, EnH)
			do i = 1, 3
				call rotMat(	U	, Aw(i,ki,:,:), Ahbar(i,:,:)	)
			end do
			call calcDmat(ki, EnH, U, Hwa, 	Dh)
			call calcAmat(ki, Ahbar, Dh, 	Ah)
			call calcFmat(ki, Fw, U, Ahbar, Dh, 	FhTens)
			call calcVmat(ki, Hwa, U, Ahbar, EnH,	Vh)
		end do
		!
		!
		return
	end subroutine







	subroutine calcFmat(ki, Fw, U, Ahbar, Dh , FhTens)
		!	the curvature is genarted by summing up 4 contributions (PRB 74, 195118 (2006))
		!
		!	F^H_ab = 	Fbar^H_ab				(1)
		!				-	[Dh_a, Ahbar_b]		(2)
		!				+	[Dh_b, Ahbar_a]		(3)
		!				- i	[Dh_a, Dh_b]		(4)
		!
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: Fw(:,:,:,:,:), U(:,:), Ahbar(:,:,:), Dh(:,:,:)		!Fw(3,3,nWfs,nWfs)
		complex(dp),	intent(out)		:: FhTens(:,:,:,:,:)
		complex(dp),	allocatable		:: Matcomm(:,:)
		integer							:: a,b
		!
		allocate(	Matcomm(nWfs,nWfs)	)
		!
		do b = 1, 3
			do a = 1, 3
				!(1):
				call rotMat(U, Fw(a,b,ki,:,:), FhTens(a,b,ki,:,:)	)
				!
				!(2):
				call myCommutat(	Dh(a,:,:), Ahbar(b,:,:)	, Matcomm(:,:))
				FhTens(a,b,ki,:,:)	= FhTens(a,b,ki,:,:) - Matcomm(:,:)
				!
				!(3):
				call myCommutat(	Dh(b,:,:), Ahbar(a,:,:), Matcomm(:,:))
				FhTens(a,b,ki,:,:)	= FhTens(a,b,ki,:,:) + Matcomm(:,:)
				!
				!(4):
				call myCommutat(	Dh(a,:,:), Dh(b,:,:),	Matcomm(:,:)	)
				FhTens(a,b,ki,:,:)	= FhTens(a,b,ki,:,:) - i_dp * Matcomm(:,:)
			end do
		end do
		!
		return
	end subroutine






	subroutine HamSolve(ki, Hw, U, En)
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: Hw(:,:,:) 		!Hw(nKw, nWfs, nWfs)
		complex(dp),	intent(out)		:: U(:,:)
		real(dp),		intent(out)		:: En(:,:)			!EnH(ki, nWfs)
		!
		!COPY Hw to U 							
		if (	size(U,1)==size(Hw,2) .and. size(U,2)==size(Hw,3)	) then
			U = Hw(ki,:,:)
		else
			U = dcmplx(0.0_dp)
			write(*,*)"[gaugeTrafo/HamSolve]: ranks of matrix do not match"
		end if
		!
		!SOLVE EIGVALUE PROBLEM (using U)
		call eigSolver(U, En(ki,:))
		!
		return
	end subroutine


	subroutine calcDmat(ki, EnH, U, Hwa, Dh)
		integer,		intent(in)		:: ki
		real(dp),		intent(in)		:: EnH(:,:)
		complex(dp),	intent(in)		:: U(:,:), Hwa(:,:,:,:)
		complex(dp),	intent(out)		:: Dh(:,:,:)
		complex(dp),	allocatable		:: HaBar(:,:,:)
		integer							:: n, m, i
		!
		!ROT H^W_a to Ham gauge
		allocate( 	HaBar(3,nWfs,nWfs)	)
		do i = 1, 3
			call rotMat( U(:,:), Hwa(i,ki,:,:), Habar(i,:,:)	) 
		end do
		!use Kubo formula to calc Dh
		do m = 1, nWfs
			do n = 1, nWfs
				if( n==m  ) then
					Dh(:,n,m)	= dcmplx(0.0_dp)
				else
					do i = 1, 3
						Dh(i,n,m)	= HaBar(i,n,m)	/ dcmplx( EnH(ki,m) - EnH(ki,n)	)
					end do
				end if
			end do
		end do
		!
		!
		return
	end subroutine

	

	subroutine calcAmat(ki, Ahbar, Dh, Ah)		
		!	A^H_i 	= U' A^W_i U 	+ i D^H  
		!			= Ahbar			+ i D^H
		!
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: Ahbar(:,:,:), Dh(:,:,:)
		complex(dp),	intent(out)		:: Ah(:,:,:,:)
		integer							:: i
		!
		do i = 1, 3
			Ah(i,ki,:,:)	= Ahbar(i,:,:)	+	i_dp * Dh(i,:,:)
		end do
		return
	end subroutine


	subroutine calcVmat(ki, Hwa, U, Ahbar, EnH,	Vh)
		!calulates the velocity matrix (PRB 74, 195118 (2006): eq(31)	)
		!
		!	V^(H)_nm,a = Hbar^(H)_nm,a - i (E^(H)_m - E^(H)_n) Abar^(H)_nm,a
		!
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: Hwa(:,:,:,:), U(:,:), Ahbar(:,:,:)
		real(dp),		intent(in)		:: EnH(:,:)
		complex(dp),	intent(out)		:: Vh(:,:,:,:)
		complex(dp),	allocatable		:: HhbarA(:,:,:)
		integer							:: n, m, a
		!
		allocate( HhBarA(3,nWfs,nWfs))
		!
		!CALC H^(H)_a
		do a = 1,3
			call rotMat( U(:,:), Hwa(a,ki,:,:), HhBarA(a,:,:) )
		end do
		!
		!VELCOITIES
		do m = 1, nWfs
			do n = 1, nWfs
				do a = 1, 3
					Vh(a,ki,n,m) = HhBarA(a,n,m) - i_dp * ( EnH(ki,m)-EnH(ki,n) ) * Ahbar(a,n,m)
				end do
			end do
		end do
		!
		!
		return
	end subroutine









end module gaugeTrafo