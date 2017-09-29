module peierls
	use mathematics,	only:	dp, PI_dp, i_dp, myExp, crossP, nIntegrate, eigSolver
	use sysPara
	use wannier,		only:	genKham, genUnkW
	use gaugeTrafo,		only:	calcConnCurv, testIfReal
	use	polarization,	only:	calcPolViaA
	use blochWf,		only:	 genUnk
	implicit none
	

	private
	public ::	peierlsSub


	contains







!public:
	subroutine	peierlsMethod(tHopp, pPei)
		complex(dp),	intent(in)		:: tHopp(:,:,:)	! tHopp(nWfs,nWfs,nSC)
		real(dp),		intent(out)		:: pPei(3)
		complex(dp),	allocatable		:: Ham(:,:)
		real(dp),		allocatable		:: EnPei(:,:)
		integer							:: qi
		!
		allocate( 	Ham(nWfs,nWfs)		)
		allocate(	EnPei(nWfs,nQ)		)
		!
		!
		call peierlsSub(tHopp)

		do qi = 1, nQ
			call genKham(qi, tHopp, Ham)
			call eigSolver(Ham, EnPei(:,qi))
			!ToDo
			!	calc bwfs, or unks
			!	calc connection
		end do
		!call calcPolViaA()


		return
	end subroutine


	subroutine peierlsSub(tHopp)

		return
	end subroutine

	!subroutine peierlsSub(wnF, unk, Ah, Fh, Vh, EnH, pPei )
	!	complex(dp),	intent(inout)	:: wnF(:,:,:)
	!	complex(dp),	intent(out)		:: unk(:,:,:), Ah(:,:,:,:), Fh(:,:,:,:), Vh(:,:,:,:)
	!	real(dp),		intent(out)		:: EnH(:,:), pPei(3)
	!	!
	!	!SHIFT WANNIER FUNCTIONS & GENERATE UNKs from shifted WANN
	!	call shiftWannF(wnF)
	!	call genUnkW(wnF, unk)
	!	!
	!	!CALC CONNECTION & CHECK IF REAL
	!	call calcConnCurv(unk, wnF, Ah, Fh, Vh, EnH)
	!	call testIfReal(Ah, Fh)
	!	!
	!	!INTEGRATE CONNECTION TO GET POL
	!	call calcPolViaA(dreal(Ah), pPei)
	!	!
	!	!
	!	return
	!end subroutine
!

	subroutine 	shiftWannF(wnF)
		!>shifts the wannier functions as following
		!!	wnf'(r)	= exp[i*(q/hbar)*\int^r_R A(r') \cdot dr'] * wnf(r)
		!
		complex(dp),	intent(inout)	:: wnF(:,:,:)
		real(dp)						:: aInt		
		integer							:: n, sc, ri
		!
		do n = 1, nWfs
			do sc = 1, nSC
				do ri = 1, nR
					aInt			= 1.0_dp * integrateA( Rcell(:,sc), rpts(:,ri) )	!1.0_dp = e/hbar in a.u.
					wnF(ri,sc,n) 	= wnF(ri,sc,n) * myExp(aInt)
				end do
			end do
		end do
		!
		!
		return
	end subroutine




real(dp) function integrateA(rMin, rMax)
	!integrates the vector potential A(r) analytically
	!	A(r)	= -0.5 cross_p[r,B]
	!	return 	= Integrate^\vec{rMax]_\vec{rMin} 		\vec{A(r')}.\vec{dr'}
	!			= -0.5 cross_p[rMin,B].(rMax - rMin)
	real(dp),		intent(in)		:: rMin(2), rMax(2)
	real(dp)						:: rU(3), rL(3)
	rU(1:2)		= rMax(1:2)
	rU(3)		= 0.0_dp
	rL(1:2)		= rMin(1:2)
	rL(3)		= 0.0_dp
	!
	integrateA	= -0.5_dp * dot_product( crossP(rL,Bext)	, (rU-rL)	)
	!
	return
end function



!privat:
	subroutine calcVecPot(B, Amag)
		!calculates the vector potential for a uniform magnetic field B
		!
		!	A = -0.5 cross_product( r, B )
		!
		real(dp),		intent(in)		:: B(3) 
		real(dp),		intent(out)		:: Amag(:,:)
		real(dp)						:: rpt(3)
		integer							:: ri
		!
		rpt(3)	= 0.0_dp
		do ri = 1, nR
			rpt(1:2)	= rpts(:,ri) 
			Amag(:,ri)	= crossP(	rpt(:) , B(:)	)
		end do
		!
		Amag = -0.5_dp * Amag
		!
		!
		return
	end subroutine


	subroutine peierlSubst(Amag, Htb)
		!	applies Peierls substitution of the hoping elements Htb due to vector field Amag
		!
		!
		!
		real(dp),		intent(in)		:: Amag(:,:) 	!Amag(3,nR)
		complex(dp),	intent(inout)	:: Htb(:,:,:)		!Htb(nR,nWfs,nWfs)
		integer							:: R, n, m
		real(dp)						:: aInt(3)


		do R = 1, nSC
			!CALC CORRECTION FACTOR 

			!UPDATE FOR ALL STATES
			do m = 1, nWfs
				do n = 1, nWfs
					Htb(R,n,m)	= Htb(R,n,m) !* myExp(aInt)
				end do
			end do
		end do
		!
		!
		return
	end subroutine


	subroutine calcKspaceHam(qi, Htb, Hk)
		!Fourier Trafo of tight binding Hamiltonian to k space
		integer,		intent(in)		:: qi
		complex(dp),	intent(in)		:: Htb(:,:,:)
		complex(dp),	intent(out)		:: Hk(:,:)
		integer							:: R, n, m
		real(dp)						:: cellP
		!
		do R = 1, nSC
			cellP	= dot_product( qpts(:,qi), Rcell(:,R) )
			do m = 1, nWfs
				do n = 1, nWfs
					Hk(n,m)	= Hk(n,m) + myExp(cellP) * Htb(R,n,m) 
				end do
			end do
		end do
		!
		!
		return
	end subroutine




	!subroutine peierlsPol(wnF,pPei)
	!	!	via the Wannier functions a tight binding hamiltonian is set up
	!	!	then peierls substitution is applied to the hopping parameters t_ij
	!	!
	!	!	t_ij -> t_ij * exp{	i e \integrate_{ij}	Aext dr}
	!	!
	!	!	see Ibanez-Azpiroz & Modugno 	PRA 90, 033609 (2014)
	!	!
	!	!	after solving the Ham, the connection & polarization is calculated from the new unks
	!	complex(dp),	intent(in)		:: wnF(:,:,:)		!wnF(nR, nSC, nWfs	)		
	!	real(dp),		intent(out)		:: pPei(3)
	!	real(dp),		allocatable		:: Amag(:,:), Epei(:,:), Aconn(:,:,:,:)
	!	complex(dp),	allocatable		:: Htb(:,:,:), Hk(:,:), bwf(:,:), unk(:,:,:)
	!	integer							:: ki
	!	!
	!	allocate(	Amag(3,nR)				)
	!	allocate(	Epei(nK,nWfs)			)
	!	allocate(	Htb(nR,nWfs,nWfs)		)
	!	allocate(	Hk(nWfs,nWfs)			)
	!	allocate(	bwf(nR,nWfs)			)
	!	allocate(	unk(nR,nK,nWfs)			)	
	!	allocate(	Aconn(3,nK,nWfs,nWfs)	)	
	!	!
	!	!CALC VECTOR POTENTIAL OF EXT. MAG FIELD
	!	call calcVecPot(Bext, Amag)
	!	!
	!	!SET UP AND SOLVE TB MODEL AT EACH K PNT
	!	do ki = 1, nK
	!		!SET UP TB HAMILTONIAN
	!		call genTBham(ki, wnF, Htb)
	!		!PEIERL SUBSTITUTION
	!		call peierlSubst(Amag, Htb)
	!		!TRAFO HAM TO K SPACE
	!		call calcKspaceHam(ki, Htb, Hk)
	!		!SOLVE
	!		call eigSolver(Hk, Epei(ki,:))
	!		call genBlochWf(ki,Hk, bWf)
	!		call genUnk(ki, bWf, unk(:,ki,:))
	!	end do
	!	!
	!	!
	!	!CALC CONN & POL
	!	pPei = 0.0_dp
	!	!call calcConn(unk, nKx, nKy, Aconn)
	!	!call calcPolViaA(Aconn, pPei)
	!	!
	!	!write Epei to file for comparisson
!
!	!	!
!	!	return
	!end



end module peierls