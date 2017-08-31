module peierls
	use mathematics,	only:	dp, PI_dp, i_dp, myExp, crossP, nIntegrate, eigSolver
	use sysPara, 		only: 	readInp, insideAt, getKindex, getRindex, &
									dim, aX, aY,vol, nAt, atR, atPos, atPot,&
									nG, nG0, Gcut, &
									nK, nKx, nKy, nKw, nKxW, nKyW,  nWfs, nSC, nSCx, nSCy, dkx, dky, &
									nR, nRx, nRy, R0, dx, dy, &
									Gvec, atPos, atR, kpts, kptsW, rpts,Rcell, gaugeSwitch, trialOrbSw, trialOrbVAL, Zion,&
									Bext
	use wannier,		only:	genTBham
	use berryWF,		only:	calcConn, calcPolViaA
	implicit none
	

	private
	public ::	peierlsPol


	contains







!public:
	subroutine peierlsPol(wnF,pPei)
		complex(dp),	intent(in)		:: wnF(:,:,:)		!wnF(nR, nSC, nWfs	)		
		real(dp),		intent(out)		:: pPei(3)
		real(dp),		allocatable		:: Amag(:,:)
		real(dp),		allocatable		:: Epei(:,:)
		complex(dp),	allocatable		:: Htb(:,:)
		integer							:: ki

		allocate(	Amag(3,nR)		)
		allocate(	Epei(nK,nWfs)	)
		allocate(	Htb(nWfs,nWfs)	)

		!CALC VECTOR POTENTIAL OF EXT. MAG FIELD
		call calcVecPot(Bext, Amag)


		do ki = 1, nK


			!SET UP TB HAMILTONIAN
			call genTBham(ki, wnF, Htb)
			!PEIERL SUBSTITUTION
			call peierlSubst(Amag, Htb)
			!SOLVE
			call eigSolver(Htb, Epei(ki,:))
		end do



		!CALC CONN & POL
		pPei = 0.0_dp
		!call calcConn(unk, nKx, nKy, Aconn)
		!call calcPolViaA(Aconn, pPei)
		!
		!
		return
	end







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
		!
		return
	end


	subroutine peierlSubst(Amag, Htb)
		real(dp),		intent(in)		:: Amag(:,:) 	!Amag(3,nR)
		complex(dp),	intent(inout)	:: Htb(:,:)		!Htb(nWfs,nWfs)


		return
	end


end module peierls