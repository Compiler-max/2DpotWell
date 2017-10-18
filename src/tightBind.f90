module tightBind
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, machineP, myExp, myLeviCivita, eigSolver, nIntegrate, isUnit, isHermitian
	use sysPara
	use blochWf,		only:	genBwfVelo, testNormUNK							
	use output,			only:	printMat, writeEnAndUNK
	implicit none	
	
	private
	public ::					solveELECTRstruct



		

	


	contains
	

	subroutine solveELECTRstruct(unk, En)
		complex(dp),		intent(out)		:: unk(:,:,:)
		real(dp),			intent(out)		:: En(:,:)
		complex(dp),		allocatable		:: Hmat(:,:), tHopp(:,:,:,:,:)

		integer								:: qi, nOrb

		allocate(	Hmat(nOrb*nAt,nOrb*nAt)			)
		allocate(	tHopp(nOrb,nOrb,nAt,nAt,nSC)	)
		!call setUpHopping(tHopp)

		do qi = 1, nQ

			call setUpHam(qi, tHopp, Hmat)
			call eigSolver(Hmat,En(:,qi))
		
			call genBwfVelo(qi, Hmat, unk(:,:,qi))	!omp 
		end do


	end subroutine






	subroutine setUpHam(qi, tHopp, Hmat)
		integer,			intent(in)		:: qi
		complex(dp),		intent(in)		:: tHopp(:,:,:,:,:)
		complex(dp),		intent(out)		:: Hmat(:,:)
		integer								:: j, i, nu, mu, nOrb
		nOrb = 1
		Hmat = dcmplx(0.0_dp)
			do j = 1, nAt
				do i = 1, nAt
					do nu = 1, nOrb
						do mu = 1, nOrb

						end do
					end do
				end do
			end do

			!call addSS(Hmat)
			!call addSP(Hmat)
			!call addPP(Hmat)
		return
	end subroutine

















end module tightBind