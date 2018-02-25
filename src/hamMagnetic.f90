module ham_Magnetic
	!add Zeeman like term to hamitlonian
	use mpi
	use omp_lib
	use util_math,		only:	dp, PI_dp, machineP, myExp
	use util_sysPara				
	implicit none	
	!#include "mpif.h"

	private
	public ::			add_magHam			



		

	


	contains





	!EXTERNAL MAGNETIC FIELD ( OSCILLATING )
	subroutine add_magHam(qLoc, Hmat)
		!Performs Peierls substitution on plane wave basis.
		!
		!neglects terms which are second order in the external field
		!
		!only contribtutions if i/=j, i.e. to off-diagonal terms
		! no contribution for dGx = 0
		integer,		intent(in)		::	qLoc
		complex(dp),	intent(inout)	::	Hmat(:,:)
		integer							::	i, j
		real(dp)						::	qX_Period, dGx, dGy
		complex(dp)						::	H_prefact
		!
		!period of oscillating B field
		qX_period	=	2.0_dp * PI_dp / aX 
		H_prefact 	= 	dcmplx(	0.5_dp * Bext(3) / qX_period	)
		!
		!Add to Hamiltonian
		do j = 1, nGq(qLoc)
			do i = 1, nGq(qLoc)
				!
				!
				!CASE 1: gives 0 
				if( i/=j ) then
					dGx		= Gvec(1,j,qLoc) - Gvec(1,i,qLoc) 
					dGy		= Gvec(2,j,qLoc) - Gvec(2,i,qLoc) 
					!
					!
					!
					!CASE 2: dGx == 0 gives zero
					if( abs(dGx) > machineP  ) then
						!
						!
						!
						!CASE 3: dGy == 0
						if( abs(dGy) < machineP ) then
							Hmat(i,j)	= Hmat(i,j) + H_prefact * h1_gyZero( dGx, qX_period)
						
						!-------------------------------------------------------------------------------------------------------------
						!
						!
						!
						!CASE 4: dGy /= 0
						else if( abs(dGy) >= machineP ) then
							Hmat(i,j)	= Hmat(i,j) + H_prefact * h1_gyFull( dGx, dGy, qX_period)	
						end if					

						!-------------------------------------------------------------------------------------------------------------								!						!DEFAULT						else							stop "[add_magHam]: reached forbidden region in dGx, dGy discussion"						end if 						!-------------------------------------------------------------------------------------------------------------
					end if
					!
					!-----------------------------------------------------------------------------------------------------------------
				end if
				!
				!---------------------------------------------------------------------------------------------------------------------
			end do
		end do
		!
		!DEBUG
		if( myID == root .and. qLoc == 1) 	write(*,'(a,i3,a)')	"[#",myID,";addMagHam]: hello there"		
		if( norm2(Bext(1:2)) > machineP ) 	write(*,'(a,i3,a)') "[#",myID,";addMagHam]: WARNING found non zero x or y component. External field should be along z direction! "
		!
		return
	end subroutine



	complex(dp)	function h1_gyZero( dGx, qX )
		!matrix element for dGy==0
		real(dp),		intent(in)		:: dGx, qX
		!
		h1_gyZero	= 	qX 		*	( myExp(dGx*aX) - 1.0_dp )		*	aY
		!
		return
	end function

	complex(dp)	function h1_gyFull( dGx, dGy, qX)
		!matrix element for dGy/=0
		real(dp),		intent(in)		:: dGx, dGy, qX
		!
		h1_gyFull	=	qX		*	( myExp(dGx*aX) - 1.0_dp )		* 	( myExp(dGy*aY) - 1.0_dp )
		!
		return
	end function



end module ham_Magnetic