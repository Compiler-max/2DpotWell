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
	subroutine add_magHam(nG_qi, Gvec, Hmat)
		!adds oscillating magnetic field to hamiltonian via peierls sub
		!
		!
		!	H 	= 	B_z / (2 V q) i  2pi aX 	dGy	(1-exp[i aX dGx])(1-exp[i aY dGy])		/		(dGy*(4pi^2-aX^2*dGx^2))
		!
		!		=: H_prefact * integral
		!
		!	with
		!		->		H_prefact = B_z / (2 V q) i  2pi aX		= 		 i pi B_z aX / ( V q)
		!	
		!		->		integral = dGy	(1-exp[i aX dGx])(1-exp[i aY dGy])		/		(dGy*(4pi^2-aX^2*dGx^2))
		!
		!							= (1-exp[i aX dGx])(1-exp[i aY dGy])	/	(4pi^2-aX^2*dGx^2)
		!	units:
		!		[p]	= hbar / a0
		!		[B] = hbar / (e a0^2)
		!		[q]	= 1 / a0
		!		[A] = e [B] / [q] e hbar / (e a0) = hbar / a0 = [p]
		!
		!only contribtutions if i/=j, i.e. to off-diagonal terms
		! no contribution for dGx = 0
		integer,		intent(in)		::	nG_qi
		real(dp),		intent(in)		::	Gvec(:,:)
		complex(dp),	intent(inout)	::	Hmat(:,:)
		integer							::	i, j
		real(dp)						::	qX_Period, dGx, dGy
		complex(dp)						::	H_prefact, denom, numer
		!
		!period of oscillating B field
		qX_period	=	2.0_dp * PI_dp / aX 
		H_prefact 	= 	dcmplx( 0.0_dp, 		PI_dp *	Bext(3) * aX / (qX_period*vol)	)		!purely imaginary due to i in formula
		!
		!Loop elements
		do j = 1, nG_qi
			do i = 1, nG_qi
				!
				!
				!only off diagonal elements can contribute
				if( i/=j ) then
					dGx		= Gvec(1,j) - Gvec(1,i) 
					dGy		= Gvec(2,j) - Gvec(2,i) 
					!
					!
					if( abs(dGx) > machineP  .and. abs(dGy) > machineP ) then
						!
						numer		= (	myExp(aX*dGx) - 1.0_dp) * (	myExp(aY*dgY) - 1.0_dp)
						!
						!
						denom		= dcmplx( 4.0_dp*PI_dp - aX**2 * dGx**2	)
						!
						!
						Hmat(i,j)	= Hmat(i,j) + H_prefact * numer / denom
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