module ham_Zeeman
	!add Zeeman like term to hamitlonian
	use mpi
	use omp_lib
	use util_math,	only:	dp, PI_dp,i_dp, machineP, myExp
	use util_sysPara				
	implicit none	
	!#include "mpif.h"

	private
	public ::			add_Zeeman			



		

	


	contains




	subroutine add_Zeeman(qLoc, Hmat)
		!adds the operator 
		!				H_alpha = alpha ( op(x) op(p_y)	-	op(y) op(p_x)	)
		!
		!			this introduces a zeeman like spliting into Landau levels (breaks time reversal)
		!			as a result of a coupling B cdot L
		!			where B is the field breaking the symm
		!			L is the orbital momentum inside the wells
		!
		!	to compute analytic integrals the following cases have to be discussed
		!	case		dGx 		dGy
		!		1			0			0
		!		2			0			1
		!		3			1			0
		!		4			1			1
		!---------------------------------------------------------------------------
		integer,		intent(in)		::	qLoc
		complex(dp),	intent(inout)	::	Hmat(:,:)
		integer							::	i, j
		real(dp)						::	dGx, dGy
		complex(dp)						::	nom1, nom2, denom, prefact, integral
		!
		!
		do j = 1, nGq(qLoc)
			do i = 1, nGq(qLoc)
				!CASE 1 (i==j) gives 0
				if( i /= j )	then
					dGx		= Gvec(1,j,qLoc) - Gvec(1,i,qLoc) 
					dGy		= Gvec(2,j,qLoc) - Gvec(2,i,qLoc) 
					

					!todo: loop wells, get xR, xL, yL, yR
					!
					!
					!CASE 2 (dGx==0)
					if( abs(dGx) < machineP ) then
						prefact		=	dcmplx(		Gvec(1,j,qLoc)*(xL-xR) 		)
						denom		=	dcmplx(		2.0_dp * dGy**2				)
						nom1		=	i_dp * 	myExp(dGy*yL) * ( 	2.0_dp*i_dp	 	+ 		(yL-yR)*dGy 		)
						nom2		=			myExp(dGy*yR) * ( 	2.0_dp 			+ i_dp* (yL-yR)*dGy			)

						!
						integral	=	prefact * ( nom1 + nom2 ) / denom
					!------------------------------------------------------------------
					!
					!
					!CASE 3 (dGy==0)
					else if(	abs(dGy) < machineP	) then

					!------------------------------------------------------------------
					!
					!
					!CASE 4 (dGx/=0 and dGy/=0)
					else if(	abs(dGx) >= machineP	 .and. 		abs(dGy) >= machineP	) then

					!------------------------------------------------------------------
					!
					!
					!DEFAULT
					else
						stop "[add_Zeeman]: reached forbidden default"
					end if

					Hmat(i,j)	= Hmat(i,j)		+ 	(alphaZee/vol) * nom/denom
 
				end if
			end do
		end do


		return
	end subroutine




end module ham_Zeeman