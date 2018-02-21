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
		real(dp)						::	dGx, dGy, xL, xR, yL, yR
		complex(dp)						::	num1, num2, denom, prefact, integral, &
											num1a, num1b, num2b, num2a
		!
		!
		do j = 1, nGq(qLoc)
			do i = 1, nGq(qLoc)
				!CASE 1 (i==j) gives 0
				if( i /= j )	then
					dGx		= Gvec(1,j,qLoc) - Gvec(1,i,qLoc) 
					dGy		= Gvec(2,j,qLoc) - Gvec(2,i,qLoc) 
					!
					!for each atom/well
					do at = 1, nAt
						xL	=	atPos(1,at) - atR(1,at)
						xR	=	atPos(1,at)	+ atR(1,at)
						yL	=	atPos(2,at)	- atR(2,at)
						yR	=	atPos(2,at)	+ atR(2,at)
						!
						!
						!CASE 2 (dGx==0)
						if( abs(dGx) < machineP ) then
							prefact		=	dcmplx(		Gvec(1,j,qLoc)*(xL-xR) 		)
							denom		=	dcmplx(		2.0_dp * dGy**2				)
							num1		=	i_dp * 	myExp(dGy*yL) * ( 	2.0_dp * i_dp	+ 			(yL-yR)*dGy 		)
							num2		=			myExp(dGy*yR) * ( 	2.0_dp 			+ i_dp * 	(yL-yR)*dGy			)
							!
							!
							integral	=	prefact * ( num1 + num2 ) / denom
						!---------------------------------------------------------------------------------------------------------------
						!
						!
						!CASE 3 (dGy==0)
						else if(	abs(dGy) < machineP	) then
							prefact		=	dcmplx(		Gvec(2,j,qLoc)*(yL-yR)		)
							denom		=	dcmplx(		2.0_dp * dGx**2				)
							num1		=			myExp(dGx*xL) * ( 	2.0_dp			+ i_dp *	(xR-xL)*dGx			)
							num2		=	i_dp *	myExp(dGx*xR) * (	2.0_dp * i_dp	+			(xR-xL)*dGx			)
							!
							!
							integral	=	prefact * ( num1 + num2 ) / denom
						!---------------------------------------------------------------------------------------------------------------
						!
						!
						!CASE 4 (dGx/=0 and dGy/=0)
						else if(	abs(dGx) >= machineP	 .and. 		abs(dGy) >= machineP	) then
							prefact		=	dcmplx(		1.0_dp						)
							denom		=	dcmplx(		2.0_dp	* dGx**2 * dGy**2	)
							!
							!
							num1a		=	-1.0_dp * myExp(dGy*yL) * (			Gvec(2,j,qLoc) 	*	( -2.0_dp*i_dp + (xL-xR)*dGx )	*	dGy	 &	
																			+	Gvec(1,j,qLoc) 	*	(  2.0_dp*i_dp + (yL-yR)*dGy )	*	dGx	 )
							!-----------------------------
							num1b		=	+1.0_dp * myExp(dGy*yR) * (			Gvev(2,j,qLoc) 	*	( -2.0_dp*i_dp + (xL-xR)*dGx )	*	dGy	&
																			+	Gvec(1,j,qLoc)	*	(  2.0_dp*i_dp - (yL-yR)*dGy )	*	dGx	)
							!-----------------------------
							num2a		=	+1.0_dp * myExp(dGy*yL) * (			Gvec(2,j,qLoc)	*	( -2.0_dp*i_dp - (xL-xR)*dGx )	*	dGy &
																			+	Gvec(1,j,qLoc)	*	(  2.0_dp*i_dp + (yL-yR)*dGy )	*	dGx )
							!-----------------------------
							num2b		=	-1.0_dp * myExp(dGy*yR) * (			Gvec(2,j,qLoc)	*	( -2.0_dp*i_dp - (xL-xR)*dGx )	*	dGy &
																			+	Gvec(1,j,qLoc)	*	(  2.0_dp*i_dp - (yL-yR)*dGy )	*	dGx )
							!-----------------------------
							!
							!
							num1		=	myExp(dGx*xR)	* 	(	num1a + num1b )
							num2		=	myExp(dGx*xL)	*	(	num2a + num2b )
							!
							!
							integral	=	prefact * ( num1 + num2 ) / denom	
	
						!---------------------------------------------------------------------------------------------------------------
						!
						!
						!DEFAULT
						else
							integral 	= 	dcmplx(0.0_dp)
							denom		= 	dcmplx(1.0_dp)
							stop "[add_Zeeman]: reached forbidden default"
						end if
						!
						!
						Hmat(i,j)	= Hmat(i,j)		+ 	(alphaZee/vol) * integral
						!
						!
					end do
				end if
				!
				!
			end do
		end do


		return
	end subroutine




end module ham_Zeeman