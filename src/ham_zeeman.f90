module ham_Zeeman
	!add Zeeman like term to hamitlonian
	use util_math,		only:	dp, PI_dp,i_dp, machineP, myExp
	use util_sysPara				
	implicit none	
	!#include "mpif.h"

	private
	public ::			add_Zeeman			



		

	


	contains




	subroutine add_Zeeman(nG_qi, Gvec, Hmat)
		!adds the operator 
		!				H_alpha = alpha ( op(x) op(p_y)	-	op(y) op(p_x)	)
		!
		!			V_zeeman 	= 	- magMom \cdot B		  [hatreee]
		!			magMom		=	mu_bohr * g_L * L / hbar    [a.u.=	e hbar/(2 m_e)]
		!			g_L			= 	1.0
		!		units:
		!			[L]			=	hbar
		!			[B]			= 	hbar / (e a0^2)
		!			[magMom][B]	=	hbar^2 / ( m_e a0^2) = m_e c^2 alpha^2 = hartree
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
		integer,		intent(in)		::	nG_qi
		real(dp),		intent(in)		::	Gvec(:,:)
		complex(dp),	intent(inout)	::	Hmat(:,:)
		integer							::	i, j, at 
		real(dp)						::	dGx, dGy, Gjx, Gjy, x0, y0, at_rad, magMom, minI, maxI
		complex(dp)						::	ham, alphaZee, integral, &
											num, denom, prefact, numA, numB
																!UNITS:
		magMom		=	+ 0.5_dp 								!a.u.= e hbar/( 2 m_e)
		alphaZee	=	dcmplx(	-1.0_dp * magMom * Bext(3)	)	! [magMom] [Bext]
		!
		!
		do j = 1, nG_qi
			do i = 1, nG_qi
				!CASE 1 (i==j) gives 0
				if( i /= j )	then
					dGx		= Gvec(1,j) - Gvec(1,i) 
					dGy		= Gvec(2,j) - Gvec(2,i) 
					Gjx		= Gvec(1,j)
					Gjy		= Gvec(2,j)
					!
					!for each atom/well
					do at = 1, nAt
						x0		= atPos(1,at)
						y0		= atPos(2,at)
						at_rad 	= atR(1,at)
						if( atR(2,at) /= at_rad)		stop '[hamZeeman]: only quadratic wells are supported'
						!
						!
						!
						!CASE 2 (dGx==0)
						if( abs(dGx) < machineP ) then
							prefact		=	i_dp * dcmplx(4.0_dp * at_rad * Gjx												) * myExp(dGy*y0)
							num			=	dcmplx(		at_rad * dGy * cos(at_rad*dGy) 		- 		sin(at_rad*dGy)			)
							denom		=	dcmplx(										dGy**2								)
							!
							integral	=	prefact * num / denom
						!---------------------------------------------------------------------------------------------------------------
						!
						!
						!
						!CASE 3 (dGy==0)
						else if(	abs(dGy) < machineP	) then
							prefact		=	i_dp * dcmplx(-4.0_dp * at_rad * Gjy											) * myExp(dGx*x0)					
							num			=	dcmplx(		at_rad * dGx * cos(at_rad*dGx)		-		sin(at_rad*dGx)			)		
							denom		=	dcmplx(										dGx**2								)
							!
							integral	=	prefact * num / denom
						!---------------------------------------------------------------------------------------------------------------
						!
						!
						!
						!CASE 4 (dGx/=0 and dGy/=0)
						else if(	abs(dGx) >= machineP	 .and. 		abs(dGy) >= machineP	) then
							prefact		=	i_dp * dcmplx(4.0_dp 															) * myExp(dGx*x0+dGy*y0) 
			
							numA 		=	sin(at_rad*dGx) *		( 		at_rad * Gjx * dGx * dGy * cos(at_rad*dGy)		)
							numB		= -	sin(at_rad*dGy) *		(		at_rad * Gjy * dGx * dGy * cos(at_rad*dGx)	&
																		+	( Gjx*dGx - Gjy*dGy )	 * sin(at_rad*dGx)		)
							
							denom		=	dcmplx(		dGx**2 * dGy**2														)
							!
							!
							integral	=	prefact * ( numA + numB ) / denom	
						!---------------------------------------------------------------------------------------------------------------
						!
						!
						!DEFAULT
						else
							integral 	= 	dcmplx(0.0_dp)
							stop "[add_Zeeman]: reached forbidden default"
						end if
						!
						!add Ham
						ham = 	(alphaZee/dcmplx(vol)) * integral 
						Hmat(i,j)	= Hmat(i,j)		+ 	ham
						!
						!
						!GET MAX MIN CONTRIBUTION
						if( i==2 .and. j==1 ) then
							minI = abs(ham)
							maxI = minI
						end if
						!
						minI = min(	minI, abs(ham))
						maxI = max( maxI, abs(ham))
						!
						!
					end do
				end if
				!
				!
			end do
		end do

		!write(*,'(a,i3,a,e16.8,a,e16.8,a)')	"[#",myID,";add_Zeeman]: min/max contribution: [", minI," : ",maxI,"]."

		return
	end subroutine




end module ham_Zeeman