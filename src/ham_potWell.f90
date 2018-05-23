module ham_PotWell
	!add Zeeman like term to hamitlonian
	use mpi
	use omp_lib
	use util_math,	only:	dp, PI_dp,i_dp, machineP, myExp
	use util_sysPara				
	implicit none	
	!#include "mpif.h"

	private
	public ::			add_potWell		



		

	


	contains


	!POTENTIAL WELLS
	subroutine add_potWell(nG_qi, Gvec, Hmat)
		!adds kinetic energy and potential Well term to Hamiltonian
		integer,		intent(in)		::	nG_qi
		real(dp),		intent(in)		::	Gvec(:,:)
		complex(dp),	intent(inout)	::	Hmat(:,:)
		integer							::	i, j
		!
		do j = 1, nG_qi
			do i = 1, nG_qi
				if( doVdesc )	then
					Hmat(i,j)	= Hmat(i,j)		+ Vdesc(Gvec, i,j)
				else
					Hmat(i,j)	= Hmat(i,j)		+ Vconst(Gvec, i,j)
				end if
			end do
		end do
		!
		return
	end subroutine




	complex(dp) function Vconst(Gvec, i,j)
		!calc potential matrix elements for constant potential inside each well
		!
		!the integrals were solved analytical and are hard coded in this function
		integer,		intent(in)	::	i, j
		real(dp),		intent(in)	::	Gvec(:,:)
		integer						::	at
		complex(dp)					::	Vpot
		real(dp)					::  xL, yL, xR, yR, dGx, dGy
		!
		Vconst 	= dcmplx(0.0_dp)
		dGx		= Gvec(1,j) - Gvec(1,i) 
		dGy		= Gvec(2,j) - Gvec(2,i) 
		!
		do at = 1, nAt
			Vpot	=	dcmplx(atPot(at))		!in hartree
			xL	=	atPos(1,at) - atR(1,at)
			xR	=	atPos(1,at) + atR(1,at) 
			yL	=	atPos(2,at) - atR(2,at)
			yR	=	atPos(2,at) + atR(2,at) 
			!
			!
			!
			!CASE 1 (dGx==0,dGy==0)
			if( i == j) then		
				Vconst	= Vconst + Vpot 			*	( xR - xL ) * 	( yR - yL )			 			/ vol
			
			!----------------------------------------------------------------------------------------------------------------------------------
			!
			!
			!CASE 2 (dGx==0)
			else if( abs(dGx) < machineP ) then	
				Vconst	= Vconst + Vpot  * i_dp 	* 	( xR - xL ) * ( myExp(dGy*yL) - myExp(dGy*yR) )	/ ( vol * dGy )
			
			!----------------------------------------------------------------------------------------------------------------------------------
			!
			!
			!CASE 3 (dGy==0)
			else if( abs(dGy) < machineP ) then
				Vconst	= Vconst + Vpot * i_dp	 	* 	( yR - yL) 	* ( myExp(dGx*xL) - myExp(dGx*xR) ) / (vol * dGx )
			
			!----------------------------------------------------------------------------------------------------------------------------------
			!
			!
			!CASE 4 (dGx/=0 ; dGy/=0)
			else if( abs(dGx) >= machineP 	.and. 		abs(dGy) >= machineP )	then
				Vconst	= Vconst -  Vpot 	 * ( myExp(dGx*xL) - myExp(dGx*xR) ) * ( myExp(dGy*yL) - myExp(dGy*yR) ) / (vol * dGx * dGy )
			
			!----------------------------------------------------------------------------------------------------------------------------------
			!
			!
			!DEFAULT
			else
				stop	"[Vconst]: reached forbidden region in dGx, dGy discussion"
			end if
		end do
		!
		return
	end function


	complex(dp) function Vdesc(Gvec, i,j)
		!potential integration for a well constant gradient in x direction (uniform electric field)
		!it starts from V0 at xL till V0-dV at xR
		integer,		intent(in)	::	i, j
		real(dp),		intent(in)	::	Gvec(:,:)
		integer						::	at
		complex(dp)					::	Vpot
		real(dp)					::  xL, yL, xR, yR, dGx, dGy, dV, fact
		!
		Vdesc 	= dcmplx(0.0_dp)
		dGx		= Gvec(1,j) - Gvec(1,i) 
		dGy		= Gvec(2,j) - Gvec(2,i) 
		!
		do at = 1, nAt
			Vpot=	dcmplx(atPot(at))
			dV	=	dVpot(at)
			xL	=	atPos(1,at) - atR(1,at)
			xR	=	atPos(1,at) + atR(1,at) 
			yL	=	atPos(2,at) - atR(2,at)
			yR	=	atPos(2,at) + atR(2,at) 
			!write(*,*)	"[",myId,"]dV=",dV
			!
			!
			!
			!
			!CASE 1 (dGx==0,dGy==0)
			if( i == j) then
				fact	= (2.0_dp*Vpot - dV)	 / (2.0_dP*vol)	
				Vdesc	= Vdesc 	+			fact *	(xL-xR) * 	(yL-yR)			
			
			!----------------------------------------------------------------------------------------------------------------------------------
			!
			!
			!CASE 2 (dGx==0)
			else if( abs(dGx) < machineP ) then	
				fact	= (2.0_dp*Vpot - dV) / ( 2.0_dP* vol * dGy )
				Vdesc	= Vdesc 	- 	i_dp  * fact * (xL-xR) * ( myExp(dGy*yL) - myExp(dGy*yR) )
			
			!----------------------------------------------------------------------------------------------------------------------------------
			!
			!
			!CASE 3 (dGy==0)
			else if( abs(dGy) < machineP ) then
				fact	=  1.0_dp / ( dGx**2 * vol * (xL-xR) )
				Vdesc	= Vdesc 	+ 	i_dp * fact * (yR-yL) * myExp(dGx*xL) * (dGx * Vpot			* (xL-xR) + i_dp * dV) 	
				Vdesc	= Vdesc 	-	i_dp * fact * (yR-yL) * myExp(dGx*xR) * (dGx * (Vpot-dV) * (xL-xR) + i_dp * dV)	
			
			!----------------------------------------------------------------------------------------------------------------------------------
			!
			!
			!CASE 4 (dGx/=0 ; dGy/=0)
			else if( abs(dGx) >= machineP 	.and. 		abs(dGy) >= machineP )	then
				fact	= -1.0_dp * ( myExp(dGy*yL)-myExp(dGy*yR) ) 	/ ( dGx**2 * dGy * vol * (xL-xR) )
				Vdesc	= Vdesc		+			fact * 			myExp(dGx*xL) * (dGx * Vpot			* (xL-xR) + i_dp * dV) 
				Vdesc	= Vdesc		-			fact *			myExp(dGx*xR) * (dGx * (Vpot-dV)	* (xL-xR) + i_dp * dV)
			
			!----------------------------------------------------------------------------------------------------------------------------------
			!
			!
			!DEFAULT
			else
				stop	"[Vdesc]: reached forbidden region in dGx, dGy discussion"
			end if
		end do
		!
		!
		return
	end function











end module ham_PotWell