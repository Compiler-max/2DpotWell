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
		!
		call Vconst(nG_qi, Gvec, Hmat)
		!
		if( doVdesc ) then
			call Vdesc(nG_qi, Gvec, Hmat)
		end if

		!old: disadvantages: function call overhead, the const potential contribution was calculated in different way
		!do j = 1, nG_qi
		!	do i = 1, nG_qi
		!		if( doVdesc )	then
		!			Hmat(i,j)	= Hmat(i,j)		+ Vdesc(Gvec, i,j)
		!		else
		!			Hmat(i,j)	= Hmat(i,j)		+ Vconst(Gvec, i,j)
		!		end if
		!	end do
		!end do
		!
		return
	end subroutine




	subroutine Vconst(nG_qi, Gvec, Hmat)
		!calc potential matrix elements for constant potential inside each well
		!
		!the integrals were solved analytical and are hard coded in this function
		integer,		intent(in)		::	nG_qi
		real(dp),		intent(in)		::	Gvec(:,:)
		complex(dp),	intent(inout)	::	Hmat(:,:)
		integer							::	at, i, j
		complex(dp)						::  Vpot, xL, yL, xR, yR, dGx, dGy, cvol
		!
		cvol	= dcmplx(	vol	)
		!
		do at = 1, nAt
			Vpot=	dcmplx(		atPot(at)					)		!in hartree
			xL	=	dcmplx(		atPos(1,at) - atR(1,at)		)  
			xR	=	dcmplx(		atPos(1,at) + atR(1,at)		) 
			yL	=	dcmplx(		atPos(2,at) - atR(2,at)		)
			yR	=	dcmplx(		atPos(2,at) + atR(2,at)		) 
			!
			!
			do j = 1, nG_qi
				do i = 1, nG_qi
					dGx		= dcmplx(	Gvec(1,j) - Gvec(1,i) )
					dGy		= dcmplx(	Gvec(2,j) - Gvec(2,i) )
					!
					!CASE 1 (dGx==0,dGy==0)
					if( i == j) then		
						Hmat(i,j)	= Hmat(i,j) + Vpot 			*	( xR - xL ) * 	( yR - yL )			 			/ cvol
					
					!----------------------------------------------------------------------------------------------------------------------------------
					!
					!
					!CASE 2 (dGx==0)
					else if( abs(dGx) < machineP ) then	
						Hmat(i,j)	= Hmat(i,j) + Vpot  * i_dp 	* 	( xR - xL ) * ( myExp(dreal(dGy*yL)) - myExp(dreal(dGy*yR)) )	/ ( cvol * dGy )
					
					!----------------------------------------------------------------------------------------------------------------------------------
					!
					!
					!CASE 3 (dGy==0)
					else if( abs(dGy) < machineP ) then
						Hmat(i,j)	= Hmat(i,j) + Vpot * i_dp	 	* 	( yR - yL) 	* ( myExp(dreal(dGx*xL)) - myExp(dreal(dGx*xR)) ) / (cvol * dGx )
					
					!----------------------------------------------------------------------------------------------------------------------------------
					!
					!
					!CASE 4 (dGx/=0 ; dGy/=0)
					else if( abs(dGx) >= machineP 	.and. 		abs(dGy) >= machineP )	then
						Hmat(i,j)	= Hmat(i,j) -  Vpot 	 * ( myExp(dreal(dGx*xL)) - myExp(dreal(dGx*xR)) ) * ( myExp(dreal(dGy*yL)) - myExp(dreal(dGy*yR)) ) / (cvol * dGx * dGy )
					
					!----------------------------------------------------------------------------------------------------------------------------------
					!
					!
					!DEFAULT
					else
						stop	"[Vconst]: reached forbidden region in dGx, dGy discussion"
					end if
					!
					!
				end do
			end do
			!
			!
		end do
		!
		return
	end subroutine




	subroutine Vdesc(nG_qi, Gvec, Hmat)
		!potential integration for a well constant gradient in x direction (uniform electric field)
		!it starts from V0 at xL till V0-dV at xR
		integer,		intent(in)		::	nG_qi
		real(dp),		intent(in)		::	Gvec(:,:)
		complex(dp),	intent(inout)	::	Hmat(:,:)
		integer							::	at, i, j
		complex(dp)						::  xL, yL, xR, yR, dGx, dGy, dV, fact, cvol, int1, int2
		!
		cvol = dcmplx( vol )
		!
		do at = 1, nAt
			dV	=	dcmplx(		dVpot(at)					)	
			xL	=	dcmplx(		atPos(1,at) - atR(1,at)		)
			xR	=	dcmplx(		atPos(1,at) + atR(1,at) 	)
			yL	=	dcmplx(		atPos(2,at) - atR(2,at) 	)
			yR	=	dcmplx(		atPos(2,at) + atR(2,at) 	)
			!write(*,*)	"[",myId,"]dV=",dV
			!
			!
			do j = 1, nG_qi
				do i = 1, nG_qi
					dGx		= dcmplx(	Gvec(1,j) - Gvec(1,i) )
					dGy		= dcmplx(	Gvec(2,j) - Gvec(2,i) )
					!
					!
					!CASE 1 (dGx==0,dGy==0)
					if( i == j) then
						fact	= dV	 / (2.0_dP*cvol)
						int1	= -1.0_dp*(xL-xR) * 	(yL-yR)	
						!	
						Hmat(i,j)	= Hmat(i,j) 	+	fact *	 int1		
					
					!----------------------------------------------------------------------------------------------------------------------------------
					!
					!
					!CASE 2 (dGx==0)
					else if( abs(dGx) < machineP ) then	
						fact	= dV / ( 2.0_dP* cvol * dGy )
						int1	= (xL-xR) * i_dp * ( myExp(dreal(dGy*yL)) - myExp(dreal(dGy*yR)) )
						!
						Hmat(i,j)	= Hmat(i,j) 	+ 	fact * int1
					
					!----------------------------------------------------------------------------------------------------------------------------------
					!
					!
					!CASE 3 (dGy==0)
					else if( abs(dGy) < machineP ) then
						fact	=  dV / ( dGx**2 * cvol * (xL-xR) )
						int1	= (yR-yL)  * i_dp *   ( 	myExp(dreal(dGx*xL)) 	+ 		i_dp * myExp(dreal(dGx*xR)) * (i_dp+dGx*(xR-xL))	)
						!
						Hmat(i,j)	= Hmat(i,j) 	+ 	fact * int1
									
					!----------------------------------------------------------------------------------------------------------------------------------
					!
					!
					!CASE 4 (dGx/=0 ; dGy/=0)
					else if( abs(dGx) >= machineP 	.and. 		abs(dGy) >= machineP )	then
						fact	=	dV / (dGx**2 * dGy * cvol * (xL-xR) )
						int1	=	myExp( dreal(dGy*yL)	) - myExp( dreal(dGy*yR)	)
						int2	=	-i_dp * myExp( dreal(dGx*xL)	) + myExp( dreal(dGx*xR)	) * (i_dp+ dGx*(xR-xL))				
						!
						Hmat(i,j)	= Hmat(i,j)		+		fact * int1 * int2 	

					
					!----------------------------------------------------------------------------------------------------------------------------------
					!
					!
					!DEFAULT
					else
						stop	"[Vdesc]: reached forbidden region in dGx, dGy discussion"
					end if
					!
				end do
			end do
			!
			!
		end do
		!
		!
		return
	end subroutine











end module ham_PotWell