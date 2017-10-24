module polarization
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, nIntegrate
	use sysPara

	implicit none

	private
	public :: calcPolWannCent, calcPolViaA, calcIonicPol

	contains





!public:
	subroutine calcPolWannCent(wCent,pE)
		!calcuates 0 order polarization from the wannier centers
		!	this is only defined up to a uncertainty quantum of e*\vec{a}/V0, 
		!	where e is electron charge, \vec{a} a bravais lattice vector and V0 the volume of the unit cell
		!	uncertainty is resolved by projecting the wannier centers into the first unit cell at (0,0) 
		real(dp), intent(in)	:: wCent(:,:)
		real(dp), intent(out)	:: pE(2)
		integer 				:: n
		real(dp) 				:: cent(2), bondC(2)
		!
		
		pE	 	= 0.0_dp
		!
		!ATOM LIKE
		do n = 1,size(wCent,2)
			cent(1) = dmod(wCent(1,n),aX) !get current center by projection into first unit cell
			cent(2)	= dmod(wCent(2,n),aY)
			!
			
			!if( mod(n,2) == 0 ) then
			!	cent(:) = cent(:) - atPos(:,2)
			!else
			!	cent(:) = cent(:) - atPos(:,1)
			!end if

			write(*,'(a,f8.5,a,f8.5,a,f8.6,a,f8.6,a)')"[calc0ElPol]: Wcent = (",wCent(1,n),", ",wCent(2,n),") modified cent = (", cent(1),", ",cent(2),")"
			pE = pE + cent(:)				
		end do



	

		!2BAND
		!bondC(:)	= ( atPos(:,1) + atPos(:,2) ) / 2
		!do n = 1 , size(wCent,2)
		!	pE 	= pE + wCent(:,n) - bondC(:)
		!end do

		!MOD QUANTUM
		!pE(1) = dmod(pE(1),aX/vol)
		!pE(2) = dmod(pE(2),aY/vol)
		

		!
		!NORMALIZE
		!pE = pE / vol
		!pI = pI / vol
		!
		!SHIFT WITH RESPECT TO CENTER OF UNIT CELL
		!cent(1)	= aX * 0.5_dp
		!cent(2)	= aY * 0.5_dp
		!pE = (pE - cent ) / vol
		!

		return
	end subroutine



	subroutine calcPolViaA(A, pElA)
		!calculates the polarization by integrating connection over the brillouin zone
		! r_n 	= <0n|r|0n> 
		!		=V/(2pi)**2 \integrate_BZ <unk|i \nabla_k|unk>
		!		=V/(2pi)**2 \integrate_BZ A(k)
		complex(dp),		intent(in)		:: A(:,:,:,:)			!A(2,	 nWfs, nWfs, nQ	)	
		real(dp),			intent(out)		:: pElA(:)
		complex(dp)	,		allocatable		:: val(:)
		real(dp)							:: machine, bondC(2)
		integer								:: n, qi
		!
		allocate(	val( size(A,1) )	)
		val		= dcmplx(0.0_dp)
		machine	= 1e-15_dp
		
		!
		!
		!SUM OVER  STATES
		pElA	= 0.0_dp
		do n 	= 1, size(A,2)
			val	= dcmplx(0.0_dp)
			!INTEGRATE
			do qi = 1, size(A,4)
				val(1:2)	= val(1:2) + A(1:2,n,n,qi) / real(size(A,4),dp)
			end do
			write(*,'(a,i3,a,f8.4,a,f8.4,a)')	"[calcPolViaA]: n=",n,"dreal(p_n)=",dreal(val(1)),",",dreal(val(2)),")."

			if( abs(dimag(val(1))) > acc .or. abs(dimag(val(2))) > acc	) then
				write(*,*)	"[calcPolViaA]: found non zero imaginary contribution from band n=",n 
			end if
			!relative to bond center
			!if( mod(n,2) == 0 ) then
			!	val(:) = val(:) - atPos(:,2)
			!else
			!	val(:) = val(:) - atPos(:,1)
			!end if

			pelA(:)	= pElA(:) + dreal(val(:)) 
		end do
		

		!MOD QUANTUM
		!pelA(1)	= dmod(pElA(1),aX/vol)	
		!pelA(2)	= dmod(pElA(2),aY/vol)	

		!
		if( dimag(val(1)) > acc .or. dimag(val(2)) > acc ) then
			write(*,*)	"[calcPolViaA]: warning, connection has imaginary part none zeroDoGaugeTrafo(unkW, tHopp, EnH, AconnH, FcurvH, veloH)"
		end if 
		
		return
	end subroutine


	subroutine calcIonicPol(pI)
		!ionic contribution to total polarization
		real(dp),		intent(out)		:: pI(2)
		integer							:: at
		!
		pI	= 0.0_dp
		do at = 1, nAt
			pI	= pI + Zion(at) * atPos(:,at) 
		end do
		!
		!
		return
	end subroutine



end module polarization