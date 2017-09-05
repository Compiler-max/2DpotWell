module polarization
use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, nIntegrate
	use sysPara

	implicit none

	private
	public :: calcPolWannCent, calcPolViaA

	contains





!public:
	subroutine calcPolWannCent(wCent,pE, pI, pT)
		!calcuates 0 order polarization from the wannier centers
		!	this is only defined up to a uncertainty quantum of e*\vec{a}/V0, 
		!	where e is electron charge, \vec{a} a bravais lattice vector and V0 the volume of the unit cell
		!	uncertainty is resolved by projecting the wannier centers into the first unit cell at (0,0) 
		real(dp), intent(in)	:: wCent(:,:)
		real(dp), intent(out)	:: pE(2), pI(2), pT(2)
		integer 				:: n, at
		real(dp) 				:: vol, cent(2)
		!
		
		vol 	= aX * aY !1D
		pT		= 0.0_dp
		pE	 	= 0.0_dp
		pI		= 0.0_dp
		!
		!ELECTRONIC
		do n = 1,size(wCent,2)
			cent(1) = dmod(wCent(1,n),aX) !get current center by projection into first unit cell
			cent(2)	= dmod(wCent(2,n),aY)
			!!
			!write(*,'(a,f8.5,a,f8.5,a,f8.6,a,f8.6,a)')"[calc0ElPol]: Wcent = (",wCent(1,n),", ",wCent(2,n),") modified cent = (", cent(1),", ",cent(2),")"
			pE = pE + cent				
		end do
		!IONIC
		do at = 1, nAt
			pI = pI + Zion(at) * atPos(:,at) 
		end do
		!NORMALIZE
		!pE = pE / vol
		!pI = pI / vol
		!
		!SHIFT WITH RESPECT TO CENTER OF UNIT CELL
		cent(1)	= aX * 0.5_dp
		cent(2)	= aY * 0.5_dp
		pE = (pE - cent ) / vol
		pI = (pI - cent ) / vol 
		!TOTAL
		pT = pI + pE												
		!
		return
	end



	subroutine calcPolViaA(A, pElA)
		!calculates the polarization by integrating connection over the brillouin zone
		! r_n 	= <0n|r|0n> 
		!		=V/(2pi)**2 \integrate_BZ <unk|i \nabla_k|unk>
		!		=V/(2pi)**2 \integrate_BZ A(k)
		real(dp),		intent(in)		:: A(:,:,:,:)			!A(2,	 nK, nWfs, nWfs	)	
		real(dp),		intent(out)		:: pElA(:)
		complex(dp)	,	allocatable		:: val(:)
		real(dp)						:: machine
		integer							:: n, ki
		!
		allocate(	val( size(A,1) )	)
		val		= dcmplx(0.0_dp)
		pElA	= 0.0_dp
		machine	= 1e-15_dp
		!
		!SUM OVER K SPACE AND OVER STATES
		do n 	= 1, size(A,3)
			do ki = 1, size(A,2)
				val(:)	= val(:) + A(:,ki,n,n)
			end do
		end do
		!
		!NORMALIZE
		val		= val / real(size(A,2),dp)
		!
		!HARVEST
		pElA	= val !/ (aX*aY)
		!
		!MOD TO FIRST UNIT CELL
		pElA(1)	= mod( pElA(1), aX )
		pElA(2) = mod( pElA(2), aY )
		!
		!DEBUGGING
		if(		dimag( val(1) ) > machine 	.or. 	dimag( val(2) ) > machine		) then
			write(*,*)"[calcPolViaA]: non zero imaginary part of polarization"
		end if
		return
	end



end module polarization