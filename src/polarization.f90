module polarization
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp
	use sysPara,		only:	nAt, atPos, Zion

	implicit none

	private
	public :: calcPolWannCent, calcPolViaA

	contains





!public:
	subroutine calcPolWannCent(wCent,pE)
		!calcuates 0 order polarization from the wannier centers
		!	this is only defined up to a uncertainty quantum of e*\vec{a}/V0, 
		!	where e is electron charge, \vec{a} a bravais lattice vector and V0 the volume of the unit cell
		!	uncertainty is resolved by projecting the wannier centers into the first unit cell at (0,0) 
		real(dp), intent(in)	:: wCent(:,:)
		real(dp), intent(out)	:: pE(3)
		integer 				:: n
		real(dp) 				:: cent(3)
		!
		pE	 	= 0.0_dp
		!
		do n = 1,size(wCent,2)
			!grep center
			cent(1:3)	= wCent(1:3,n)
			!substract atom center from pol
			call substractAtPos(n,cent)
			!debug message
			write(*,'(a,f12.5,a,f12.5,a,f12.5,a,f12.5,a)')"[calc0ElPol]: Wcent = (",&
								wCent(1,n),", ",wCent(2,n),") modified cent = (", cent(1),", ",cent(2),")"
			!sum total
			pE = pE + cent(:)				
		end do
		!
		!
		return
	end subroutine



	subroutine calcPolViaA(A_mat, centers)
		!calculates the polarization by integrating connection over the brillouin zone
		! r_n 	= <0n|r|0n> 
		!		=V/(2pi)**2 \integrate_BZ <unk|i \nabla_k|unk>
		!		=V/(2pi)**2 \integrate_BZ A(k)
		complex(dp),		intent(in)		:: A_mat(:,:,:,:)			!A(2,	 nWfs, nWfs, nQ	)	
		real(dp),			intent(out)		:: centers(:,:)
		complex(dp)	,		allocatable		:: val(:)
		integer								:: n
		!
		allocate(	val(  size(A_mat,1) )	)
		!
		do n 	= 1, size(A_mat,2)
			!
			!INTEGRATE
			val(1) = sum(A_mat(1,n,n,:)) / size(A_mat,4)
			val(2) = sum(A_mat(2,n,n,:)) / size(A_mat,4)
			if(size(A_mat,1)==3)	val(3) = sum(A_mat(3,n,n,:)) / size(A_mat,4)
			!
			!COLLECT REAL PART
			centers(:,n) 	= dreal(val(:))
			!
			!DEBUG MESSAGE
			write(*,'(a,i3,a,f8.4,a,f8.4,a)')	"[calcPolViaA]: n=",n,"p_n=",dreal(val(1)),",",dreal(val(2)),")."
			if( abs(dimag(val(1))) > acc .or. abs(dimag(val(2))) > acc .or. abs(dimag(val(3))) > acc	) then
				write(*,*)	"[calcPolViaA]: found non zero imaginary contribution from band n=",n 
			end if	
			!
		end do
		!
		!		
		return
	end subroutine


!private
	subroutine substractAtPos(n, cent)
		integer,		intent(in)		:: n
		real(dp),		intent(inout)	:: cent(:)
		!
		!!SINGLE ATOM
		if( nAt == 1 ) then
			cent(1:2)	= cent(1:2) - atPos(1:2,1)		!calc center w.r.t. atom center
		!!DOUBLE ATOM
		else if( nAt == 2 ) then
			if( mod(n,2)== 0 ) then
				cent(1:2)	= cent(1:2) - atPos(1:2,2)
			else	
				cent(1:2)	= cent(1:2) - atPos(1:2,1)
			end if
		!!DEFAULT
		else
			write(*,*)	"[calcPolWannCent]: to many atoms in unit cell (more then 2)!"
		end if
		!
		return
	end subroutine

end module polarization