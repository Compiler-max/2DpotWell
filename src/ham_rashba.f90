module ham_Rashba
	!add Zeeman like term to hamitlonian
	use util_math,		only:	dp, PI_dp, machineP, myExp
	use util_sysPara				
	implicit none	
	!#include "mpif.h"

	private
	public ::			add_rashba			



		

	


	contains



	subroutine add_rashba(nG_qi, Gvec, Hmat)
		!
		!	adds term Gi_y to the Hamiltonian
		!
		integer,		intent(in)		::	nG_qi
		real(dp),		intent(in)		::	Gvec(:,:)
		complex(dp),	intent(inout)	::	Hmat(:,:)
		complex(dp)						::	rashHam
		integer							::	gi
		!
		!write(*,'(a,i3,a,e14.4,a)')	"[#",myID,";add_rashba]: hello there, aRashba=",aRashba," a.u."
		!if( qLoc == 1 ) write(*,*)	"[#",myID,";add_rashba]: use_px_rashba = ",use_px_rashba
		!
		do gi = 1, nG_qi	
			!
			if(	use_px_rashba )	then
				rashHam	=  	dcmplx(	aRashba * Gvec(1,gi) 	)
			else
				rashHam =	dcmplx( aRashba * Gvec(2,gi)	)
			end if
			!
			Hmat(gi,gi) = Hmat(gi,gi) + rashHam 
		end do
		!
		!
		return
	end subroutine




end module ham_Rashba