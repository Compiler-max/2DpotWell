module blochWf
	!generates bloch and lattice periodidc functions given a basCoeff matrix
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, acc, machineP,& 

								myExp, myLeviCivita, nIntegrate
	use sysPara

	implicit none

	private
	public	::	calcBasis, UNKoverlap, calcVeloGrad


	contains







!public
	complex(dp) function UNKoverlap(n, m, qi, knb, gShift, ck)
		!calculates the overlap between unk at qi and at a neigbhouring k point knb
		!
		integer,		intent(in)		:: n, m, qi, knb
		real(dp),		intent(in)		:: gShift(2)
		complex(dp),	intent(in)		:: ck(:,:,:)  !ck(			nG		,	nBands  	,	nQ	)		
		integer							:: gi, gj, cnt
		real(dp)						:: delta(2)
		logical							:: notFound
		!
		UNKoverlap	= dcmplx(0.0_dp)
		cnt	= 0
		do gi = 1, nGq(qi)
			notFound 	= .true.
			gj			= 1
			do while( gj<= nGq(knb) .and. notFound ) 
				delta(1:2)	=  ( Gvec(1:2,gi,qi)-qpts(1:2,qi) ) 	-  		( Gvec(1:2,gj,knb)-qpts(1:2,knb)-gShift(1:2) )
				if( norm2(delta) < machineP )	then
					UNKoverlap	= UNKoverlap +  dconjg( ck(gi,n,qi) ) * ck(gj,m,knb) 
					cnt = cnt + 1
					notFound = .false.
				end if
				gj = gj + 1
			end do
			!if( gj>= nGq(knb) .and. notFound	) write(*,'(a,i3,a,i3)')	"[UNKoverlap]: no neighbour for gi=",gi," at qi=",qi
		end do
		!
		if( cnt > nGq(qi)	)	write(*,'(a,i8,a,i8)')	"[UNKoverlap]: warning, used ",cnt," where nGmax(qi)=",nGq(qi)
		if( cnt < nGq(qi) / 2.0_dp)	write(*,'(a,i8,a,i8)')	"[UNKoverlap]: warning, used  only",cnt," where nGmax(qi)=",nGq(qi)
		!
		return
	end function


	subroutine calcBasis(qi, ri, basVec)
		!calculates the basis vectors e^i(k+G).r
		!	if |k+G| is larger then the cutoff the basis vector is set to zero
		!	the cutoff enforces a symmetric base at each k point
		integer,	 intent(in)		:: qi, ri
		complex(dp), intent(out)	:: basVec(:)
		integer 				 	:: i 
		!
		basVec	= 0.0_dp
		do i =1, nGq(qi)
			basVec(i) 		= myExp( dot_product( Gvec(1:2,i,qi), rpts(1:2,ri) )		)  !/ dsqrt(vol)
		end do
		!
		!
		return
	end subroutine


	subroutine calcVeloGrad(ck, v_mat )
		!calculates the velocity operator matrix
		!	Psi_n v Psi_m	= i/hbar Psi_n grad_r Psi_m
		!					= - 1 / hbar sum_G ckn^dag ckm G
		complex(dp),	intent(in)		:: 	ck(:,:,:)
		complex(dp),	intent(out)		::	v_mat(:,:,:,:)
		integer							::	qi, m, n, gi
		!
		v_mat = dcmplx(0.0_dp)
		!
		!
		if(	size(ck,3)/=size(v_mat,4)	) then
			write(*,*)	"[calcVeloGrad]: coeff and velo defined on different k meshes, stop now"
		else
			do qi = 1, size(ck,3)
				do m = 1, nWfs
					do n = 1, nWfs
						!SUM OVER BASIS FUNCTIONS
						do gi = 1 , nGq(qi)
							v_mat(1:2,n,m,qi) = v_mat(1:2,n,m,qi) -  dconjg(ck(gi,n,qi)) *  ck(gi,m,qi) *  Gvec(1:2,gi,qi)
						end do
					end do
				end do
			end do
		end if
		!
		return
	end subroutine




end module blochWf 


