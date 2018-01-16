module peierls
	use mathematics,	only:	dp, PI_dp, i_dp, myExp, crossP, eigSolverFULL, isUnit, aUtoAngstrm, aUtoEv
	use sysPara
	use blochWf,		only:	calcConnOnCoarse
	use	polarization,	only:	calcPolViaA
	use output,			only:	writePeierls, writeHtb
	implicit none
	

	private
	public ::	peierlsMethod


	integer							::	num_wann, nrpts, R_null

	integer,		allocatable		:: 	wigStzDegn(:), R_vect(:,:)
	real(dp)						::	recip_latt(3,3)
	real(dp),		allocatable		:: 	R_real(:,:), wCent(:,:)
	complex(dp),	allocatable		::	H_tb(:,:,:), r_tb(:,:,:,:)


	contains







!public:
	subroutine	peierlsMethod(ck,  pPei)
		complex(dp),	intent(in)		::	ck(:,:,:)	! tHopp(nWfs,nWfs,nSC)
		real(dp),		intent(out)		::	pPei(3)
		logical							::	foundFile
		!
		!GET tHopp
		call readTBsingle( foundFile )
		!
		!PEIERLS SUB
		if( foundFile ) then
			call doPeierls(ck, pPei)
		else
			write(*,*)	"[peierlsMethod]: did not find input file seedname_tb.dat"
			pPei	= 0.0_dp
		end if
		!
		!
		write(*,*)	"[peierlsMethod]: calculated polarization, by.."
		!
		!
		return
	end subroutine










!privat:
	real(dp) function shift(rj)
		!integrates the vector potential A(r) analytically   (assumes uniform field)
		!	A(r)	= -0.5 cross_p[r,B]
		!	return 	= Integrate^\vec{rMax]_\vec{rMin} 		\vec{A(r')}.\vec{dr'}
		!			= -0.5 cross_p[rMin,B].(rMax - rMin)
		integer,		intent(in)		:: rj
		real(dp)						:: rMin(2), rMax(2), rU(3), rL(3), integrateA
		rMin(:)		= 0.0_dp
		rMax(:)		= Rcell(:,rj)
		!
		rU(1:2)		= rMax(1:2)
		rU(3)		= 0.0_dp
		rL(1:2)		= rMin(1:2)
		rL(3)		= 0.0_dp
		!
		!old (symmetric gauge - A_vec =  -0.5 * crossp(R_vec,B_vec)		)
		!integrateA	= -0.5_dp * dot_product( crossP(rL,Bext)	, (rU-rL)	)
		!
		!new (Landau Gauge - A_vec = y_hat * abs(B) R_vec(1)	)
		integrateA  = 0.5_dp * Bext(3) * Rcell(1,rj) * Rcell(2,rj) 
		!
		!
		shift		= myExp( integrateA )
		!
		return
	end function


	subroutine doPeierls(ck, pPei)
		!performs peierls sub & calcs pol
		complex(dp),	intent(in)		::	ck(:,:,:)
		real(dp),		intent(out)		::	pPei(3)
		complex(dp),	allocatable		::	Hp(:,:), ckP(:,:,:), tshift(:,:,:), AconnP(:,:,:,:)
		real(dp),		allocatable		::	EnP(:,:)
		integer							::	Ri, qi, gi
		complex(dp)						::	phase
		!
		
		allocate(			ckP(		nG		,	nWfs	,	nQ		)			)
		allocate(			tshift(		nWfs	, 	nWfs	,	nSc		)			)
		allocate(			EnP(					nWfs	,	nQ		)			)
		allocate(			AconnP(3,	nWfs	,	nWfs	,	nQ		)			)
		!
		!
		!DO PEIERLS SUBSTITUTION
		do Ri = 1, nSC
			tshift(:,:,Ri)	= H_tb(:,:,Ri) * shift(Ri)
		end do
		write(*,*)	"[peierlsMethod]: substiution of hopping parameters done"
		!
		!
		!ELECTRONIC STRUCTURE

		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(qi, gi, Ri, Hp, phase)
		allocate(	Hp(	nWfs,	nWfs)		)
		!$OMP DO SCHEDULE( STATIC )
		do qi = 1, nQ
			!SET UP k SPACE HAMILTONIAN
			Hp	= dcmplx(0.0_dp)
			do Ri = 1, nSC
				phase	= myExp( dot_product(qpts(:,qi),Rcell(:,Ri))	) !/ dsqrt(real(nSC,dp))
				Hp(:,:)	= Hp(:,:) + phase * tshift(:,:,Ri)
			end do
			!SOLVE HAM	
			call eigSolverFULL(Hp(:,:),EnP(:,qi))
			!
			if( .not. isUnit(Hp)	) write(*,*) "[peierlsMethod]: ckP not unitary at qi=",qi
			!EXTRACT EXPANSION COEFF
			ckP(:,:,qi)	= dcmplx(0.0_dp)
			do gi = 1, nGq(qi)
				ckP(gi,:,qi)	= matmul( ck(gi,:,qi), Hp(:,:))
			end do
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		!
		!
		!GENERATE CONNECTION
		AconnP	= dcmplx(0.0_dp)
		call calcConnOnCoarse(ckP, AconnP) 
		!
		!CALC POL
		call calcPolViaA(AconnP, pPei)
		!
		!WRITE UNKs & ENERGIES
		if( writeBin ) call writePeierls(ckP, EnP)
		!
		!
		return
	end subroutine









	subroutine readTBsingle( readSuccess )
		!reads the _tb.dat file given by wannier90
		logical,		intent(out)		::	readSuccess
		integer							:: 	stat, cnt, offset, R, n, m, i, mn(2), dumI(3), line15(15), inq_unit
		real(dp)						::	real2(2), real6(6), real3(3)
		character(len=10)				::	fname
		logical							::	inq_open, inq_exist
		!
		fname		= seedName//'_tb.dat'
		


		!try opening file
		open(unit=310, iostat=stat, file=fname, status='old', action='read' )
		inquire(file=fname,opened=inq_open, exist=inq_exist, number=inq_unit)
		if( .not. inq_exist	) write(*,*)	"[Pei/readTBsingle]: file ",fname, " does not exist (unit=",inq_unit,")."
		if( .not. inq_open	) write(*,*)	"[Pei/readTBsingle]: file ",fname, " is not open(unit=",inq_unit,")."
		


		if( stat /= 0)  then
			write(*,*) "[Pei/readTBsingle]: warning, file seedname_tb.dat not found"
			readSuccess 	= .false.
			recip_latt		= 0.0_dp
			R_real			= 0.0_dp
			H_tb			= dcmplx(0.0_dp)
			r_tb			= dcmplx(0.0_dp)
		else
			readSuccess	= .true.
			!
			read(310,*)
			!recip lattice (read into buffer, avoids compiler warning)
			read(310,*) 		real3(:)
			recip_latt(1,:)	= 	real3(:)
			read(310,*)			real3(:)
			recip_latt(2,:)	= 	real3(:)
			read(310,*)			real3(:)
			recip_latt(3,:)	= 	real3(:)
			!sys info
			read(310,*) num_wann
			read(310,*)	nrpts
			!
			allocate( 	wigStzDegn(									nrpts	)		)
			allocate(	R_vect(			3,							nrpts	)		)
			allocate(	R_real(			3,							nrpts	)		)
			allocate(	H_tb(				num_wann,	num_wann, 	nrpts	)		)
			allocate(	r_tb(			3,	num_wann,	num_wann,	nrpts	)		)	
			allocate(	wCent(			3,		num_wann					)		)	
			!
			!read degeneracy of each wigner seitz grid point
			cnt 	= 0
			offset 	= 0
			if(	nrpts <= 15 ) then
				read(310,*) wigStzDegn
			else
				do while ( cnt < nrpts )		!read 15 entries per line till nrpts real2ues are read 
					if( nrpts - cnt >= 15	) then
						read(310,*)		line15
						cnt	= cnt + 15
						do i = 1 , 15
							if(offset+i <= size(wigStzDegn))	wigStzDegn(offset+i)	= line15(i) 
						end do
						offset=	offset + 15
					else						!the last line might contain less then 15 entries
						read(310,*)	wigStzDegn( (offset+1):(offset+(nrpts-cnt)) )
						cnt = cnt + (nrpts-cnt)
						offset = offset + (nrpts-cnt)
					end if
				end do
			end if
			!
			!READ HOPPINGS
			do R = 1, nrpts
				!skip first line
				read(310,*)
				!second line fractional real6 of R
				read(310,*)	dumI(1:3)
				R_vect(1:3,R)	= dumI(1:3)
				if( R_vect(1,R)==0 .and. R_vect(2,R)==0 .and. R_vect(3,R)==0 ) R_null = R
				!get indices & hoppings
				do n = 1, num_wann
					do m = 1, num_wann
						read(310,*) mn(1:2), real2(1:2)
						H_tb(mn(1),mn(2),R)	= dcmplx(real2(1)) + i_dp * dcmplx(real2(2))
					end do
				end do
			end do
			!
			!READ POSITIONS
			do R = 1, nrpts
				!skip first line
				read(310,*)
				!second line fractional real6 of R
				read(310,*)	dumI(1:3)
				if( dumI(1) /= R_vect(1,R) .or. dumI(2) /= R_vect(2,R) .or. dumI(3) /= R_vect(3,R) ) then
					write(*,*)	"[Pei/readTBsingle]: warning, detected R vector ordering issue while reading positions"
				end if
				!
				do n = 1, num_wann
					do m = 1, num_wann
						read(310,*) mn(1:2), real6(1:6)
						r_tb(1,mn(1),mn(2),R)	= dcmplx(real6(1)) + i_dp * dcmplx(real6(2))
						r_tb(2,mn(1),mn(2),R)	= dcmplx(real6(3)) + i_dp * dcmplx(real6(4))
						r_tb(3,mn(1),mn(2),R)	= dcmplx(real6(5)) + i_dp * dcmplx(real6(6))
					end do
				end do
			end do
			close(310)
			!
			!CONVERT BACK TO [a.u.]
			H_tb	= H_tb /	aUtoEv 
			r_tb	= r_tb /	aUtoAngstrm
			!
			!get centers
			do n = 1, num_wann
				wCent(1:3,n)	= r_tb(1:3,n,n,R_null)
			end do
			!
			!calculate real R vector
			do R = 1 , nrpts
				R_real(1,R)	= R_vect(1,R)	* aX
				R_real(2,R)	= R_vect(2,R)	* aY
				R_real(3,R)	= R_vect(3,R)	* 0
			end do
			!
			call writeHtb(H_tb)
			!
		end if
		!
		!
		return
	end subroutine




end module peierls