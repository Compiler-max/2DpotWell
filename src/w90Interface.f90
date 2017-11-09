module w90Interface	

	use mathematics,	only:	dp, PI_dp, i_dp, machineP, aUtoAngstrm, aUtoEv, myExp
	use sysPara

	implicit none

	private
	public ::					w90Interf

	character(len=3)					:: 	seed_name
	character(len=20),	allocatable		:: 	atom_symbols(:)
	logical								:: 	gamma_only, spinors
	logical,			allocatable		::	lwindow(:,:)
	integer								:: 	mp_grid(3), num_kpts, num_bands_tot, num_atoms, num_nnmax, &
											nntot, num_bands, num_wann
	integer,			allocatable		:: 	nnlist(:,:), nncell(:,:,:), proj_l(:), proj_m(:), &
											proj_radial(:), exclude_bands(:), proj_s(:)
	real(dp)							:: 	real_lattice(3,3), recip_lattice(3,3)
	real(dp),			allocatable		:: 	kpt_latt(:,:), atoms_cart(:,:), eigenvalues(:,:), &
											proj_site(:,:), proj_z(:,:), proj_x(:,:), proj_zona(:), proj_s_qaxis(:,:),&
											wann_centres(:,:), wann_spreads(:), spread(:)
	complex(dp),		allocatable		::	M_matrix(:,:,:,:), M_matrix_orig(:,:,:,:), &	
											A_matrix(:,:,:), U_matrix(:,:,:), U_matrix_opt(:,:,:)									


	contains




!public:
	subroutine w90Interf(ck, En)
		complex(dp),	intent(in)		:: ck(:,:,:)	
		real(dp),		intent(in)		:: En(:,:)

		call writeW90input()

		write(*,*)	"[w90Interf]: wrote w90 input file"
		call w90setup()
		write(*,*)	"[w90Interf]: done with w90 setup"
		write(*,*)	"[w90Interf]: will use num_bands= ",num_bands, " bands"
		write(*,*)	"[w90Interf]: to gen   num_wann=  ",num_wann, " wnfs"

		allocate(	eigenvalues(							num_bands		,	num_kpts	)		)
		allocate(	M_matrix(		num_bands	,	num_bands	,	nntot	,	num_kpts	)		)
		allocate(	M_matrix_orig(	num_bands	,	num_bands	,	nntot	,	num_kpts	)		)
		allocate(	A_matrix(		num_bands	,			num_wann		,	num_kpts	)		)
		allocate(	U_matrix(		num_wann	,			num_wann		,	num_kpts	)		)
		allocate(	U_matrix_opt(	num_wann	,			num_wann		,	num_kpts	)		)
		allocate(	lwindow(		num_bands								,	num_kpts	)		)
		allocate(	wann_centres(	3			,			num_wann						)		)
		allocate(	wann_spreads(							num_wann						)		)
		allocate(	spread(									num_wann						)		)

		call w90prepEigVal(En)
		call w90prepMmat(ck)
		call w90prepAmat(ck)
		write(*,*)	"[w90Interf]: done preparing wannierization input matrices"

		call wannier_run(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice, &
							kpt_latt,num_bands,num_wann,nntot,num_atoms,atom_symbols, &
							atoms_cart,gamma_only,M_matrix_orig,A_matrix,eigenvalues, &
							U_matrix,U_matrix_opt,lwindow,wann_centres,wann_spreads, &
							spread)


	end subroutine







!prviate
	subroutine writeW90input()
		integer				:: stat
		integer				:: nw90it
		!
		nw90it = 100
		seed_name			= 'wf1'
		!
		!Delete old input file
		open(unit=200, iostat=stat, file=seed_name//'.win', status='old')
		if (stat == 0) close(200, status='delete')
		!Delete old wout file
		open(unit=210, iostat=stat, file=seed_name//'.wout', status='old')
		if (stat == 0) close(210, status='delete')
		!		
		!
		!
		!create new
		open(unit=100,file=seed_name//'.win',action='write', status='new')
		!
		!BASIC INFO
		write(100,*)	'num_wann  = ',nWfs
		write(100,*)	'num_bands = ',nBands
		write(100,*)	'mp_grid   = ', mp_grid(1) , ' ', mp_grid(2), ' ', mp_grid(3)
		write(100,*)	'num_iter  = ', nW90it
		write(100,*)	''
		!
		!PROJECTIONS
		write(100,*)	'Begin Projections'
		write(100,*)	'H: l=0;l=1,mr=1;l=1,mr=2'
		write(100,*)	'End Projections'

		!
		!write(100,*)	'begin unit_cell_cart'
		!write(100,*)	'bohr'
		!write(100,*)	real_lattice(1,1), ' ', real_lattice(1,2), ' ', real_lattice(1,3)	
		!write(100,*)	real_lattice(2,1), ' ', real_lattice(2,2), ' ', real_lattice(2,3)	
		!write(100,*)	real_lattice(3,1), ' ', real_lattice(3,2), ' ', real_lattice(3,3)		
		!write(100,*)	'end unit_cell_cart'
		!
		!JOBS
		write(100,*)	'write_tb = true'

		close(100)
		return
	end subroutine




	subroutine w90setup()
		integer							:: at, n
		!
		!general info
		gamma_only			= .false.
		if( nQ==1 )			gamma_only	= .true.
		spinors				= .false.
		!
		!q point grid
		mp_grid(1)			= nQx
		mp_grid(2)			= nQy
		mp_grid(3)			= 1
		num_kpts			= nQ
		num_nnmax			= 12
		!unit cell
		real_lattice		= 0.0_dp
		real_lattice(1,1)	= aX * aUtoAngstrm 
		real_lattice(2,2)	= aY * aUtoAngstrm
		real_lattice(3,3)	= 10.0_dp * (aX * aY) * aUtoAngstrm
		!reciprocal cell
		recip_lattice		= 0.0_dp
		recip_lattice(1,1)	= 2.0_dp * PI_dp / ( aX 				* aUtoAngstrm )		
		recip_lattice(2,2)	= 2.0_dp * PI_dp / ( aY 				* aUtoAngstrm )  
		recip_lattice(3,3)	= 2.0_dp * PI_dp / (real_lattice(3,3) 	* aUtoAngstrm )
		!atoms
		num_bands_tot		= nBands
		num_atoms			= nAt
		!		
		!
		allocate(	atom_symbols(					num_atoms						)			)
		allocate(	kpt_latt(			3	,		num_kpts						)			)
		allocate(	atoms_cart(			3	,		num_atoms						)			)
		allocate(	nnlist(							num_kpts	, 	num_nnmax		)			)
		allocate(	nncell(				3	,		num_kpts	,	num_nnmax		)			)
		allocate(	proj_site(			3	,		num_bands_tot					)			)
		allocate(	proj_l(							num_bands_tot					)			)
		allocate(	proj_m(							num_bands_tot					)			)
		allocate(	proj_radial(					num_bands_tot					)			)
		allocate(	proj_z(				3	,		num_bands_tot					)			)
		allocate(	proj_x(				3	,		num_bands_tot					)			)
		allocate(	proj_zona(						num_bands_tot					)			)
		allocate(	exclude_bands(					num_bands_tot					)			)
		allocate(	proj_s(							num_bands_tot					)			)
		allocate(	proj_s_qaxis(		3	,		num_bands_tot					)			)
		!
		!fill atom related arrays
		do	at = 1, num_atoms
			atom_symbols(at)	= 'H'
			atoms_cart(1:2,at)	= atPos(1:2,at)
			atoms_cart(3,at)	= 0.0_dp
		end do
		!kpt lattice (fractional co-ordinates relative to the reciprocal lattice vectors)
		kpt_latt(1,:)			= qpts(1,:) / ( recip_lattice(1,1) * aUtoAngstrm ) 
		kpt_latt(2,:)			= qpts(2,:)	/ ( recip_lattice(2,2) * aUtoAngstrm ) 
		kpt_latt(3,:)			= 0.0_dp
		!
		!			
		call wannier_setup(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice, &
								kpt_latt,num_bands_tot,num_atoms,atom_symbols,atoms_cart, &
								gamma_only,spinors,nntot,nnlist,nncell,num_bands,num_wann,proj_site, &
								proj_l,proj_m,proj_radial,proj_z,proj_x,proj_zona, &
								exclude_bands,proj_s,proj_s_qaxis)
		!
		if( nntot > 6 )	write(*,*)"[w90setup]: warning to many nearest neighborurs choosen nntot=",nntot
		do n = 1, num_bands_tot
			if( exclude_bands(n) /= 0) write(*,'(a,i4,a,i4)')"[w90setup]: warning exclude_bands(",n,")=",exclude_bands(n)
		end do
		!			
		return
	end subroutine




	subroutine w90prepEigVal(En)
		!convert eigenvalues from atomic units to eV
		real(dp),		intent(in)		:: En(:,:)
		integer							:: qi, n
		!
		if(	size(En,1) /= num_bands	) write(*,*)"[w90prepEigVal]: warning En has wrong numbers of bands"
		if(	size(En,2) /= num_kpts 	) write(*,*)"[w90prepEigVal]: warning En has wrong numbers of kpts"
		!
		do qi = 1, num_kpts
			do n = 1, num_bands
				eigenvalues(n,qi)	= En(n,qi) * aUtoEv
			end do
		end do
		!
		return
	end subroutine






!M_MATRIX ROUTINES
	subroutine w90prepMmat(ck)
		complex(dp),	intent(in)		:: ck(:,:,:) 	 !ck(			nG		,	nBands  	,	nQ	)	
		integer							:: qi, nn, n, m
		complex(dp)						:: oLap
		real(dp)						:: gShift(2)
		!
		if(	size(ck,2)/= num_bands ) write(*,*)"[w90prepEigVal]: waring ab initio exp. coefficients have wrong number of bands"
		if(	size(ck,3)/= num_kpts )  write(*,*)"[w90prepEigVal]: waring ab initio exp. coefficients have wrong numbers of kpts"
		!
		M_matrix	= dcmplx(0.0_dp)
		!
		do qi = 1, num_kpts
			do nn = 1, nntot
				do m = 1 , num_bands
					do n = 1, num_bands
						!calc overlap of unks
						if( nncell(3,qi,nn)/= 0 ) then
							oLap	= dcmplx(0.0_dp)
							if(qi==1 .and. m==1 .and. n==1 ) write(*,*)	"[w90prepMmat]: oLap set to zero for nn=",nn
						else
							gShift(1)			= nncell(1,qi,nn) * 2.0_dp * PI_dp / aX
							gShift(2)			= nncell(2,qi,nn) * 2.0_dp * PI_dp / aX
							oLap				= UNKoverlap(n,m, qi, nnlist(qi,nn), gShift, ck)
						end if
						!
						!
						M_matrix(n,m,nn,qi)		=  oLap
					end do
				end do
			end do
		end do

		return
	end subroutine


	complex(dp) function UNKoverlap(n, m, qi, knb, gShift, ck)
		!HELPER for w90prepMmat
		!calculates the overlap between unk at qi and at a neigbhouring k point knb
		!	integration only over the first unit cell
		!
		integer,		intent(in)		:: n, m, qi, knb
		real(dp),		intent(in)		:: gShift(2)
		complex(dp),	intent(in)		:: ck(:,:,:)  !ck(			nG		,	nBands  	,	nQ	)		
		integer							:: gi, gj, cnt, tot
		real(dp)						:: delta(2)
		!
		UNKoverlap	= dcmplx(0.0_dp)
		cnt	= 0
		tot	= 0


		do gi = 1, nGq(qi)
			do gj = 1, nGq(knb)
				delta(:)	=  ( Gvec(:,gi,qi)-qpts(:,qi) ) 	-  		( Gvec(:,gj,knb)-qpts(:,knb)-gShift(:) )
				if( norm2(delta) < machineP )	then
					UNKoverlap	= UNKoverlap +  dconjg( ck(gi,n,qi) ) * ck(gj,m,knb) 
					cnt = cnt + 1
				end if
				tot = tot + 1
			end do
		end do
		!
		!write(*,'(a,i8,a,i8)')	"[UNKoverlap]: used ",cnt," where nGmax(qi)=",nGq(qi)
		if( cnt > nGq(qi)	)	write(*,'(a,i8,a,i8)')	"[UNKoverlap]: warning, used ",cnt," where nGmax(qi)=",nGq(qi)
		if( cnt < nGq(qi) / 2.0_dp)	write(*,'(a,i8,a,i8)')	"[UNKoverlap]: warning, used  only",cnt," where nGmax(qi)=",nGq(qi)
		!
		return
	end function

	




!A_MATRIX ROUTINES
	subroutine w90prepAmat(ck)
		complex(dp),	intent(in)		:: ck(:,:,:) 	 !ck(			nG		,	nBands  	,	nQ	)
		integer							:: qi	
		!
		if(	size(ck,2)/= num_bands )	write(*,*)"[w90prepAmat]: waring ab initio exp. coefficients have wrong number of bands"
		if(	size(ck,3)/= num_kpts )  	write(*,*)"[w90prepAmat]: waring ab initio exp. coefficients have wrong numbers of kpts"
		if(	num_wann/nAt > 3) 			write(*,*)"[w90prepAmat]: can not handle more then 3 wfs per atom at the moment"
		!
		A_matrix	= dcmplx(0.0_dp)
		do qi = 1, nQ
			call calcAmatANA(qi,ck(:,:,qi), A_matrix(:,:,qi))
		end do
		!
		return 
	end subroutine


		subroutine calcAmatANA(qi,ckH, A)
		!analytic projection with hard coded integrals
		!	projection onto sin**2, sin cos, cos sin
		integer,		intent(in)	:: qi
		complex(dp),	intent(in)	:: ckH(:,:)
		complex(dp),	intent(out)	:: A(:,:) !A(nBands,nWfs)
		integer						:: n, m
		!
		A	= dcmplx(0.0_dp)
		!reminder: use dconjg on ckH
		!
		!
		!SINGLE ATOM
		if( nAt == 1 ) then
			do n = 1, nWfs
				do m = 1, nBands
					select case( n )
					case(1)
						A(m,n)	= g1Int(qi,m,1, ckH)
					case(2)
						A(m,n)	= g2Int(qi,m,1, ckH)
					case(3)
						A(m,n)	= g3Int(qi,m,1, ckH)
					case default
						write(*,*)"[calcAmatANA]: Warning hit default in single atom switch"
						A(m,n)	= dcmplx(0.0_dp)
					end select
				end do
			end do
		!DUAL ATOMS
		else if( nAt == 2 ) then
			do n = 1, nWfs
				do m = 1, nBands
					select case( n )
					case(1)
						A(m,n)	= g1Int(qi,m,1,ckH)
					case(2)
						A(m,n)	= g1Int(qi,m,2,ckH)
					case(3)
						A(m,n)	= g2Int(qi,m,1,ckH)
					case(4)
						A(m,n)	= g2Int(qi,m,2,ckH)
					case(5)
						A(m,n)	= g3Int(qi,m,1,ckH)
					case(6)
						A(m,n)	= g3Int(qi,m,2,ckH)
					case default
						write(*,*)"[calcAmatANA]: Warning hit default in dual atom switch"
						A(m,n)	= dcmplx(0.0_dp)
					end select
				end do
			end do
		!FALLBACK
		else
			write(*,*)	"[calcAmatANA]: more then two atoms per unit cell. I can only handle two tough"
		end if 
		

		!
		return
	end subroutine


	complex(dp) function g1Int(qi,m, at,ckH)
		integer,		intent(in)	:: qi, m, at
		complex(dp),	intent(in)	:: ckH(:,:)
		complex(dp)					:: num1, num2, denom
		real(dp)					:: kappa, xL, xR, yL, yR, Gx, Gy, xc, yc
		integer						:: gi
		!
		!TRIAL ORBITAL:
		kappa	= PI_dp / ( 2.0_dp * atR(1,at) )
		if( atR(1,at) /= atR(2,at) ) write(*,*)"[g1Int]: warning analytic projection can not handle non cubic wells"
		xc		= atPos(1,at) - atR(1,at)
		yc		= atPos(2,at) - atR(2,at)
		xL 		= atPos(1,at) - atR(1,at) 
		xR		= atPos(1,at) + atR(1,at)
		yL		= atPos(2,at) - atR(2,at)
		yR		= atPos(2,at) + atR(2,at)
		!
		!SUMMATION OVER G:
		g1Int 	= dcmplx(0.0_dp)
		do gi = 1, nGq(qi)
			Gx 		= Gvec(1,gi,qi)
			Gy		= Gvec(2,gi,qi)
			!
			!
			if( 	abs(Gx) < machineP 	.and.	 abs(Gy) < machineP 	) then
				num1	= dcos((xL-xc)*kappa) - dcos((xR-xc)*kappa)
				num2	= dcos((yL-yc)*kappa) - dcos((yR-yc)*kappa)
				denom	= kappa**2
			!
			else if( abs(Gy) < machineP ) then
				num1	= 	myExp(-Gx*(xL+xR)) 	* 	( 			dcos((yL-yc)*kappa) - 				dcos((yR-yc)*kappa)	)
				num2	= 			myExp(Gx*xR)* 	( kappa * 	dcos((xL-xc)*kappa) + i_dp * Gx * 	dsin((xL-xc)*kappa)	)
				num2	= num2  - 	myExp(Gx*xL)* 	( kappa *	dcos((xR-xc)*kappa) + i_dp * Gx *	dsin((xR-xc)*kappa)	)
				denom	= kappa * ( kappa**2 - Gx**2)
			!
			else if( abs(Gx) < machineP) then
				num1	= 	myExp(-Gy*(yL+yR))	*	(			dcos((xL-xc)*kappa) -				dcos((xR-xc)*kappa)	)
				num2	=		-	myExp(Gy*yR)*	( kappa *	dcos((yL-yc)*kappa) + i_dp * Gy *	dsin((yL-yc)*kappa)	)
				num2	= num2  +	myExp(Gy*yL)*	( kappa *	dcos((yR-yc)*kappa) + i_dp * Gy *	dsin((yR-yc)*kappa)	)
				denom	= kappa * ( Gy**2 - kappa**2 )
			!
			else
				num1 	=    	- myExp(-Gx*xL)	* 	( kappa *	dcos((xL-xc)*kappa)	+ 	i_dp * Gx * dsin((xL-xc)*Kappa) )
				num1 	= num1  + myExp(-Gx*xR) * 	( kappa * 	dcos((xR-xc)*kappa)	+	i_dp * Gx * dsin((xR-xc)*Kappa)	)
				num2 	= 		  myExp(Gy*yL) 	* 	( kappa *	dcos((yR-yc)*kappa)	+	i_dp * Gy * dsin((yR-yc)*kappa)	) 
				num2 	= num2  - myExp(Gy*yR) 	* 	( kappa *	dcos((yL-yc)*kappa)	+	i_dp * Gy * dsin((yL-yc)*kappa)	)
				denom	= myExp(Gy*(yL+yR)) * (Gy**2-kappa**2) * (Gx**2-kappa**2)
			!
			end if
			!
			!
			if( abs(denom) < machineP) then
				write(*,'(a,i4,a,f8.4,a,f8.4a,f8.4)') "[g1Int]: warning zero denom at qi=",qi," Gx=",Gx,"Gy=",Gy,"denom=",denom
				denom 	= 1e-5_dp
			end if
			!
			!
			g1Int = g1Int + dconjg(ckH(gi,m)) * num1 * num2 / (dsqrt(vol) * denom)
		end do
		!
		!
		return
	end function


	complex(dp) function g2Int(qi,m, at,ckH)
		integer,		intent(in)	:: qi, m, at
		complex(dp),	intent(in)	:: ckH(:,:)
		complex(dp)					:: num1, num2, denom
		real(dp)					:: kappa, xL, xR, yL, yR, Gx, Gy, xc, yc
		integer						:: gi
		!
		!TRIAL ORBITAL:
		kappa	= PI_dp / (2.0_dp*atR(1,at))
		if( atR(1,at) /= atR(2,at) ) write(*,*)"[g2Int]: warning analytic projection can not handle non cubic wells"
		xc		= atPos(1,at) - atR(1,at)
		yc		= atPos(2,at) - atR(2,at)
		xL 		= atPos(1,at) - atR(1,at) 
		xR		= atPos(1,at) + atR(1,at)
		yL		= atPos(2,at) - 2.0_dp * atR(2,at)
		yR		= atPos(2,at) + 2.0_dp * atR(2,at)
		!
		!SUMMATION OVER G:
		g2Int 	= dcmplx(0.0_dp)
		do gi = 1, nGq(qi)
			!
			Gx 		= Gvec(1,gi,qi)
			Gy		= Gvec(2,gi,qi)
			!
			!
			if( abs(Gx) < machineP .and. abs(Gy) < machineP ) then
				num1	= -  ( 	dcos((xL-xc)*kappa) - dcos((xR-xc)*kappa) )
				num2	= 	   	dsin((yL-yc)*kappa) - dsin((yR-yc)*kappa)
				denom	= 	kappa**2
			!
			else if( abs(Gy) < machineP ) then
				num1	=		  myExp(-Gx*(xL+xR))*	(			dsin((yL-yc)*kappa) - 				dsin((yR-yc)*kappa)	)
				num2	= 		- myExp(-Gx*xR)		*  	( kappa *	dcos((xL-xc)*kappa) +	i_dp * Gx * dsin((xL-xc)*kappa)	)
				num2	= num2 	+ myExp(Gx*xL)		*	( kappa *	dcos((xR-xc)*kappa) +	i_dp * Gx * dsin((xR-xc)*kappa)	)
				denom	= kappa * ( kappa**2 - Gx**2 )
			!
			else if( abs(Gx) < machineP ) then
				num1	= 		  myExp(-Gy*(yL+yR))*	(			dcos((xL-xc)*kappa) -				dcos((xR-xc)*kappa)	)
				num2	=		  myExp(Gy*yR)		*	( kappa *	dsin((yL-yc)*kappa) -	i_dp * Gy * dcos((yL-yc)*kappa)	)
				num2	= num2 	+ myExp(Gy*yL)		*	(-kappa *	dsin((yR-yc)*kappa) +	i_dp * Gy * dcos((yR-yc)*kappa)	)
				denom	= kappa * ( Gy**2 - kappa**2 )
			!
			else  
				num1 	=		- myExp(-Gx*xL)		* 	( kappa*	dcos((xL-xc)*kappa) + 	i_dp * Gx * dsin((xL-xc)*Kappa)	)
				num1 	= num1	+ myExp(-Gx*xR)		* 	( kappa*	dcos((xR-xc)*kappa) +	i_dp * Gx * dsin((xR-xc)*Kappa)	)
				num2 	=		+ myExp(Gy*yR) 		* 	( kappa*	dsin((yL-yc)*kappa) -	i_dp * Gy * dcos((yL-yc)*kappa)	)
				num2 	= num2	+ myExp(Gy*yL) 		* 	(-kappa* 	dsin((yR-yc)*kappa) +	i_dp * Gy *	dcos((yR-yc)*kappa)	) 
				denom	=		+ myExp(Gy*(yL+yR)) * (Gy**2-kappa**2) * (Gx**2-kappa**2)  
			end if
			!
			!
			if( abs(denom) < machineP) then
				write(*,'(a,i4,a,f8.4,a,f8.4a,f8.4)') "[g2Int]: warning zero denom at qi=",qi," Gx=",Gx,"Gy=",Gy,"denom=",denom
				denom 	= 1e-5_dp
			end if
			!
			!
			g2Int = g2Int + dconjg(ckH(gi,m)) * num1 * num2 / (dsqrt(vol) * denom)
		end do
		!
		!
		return
	end function

	complex(dp) function g3Int(qi,m, at,ckH)
		integer,		intent(in)	:: qi, m, at
		complex(dp),	intent(in)	:: ckH(:,:)
		complex(dp)					:: num1, num2, denom
		real(dp)					:: kappa, xL, xR, yL, yR, Gx, Gy, xc, yc
		integer						:: gi
		!
		kappa	= PI_dp / (2.0_dp*atR(1,at))
		if( atR(1,at) /= atR(2,at) ) write(*,*)"[g3Int]: warning analytic projection can not handle non cubic wells"
		xc		= atPos(1,at) - atR(1,at)
		yc		= atPos(2,at) - atR(2,at)
		xL 		= atPos(1,at) - 2.0_dp * atR(1,at) 
		xR		= atPos(1,at) + 2.0_dp * atR(1,at)
		yL		= atPos(2,at) - atR(2,at)
		yR		= atPos(2,at) + atR(2,at)
		g3Int 	= dcmplx(0.0_dp)
		do gi = 1, nGq(qi)
			!
			Gx 		= Gvec(1,gi,qi)
			Gy		= Gvec(2,gi,qi)
			!
			!
			if( abs(Gx) < machineP .and. abs(Gy) < machineP ) then
				num1	= - ( 	dcos((yL-yc)*kappa) - dcos((yR-yc)*kappa) 	)
				num2	= 		dsin((xL-xc)*kappa) - dsin((xR-xc)*kappa)
				denom	= kappa**2
			!
			else if( abs(Gy) < machineP ) then
				num1	= 		  myExp(-Gx*(xL+xR))*	( 			dcos((yL-yc)*kappa) - 				dcos((yR-yc)*kappa) )
				num2	= 		  myExp(Gx*xR)		*	(-kappa* 	dsin((xL-xc)*kappa) + 	i_dp * Gx * dcos((xL-xc)*kappa)	)
				num2	= num2	+ myExp(Gx*xL)		*	( kappa*	dsin((xR-xc)*kappa) - 	i_dp * Gx * dcos((xR-xc)*kappa)	)
				denom	= kappa * ( kappa**2 - Gx**2)
			!
			else if( abs(Gx) < machineP ) then
				num1	= 		  myExp(-Gy*(yL+yR))*	(			dsin((xR-xc)*kappa) -				dsin((xL-xc)*kappa)	)
				num2	=		- myExp(Gy*yR)		*	( kappa*	dcos((yL-yc)*kappa) +	i_dp * Gy * dsin((yL-yc)*kappa)	)
				num2	= num2	+ myExp(Gy*yL)		*	( kappa*	dcos((yR-yc)*kappa) +	i_dp * Gy * dsin((yR-yc)*kappa)	)
				denom	= kappa * ( Gy**2 - kappa**2 )
			!		 		
			else
				!ToDo: revisit
				num1	=  		- i_dp * myExp(Gx*xR) * Gx 		* dcos((xL-xc)*kappa)
				num1	= num1 	+ i_dp * myExp(Gx*xL) * Gx 		* dcos((xR-xc)*kappa)
				num1	= num1 	+ 	    myExp(Gx*xR) * kappa 	* dsin((xL-xc)*kappa)
				num1	= num1 	- 		myExp(Gx*xL) * kappa	* dsin((xR-xc)*kappa)
				!
				num2	=  		- 		myExp(Gy*yR) *  ( kappa	* dcos((yL-yc)*kappa) + i_dp * Gy * dsin((yL-yc)*kappa) )
				num2	= num2 +		myExp(Gy*yL) * 	( kappa * dcos((yR-yc)*kappa) + i_dp * Gy * dsin((yR-yc)*kappa) )
				!
				denom	=  + myExp( Gx*(xL+xR) + Gy*(yL+yR) )	* (Gy**2-kappa**2) * (kappa**2-Gx**2)
				!  
			end if
			!
			!
			if( abs(denom) < machineP) then
				write(*,'(a,i4,a,f8.4,a,f8.4a,f8.4)') "[g3Int]: warning zero denom at qi=",qi," Gx=",Gx,"Gy=",Gy,"denom=",denom
				denom 	= 1e-5_dp
			end if
			!
			!
			g3Int = g3Int + dconjg(ckH(gi,m)) * num1 * num2 / (dsqrt(vol) * denom)
		end do
		!
		!
		return
	end function




	





















end module
