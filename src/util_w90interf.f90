module util_w90Interf	

	use util_math,	only:	dp, PI_dp, i_dp, machineP, aUtoAngstrm, aUtoEv, myExp, isUnit, isIdentity
	use util_sysPara
	use util_basisIO,		only:	read_energies, read_Amn, read_Mmn

	implicit none

	private
	public ::					setup_w90, write_w90_matrices, &
								read_M_initial, read_FD_b_vectors, read_U_matrix, read_band_interp, read_tb_basis, read_wann_centers, & 
								printNNinfo, &
								seed_name

	!public var
	character(len=3)					:: 	seed_name
	complex(dp),		allocatable		::	U_matrix(:,:,:)
	real(dp),			allocatable		::	wann_centres(:,:), wann_spreads(:) 
	real(dp)							::	spread(3)
	!private var
	character(len=20),	allocatable		:: 	atom_symbols(:)
	logical								:: 	gamma_only, spinors
	integer								:: 	mp_grid(3), num_kpts, num_bands_tot, num_atoms, num_nnmax, &
											nntot, num_bands, num_wann
	integer,			allocatable		:: 	nnlist(:,:), nncell(:,:,:), proj_l(:), proj_m(:), &
											proj_radial(:), exclude_bands(:), proj_s(:)
	real(dp)							:: 	real_lattice(3,3), recip_lattice(3,3)
	real(dp),			allocatable		:: 	kpt_latt(:,:), atoms_cart(:,:), &
											proj_site(:,:), proj_z(:,:), proj_x(:,:), proj_zona(:), proj_s_qaxis(:,:)

		


	contains




!PUBLIC WRITE
	subroutine setup_w90(nntot_out, nnlist_out, nncell_out)
		integer,		intent(out)		::	nntot_out, nnlist_out(:,:), nncell_out(:,:,:)
		
		!
		!PREP W90 INIT
		seed_name	= seedName
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
		real_lattice(3,3)	= min(aX,aY) * 0.99_dp * aUtoAngstrm !make sure z is smallest real, results in larger distance k-space nn
		!reciprocal cell
		recip_lattice		= 0.0_dp
		recip_lattice(1,1)	= 2.0_dp * PI_dp / real_lattice(1,1) 	
		recip_lattice(2,2)	= 2.0_dp * PI_dp / real_lattice(2,2)
		recip_lattice(3,3)	= 2.0_dp * PI_dp / real_lattice(3,3) 		
		!atoms
		num_bands_tot		= nBands
		num_atoms			= nAt

		!debug
		if(	abs(dqx-dqy) > 1e-7_dp ) 	stop	"[w90setup]: ERROR - q mesh has non uniform grid spacing (required by w90)"

		
		!W90 
		call run_w90setup(nntot_out, nnlist_out, nncell_out)
		write(*,*)	"[w90Interf]: done with w90 setup"
		write(*,*)	"[w90Interf]: will use num_bands= ",num_bands, " bands"
		write(*,*)	"[w90Interf]: to gen   num_wann=  ",num_wann, " wnfs"
	

		!WRITE INPUT FILE (wann run)
		call write_W90run_input()
		write(*,*)	"[w90Interf]: wrote input file for wannierisation (.win)"
		!PREP PW90
		call writeW90KinterpMesh()
		write(*,*)	"[w90Interf]: wrote interpolation mesh file (_geninterp.kpt)"
	
		!WANNIERISE
		!call wannier_run(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice, &
		!					kpt_latt,num_bands,num_wann,nntot,num_atoms,atom_symbols, &
		!					atoms_cart,gamma_only,M_matrix_orig,A_matrix,eigenvalues, &
		!					U_matrix,U_matrix_opt,lwindow,wann_centres,wann_spreads, &
		!					spread)

		!note: the library modus gives different results compared to executing w90 manually, why is currently unclear

		!inquire(file="../thirdparty/wannier90/wannier90.x ", exist=w90success)
		!if( w90success ) then
		!	call chdir(w90_Dir)
		!	call execute_command_line ("../../thirdparty/wannier90/wannier90.x "//seed_name, exitstat=ierr)
 	 	!	print *, "Exit status of wannier90 was ", ierr
		!	call chdir("..")
		!	write(*,*)	"[w90Interf]: wannierization finished"
		!else
		!	write(*,*)	"[w90Interf]: please provide wannier90.x in /thirdpary/wannier90/ or execute wannier90 manually"
		!end if
		!
		!
	end subroutine


	subroutine write_w90_matrices()
		complex(dp),	allocatable		::	Amn(:,:,:), Mmn(:,:,:,:)
		real(dp),		allocatable		::	En(:,:)
		integer							::	qi, nn, n, m 
		!
		allocate(	Amn(nBands, nWfs, 			nQ	)	)
		allocate(	Mmn(nBands, nBands, nntot, 	nQ	)	)
		allocate(	En(	nSolve,					nQ	)	)
		!
		call read_Amn(Amn)
		call read_Mmn(Mmn)
		call read_energies(En)
		!
		!WRITE TO FILE
		!Mmn
		open(unit=120,file=w90_Dir//seed_name//'.mmn',action='write',access='stream',form='formatted', status='replace')
		write(120,*)	'overlap matrix'
		write(120,*)	num_bands, ' ', num_kpts, ' ', nntot	
		!
		do qi = 1, num_kpts
			do nn = 1, nntot
				write(120,*)	qi,' ',nnlist(qi,nn),' ',nncell(1,qi,nn),' ',nncell(2,qi,nn),' ', nncell(3,qi,nn)
				do n = 1, num_bands
					do m = 1, num_bands
						write(120,*)	dreal(Mmn(m,n,nn,qi)), ' ', dimag(Mmn(m,n,nn,qi))
					end do
				end do		
			end do
		end do
		close(120)
		!
		!Amn
		open(unit=120,file=w90_Dir//seed_name//'.amn',action='write',access='stream',form='formatted', status='replace')
		write(120,*)	'projection matrix( initial guess for U matrix)'
		write(120,*)	num_bands, ' ', num_kpts, ' ', num_wann
		do qi = 1, num_kpts
			do n = 1, num_wann
				do m = 1, num_bands	
					write(120,*)	m, ' ', n, ' ', ' ', qi, ' ', dreal(Amn(m,n,qi)), ' ', dimag(Amn(m,n,qi))
				end do
			end do
		end do
		!
		!En
		open(unit=120,file=w90_Dir//seed_name//'.eig',action='write',access='stream',form='formatted', status='replace')
		do qi = 1, nQ
			do n = 1, nBands
				write(120,*)	n, ' ', qi, ' ', aUtoEv * En(n,qi)
			end do
		end do
		close(120)
		!
		!
		return
	end subroutine



	subroutine printNNinfo(nntot, nnlist)
		integer, 	intent(in)				::	nntot, nnlist(:,:)
		integer								::	gammaPt, nn
		real(dp)							::	w_b_guess
		!
		write(*,'(a,i2,a)')	"[berryMethod]: nn info:"
		gammaPt = 1 + int(	nQx*(0.5_dp+0.5_dp*nQy)	)
		write(*,'(a,f6.2,a,f6.2,a)')	"        dqx=",dqx,"; dqy=",dqy,"."
		write(*,'(a,f6.2,a,f6.2,a)')	"        this means for qpt=(",qpts(1,gammaPt),", ",qpts(2,gammaPt),")."
		write(*,*)	" nn  | q_nn(x) | q_nn(y) | w_b "
		write(*,*)	"-------------------------------"
		do nn = 1, nntot
			w_b_guess = 2.0_dp / ( real(nntot,dp) * dqx**2 )
			write(*,'(i2,a,f6.2,a,f6.2,a,f6.2)')	nn,"  |  ",	qpts(1,nnlist(gammaPt,nn)),"  |  ",&
																				qpts(2,nnlist(gammaPt,nn)),"  |  ",w_b_guess
		end do
		!
		return
	end subroutine 












!PUBLIC READ
	subroutine read_M_initial(f_num_bands, f_num_kpts, f_nntot, f_nnlist, f_nncell, M_init)
		!read the M_matrix written by w90prepMmat	
		integer,						intent(out)	::	f_num_bands, f_num_kpts, f_nntot
		integer,		allocatable,	intent(out)	::	f_nnlist(:,:), f_nncell(:,:,:)
		complex(dp),	allocatable,	intent(out)	::	M_init(:,:,:,:)
		integer										::	qi, f_qi, nn, n, m
		real(dp)									::	realBuff(2)
		logical										:: 	foundFile
		!
		!check for file
		inquire(file=w90_Dir//seedName//'.mmn', exist=foundFile)
		if( .not. foundFile ) stop 'the .mmn file could not be found'
		!
		open(unit=120,file=w90_Dir//seedName//'.mmn',action='read',access='stream',form='formatted', status='old')
		read(120,*)
		read(120,*) f_num_bands, f_num_kpts, f_nntot		
		!		
		allocate(	f_nnlist(			f_num_kpts, 	f_nntot					)	)
		allocate(	f_nncell(	3,		f_num_kpts, 	f_nntot					)	)
		allocate(	M_init(f_num_bands, f_num_bands, 	f_nntot, f_num_kpts		)	)
		!
		f_nnlist	= 0
		f_nncell	= 0
		M_init 		= dcmplx(0.0_dp)
		!
		do qi = 1, f_num_kpts
			do nn = 1, f_nntot
				!
				read(120,*)	f_qi, f_nnlist(qi,nn), f_nncell(1:3,qi,nn)
				if(	 qi /= 	f_qi )						stop	"[read_M_initial]: WARNING q mesh ordered differently"
				!
				do n = 1, f_num_bands
					do m = 1, f_num_bands
						read(120,*)	realBuff(1:2)
						M_init(m,n,nn,qi)	= dcmplx( realBuff(1), realBuff(2) )
					end do
				end do
				!
			end do
		end do
		!
		!DEBUG
		if( f_nntot 	/= 4	)	write(*,*) 	"[read_M_initial]: WARNING more then 4 nn in .mmn file"
		if(	f_num_kpts 	/= nQ 	) 	stop 		"[read_M_initial]: .mmn file defined on different qpt mesh"
		if(	f_num_bands /= nBands)	stop		"[read_M_initial]: nBands not matching"
		write(*,*)	"[read_M_initial]: read the .mmn file"
		!
		!WRITE CLONE FOR DEBUG
		!open(unit=125,file=w90_Dir//'clone'//'.mmn',action='write',access='stream',form='formatted', status='replace')
		!write(125,*)	"clone of the "//seedName//".mmn file"
		!write(125,*)	f_num_bands, f_num_kpts, f_nntot
		!do qi = 1, f_num_kpts
		!	do nn = 1, f_nntot
		!		write(125,*)	qi, f_nnlist(qi,nn),	f_nncell(1:3,qi,nn)
		!		do n = 1, f_num_bands
		!			do m = 1, f_num_bands
		!				write(125,*)	M_init(m,n,nn,qi)
		!			end do
		!		end do
		!	end do
		!end do
		!close(125)
		!write(*,*)	"[read_M_initial]: wrote clone.mmn file"
		!
		return
	end subroutine


	subroutine read_FD_b_vectors(b_k, w_b)
		!call after w90 is finished
		!reads seedname.nnkp & seedname.wout
		real(dp),			intent(out)		::	b_k(:,:), w_b(:)
		integer								::	stat, qnn,  nn, start, end
		real(dp)							::	b_vec(3), weight
		character(len=*), 	parameter 		::	search_wout=" |                  b_k Vectors (Ang^-1) and Weights (Ang^2)                  |"
		character(len=100)					::	line
		logical								::	finished
		!
		!READ .wout
		finished = .false.
		open(unit=335, iostat=stat, file=w90_dir//seedName//'.wout', form='formatted', status='old', action='read')
		do while( .not. finished .and. stat==0)
			read(335,"(a)",iostat=stat) line
			!find the block
			if( line == search_wout) then
				read(335,*)
				read(335,*)
				read(335,*)
				do qnn = 1, size(w_b)
					!read line
					read(335,"(a)",iostat=stat) line
					!remove "|" from line  (this is neccessary to use pattern matching on string)
					start 	= scan(line,"|" ) 			+1
					end 	= scan(line,"|", .true.)	-1
					!extract values from substring
					read(line(start:end),*)	nn, b_vec, weight 
					!convert to a.u.
					b_k(1:3, nn) 	= b_vec(1:3)!*aUtoAngstrm			![b_k] = (Ang^-1)
					w_b(nn)			= weight !/ aUtoAngstrm**2			![w_b] = (Ang^2)
					!debug
					!write(*,'(i2,a,f6.4,a,f6.4,a,f6.4,a,f6.4)') nn," ",b_k(1,nn)," ",b_k(2,nn)," ",b_k(3,nn)," ",w_b(nn)
					if( qnn /= nn ) write(*,*) "[readFDscheme]: WARNING, possibly mismatched nearest neighbours"
				end do
				finished = .true.
			end if
			!
		end do
		close(335)
		!
		write(*,*)		" 	nn | 	b[Å^-1]	| 	w_b[Å^2] "
		do nn = 1, size(w_b)
			write(*,'(a,i2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)')		"  ",nn," | (",b_k(1,nn),", ",b_k(2,nn),", ",b_k(3,nn),") |    ",w_b(nn)
		end do


		!DEBUG
		if( 		finished ) write(*,*)	"[readFDscheme]: read the .wout file"
		if(	.not. 	finished ) write(*,*)	"[readFDscheme]: did not find search string: '",search_wout,"' in .wout file"
		!
		if( .not. B1condition(b_k, w_b)	)	stop	"[read_FD_scheme]: B1condition was not satisfied"
		!
		!
		return
	end subroutine




	subroutine	read_wann_centers(w_centers)
		real(dp),		intent(out)			:: 	w_centers(:,:)
		integer								::	wannF, ntot, at, stat, start
		character(len=100)					::	line
		logical								::	file_exists
		!
		!inquire file
		inquire(file=w90_dir//seedName//'_centres.xyz', exist=file_exists)
		if( .not. file_exists ) stop '[read_wann_centers]: could not find centres.xyz file'
		!
		!read file
		open(unit=340, iostat=stat, file=w90_dir//seedName//'_centres.xyz', form='formatted', status='old', action='read')
		if( stat == 0) then
			!read header
			read(340, *)	 ntot
			read(340, *)	
			if( ntot - nAt /= nWfs )	write(*,*)"[read_wann_centers]: WARNING _centres.xyz has wrong nWfs"
			!
			!read centers
			do wannF = 1,  size(w_centers,2)
				read(340,"(a)",iostat=stat) line
				start 	= scan(line,"X" ) 			+1
				!read everything after "X"
				read(line(start:100), *)	w_centers(1:3, wannF)
			end do
			!make sure correct amount of atoms is given (could also read atom positions here)
			do at = 1, nAt
				read(340, *)
			end do
		else
			write(*,*)	"[read_wann_centers]: could not read _centres.xyz file"
			w_centers = 0.0_dp
		end if
		close(340)
		!
		!convert to [a.u.]
		w_centers 	= w_centers 	!/ aUtoAngstrm
		write(*,*)	"[read_wann_centers]: read the .xyz file"
		!
		!
		return
	end subroutine



	subroutine read_U_matrix(f_num_wann, U_matrix)
		!READ U FROM W90
		integer,						intent(out)	::	f_num_wann
		complex(dp),	allocatable,	intent(out)	:: 	U_matrix(:,:,:)
		real(dp),		allocatable					:: 	krel(:,:)
		integer										:: 	f_num_kpts, stat, qi, n, m, dumI(3)
		real(dp)									:: 	val(2)
		logical										::	file_exists
		!
		!inquire file
		inquire(file=w90_dir//seedName//'_u.mat',exist=file_exists)
		if(.not. file_exists )	stop '[Umat_reader]: did not find U matrix'
		!
		!read file
		open(unit=300,iostat=stat, file=w90_dir//seedName//'_u.mat', status='old',action='read')
		if( stat /= 0)	write(*,*)	"[readUmatrix]: WARNING did not file _u.mat file"
		!
		!header 
		read(300,*)
		read(300,*) dumI(1:3)
		f_num_kpts 	= dumI(1)
		f_num_wann	= dumI(2)

		if(	dumI(3)	/= f_num_wann) 			stop	"[readUmatrix]:  ERROR - file specifies two differnt num_wann values"
		

		allocate(	krel( 		3,f_num_kpts						)	)
		allocate(	U_matrix(	f_num_wann, f_num_wann, f_num_kpts	)	)
		!
		!body
		U_matrix	= dcmplx(0.0_dp)
		do qi = 1,  size(U_matrix,3)
			read(300,*)
			read(300,*) krel(1:3,qi)
			do n = 1, size(U_matrix,2)
				do m = 1, size(U_matrix,1)
					read(300,*)	val(1:2)
					U_matrix(m,n,qi)	= dcmplx(val(1), val(2) )	
				end do
			end do
			!DEBUG
			if( abs( krel(1,qi)*2.0_dp*PI_dp/aX - qpts(1,qi)) > machineP) then
					write(*,*)	"[readUmatrix]:	kx: k_w90= ",krel(1,qi)*2.0_dp*PI_dp/aX, " qpts=",qpts(1,qi)
					write(*,*)	"[readUmatrix]:	ky: k_w90= ",krel(2,qi)*2.0_dp*PI_dp/aY, " qpts=",qpts(2,qi)
					write(*,*)	"[readUmatrix]:	WARNING k meshes are ordered diffferently"
			end if
		end do
		close(300)
		!
		!DEBUG TEST
		do qi = 1, size(U_matrix,3)
			if(			isIdentity(U_matrix(:,:,qi))		) write(*,*)	"[read_U_matrix]: U_matrix is identity at qi=", qi
			if(		.not.	 isUnit(U_matrix(:,:,qi))		) stop "[read_U_matrix]: found non unitary matrix in input file..."
		end do
		!FINALIZE
		write(*,*)	"[read_U_matrix]: read the _u.mat file"
		!
		!
		return
	end subroutine





	subroutine read_band_interp(kpts, en, v_vec )
		real(dp),		intent(out)		::	kpts(:,:), en(:,:), v_vec(:,:,:)
		real(dp)						::	buffer(7)
		integer							::	stat, qi, n, qInd
		logical							::	file_exists
		!
		!write(*,*)  seedName//'_geninterp.dat'
		!check if exists
		inquire(file=w90_dir//seedName//'_geninterp.dat', exist=file_exists)
		if( .not. file_exists ) stop 'geninterp.dat file not found'
		!read
		open(unit=320,iostat=stat, file=w90_dir//seedName//'_geninterp.dat',form='formatted', status='old',action='read')
		read(320,*)
		read(320,*)
		read(320,*)
		do qi = 1, nK
			do n = 1, nWfs
				read(320,*)	qInd, buffer
				if( n==1 )	kpts(1:3,qi) = buffer(1:3) * aUtoAngstrm  
				en(n,qInd)		= buffer(4)
				v_vec(1,n,qInd)	= buffer(5)
				v_vec(2,n,qInd)	= buffer(6)
				v_vec(3,n,qInd)	= buffer(7) 
			end do
		end do
		close(320)
		!
		!ATOMIC UNITS CONVERSION:
		en		= en	/ (	aUtoEv)
		v_vec 	= v_vec  / (aUtoEv  * aUtoAngstrm)	 ! [v_vec] = eV / Angstroem
		write(*,*)	"[readBandVelo]: read the geninterp.dat file"
		!
		!	
		return
	end subroutine



	subroutine read_tb_basis( Rcell, tHopp, rHopp)
		real(dp),		intent(out), allocatable			::	Rcell(:,:)
		complex(dp),	intent(out), allocatable			:: 	tHopp(:,:,:), rHopp(:,:,:,:)
		integer												::	stat, f_nwfs, f_nSC, readLines, line, &
																cell, n
		integer												::	cell_rel(3), index(2)
		real(dp)											::	unit_cell(3,3), compl1(2),compl3(6), rTest(3), real3(3)

		

		open(unit=330,iostat=stat, file=w90_dir//seedName//'_tb.dat',form='formatted', status='old', action='read')

		!read unit cell
		read(330,*)
		read(330,*)	real3
		unit_cell(1,1:3)	= real3
		read(330,*)	 real3
		unit_cell(2,1:3) 	= real3
		read(330,*)	 real3
		unit_cell(3,1:3)	= real3
		unit_cell = unit_cell / aUtoAngstrm


		
		if( abs(unit_cell(1,1)-aX) > 1e-8_dp) 	write(*,*)	"[read_tb_basis]: WARNING problem reading aX"
		if( abs(unit_cell(2,2)-aY) > 1e-8_dp) 	write(*,*)	"[read_tb_basis]: WARNING problem reading aY"
		if( unit_cell(3,3) > unit_cell(1,1) .or. unit_cell(3,3)  > unit_cell(2,2) )	write(*,*)	"[read_tb_basis]: WARNING aZ is to large"

		!read basis size
		read(330,*)	f_nwfs
		read(330,*)	f_nSC

		!allocate targets
		allocate(	Rcell(3,						f_nSC		)	)
		allocate(	rHopp(3,	f_nwfs, f_nwfs,		f_nSC		)	)
		allocate(	tHopp(		f_nwfs,	f_nwfs,		f_nSC		)	)


		
		!debug warnings
		if(	f_nwfs 	/= nWfs )	write(*,'(a,i4,a,i4,a)')	"[read_tb_basis]: WARNING wrong nWfs detected: got ",f_nwfs," (expected ",nWfs,")."
		if(	f_nSC	/= nSC )	write(*,'(a,i4,a,i4,a)')	"[read_tb_basis]: WARNING wrong nSC detected: got ",f_nSC," (expected ",nSC,")."

		!dummy read the degeneracy list
		readLines = ceiling( real(f_nSC,dp) / real(15,dp)		)
		do line = 1, readLines
			read(330,*)
		end do

		!read tHopp (eV)
		do cell = 1, nSC
			read(330,*)
			read(330,*) cell_rel(1:3)

			Rcell(1:3,cell)	= matmul(unit_cell,	real(cell_rel(1:3),dp)	)
			do n = 1, nWfs**2
				read(330,*)	index(1:2), compl1(1:2)
				!
				tHopp(index(1),index(2), cell)	= dcmplx(	compl1(1), compl1(2)		)
			end do
		end do
		

		!read rHopp (ang)
		do cell = 1, nSC
			read(330,*)
			read(330,*)	cell_rel(1:3)

			rTest(1:3)	= matmul(unit_cell,	real(cell_rel(1:3),dp)	)
			if(		 norm2( rTest(1:3) - Rcell(1:3,cell) )	> 1e-8_dp		) write(*,*)	"[read_tb_basis]: WARNING rHopp has different sc order then tHopp"

			do n = 1, nWfs**2
				read(330,*)	index(1:2), compl3(1:6)

				rHopp(1,	index(1), index(2), cell)	= dcmplx( compl3(1), compl3(2)	)
				rHopp(2,	index(1), index(2), cell)	= dcmplx( compl3(3), compl3(4)	)
				rHopp(3,	index(1), index(2), cell)	= dcmplx( compl3(5), compl3(6)	)
			end do
		end do 
		!
		close(330)
		!
		!convert to a.u.
		tHopp	= 	tHopp / aUtoEv
		rHopp	= 	rHopp / aUtoAngstrm
		!
		!
		return
	end subroutine 








































!prviate
!PREPARE & RUN W90
	subroutine write_W90setup_input()
		!write input file for wannier_setup call
		integer	:: at, n_per_at, nWf, i
		character(len=30)	:: 	orbitals
		!
		open(unit=100,file=seed_name//'.win',action='write', form='formatted', status='replace')
		!
		!BASIC INFO
		write(100,*)	'num_wann  = ',nWfs
		write(100,*)	'num_bands = ',nBands
		write(100,*)	'mp_grid   = ', mp_grid(1) , ' ', mp_grid(2), ' ', mp_grid(3)
		if( useBloch )	write(100,*)	'use_bloch_phases = true '
		write(100,*)	
		!
		!FINITE DIFFERENCE
		!write(100,*)	'search_shells = 200'
		if(	nShells==1 )	write(100,*)	'shell_list =',	shells(1)
		if( nShells > 1)	write(100,*)	'shell_list = ', (shells(i), i=1,size(shells) )
		write(100,*)	'skip_b1_tests = .true.'
		!

		! k point neighbours
		write(100,*)	'postproc_setup = .true.'
		write(100,*)
		!
		!PROJECTIONS
		write(100,*)	'Begin Projections'
		do at = 1, size(atom_symbols)
			!get states living at current atom
			n_per_at	= 0
			do nWf = 1, nWfs
				if( proj_at(nWf) == at ) n_per_at = n_per_at + 1
			end do

			!generate generic projection strintg
			call randomProjString(n_per_at, orbitals)

			!write to input
			write(100,*)	atom_symbols(at),' : '//orbitals
		end do
		write(100,*)	'End Projections'
		
		!
		!
		close(100)
		return
	end subroutine


	subroutine randomProjString(nStates,string)
		integer,			intent(in)			:: nStates
		character(len=30), 	intent(out)			:: string
		!
		select case(nStates)
			case (1)
				string = 'l=0' 								! s
			case (2)
				string = 'l=0; l=1, mr=2'	 				! s ; px
			case (3)
				string = 'l=0; l=1, mr=2,3'					! s ; px py
			case (4)
				string = 'l=0; l=1, mr=2,3; l=2, mr=4'		! s	; px py ; dx**2-y**2
			case (5)
				string = 'l=0; l=1, mr=2,3; l=2, mr=4,5'	! s ; px py ; dx**2-y**2 dxy
			case (6)
				string = 'l=0; l=1; l=2, mr=4,5'			!random
			case (7)
				string = 'l=0; l=1; l=2, mr=3,4,5'			!random
			case (8)
				string = 'l=0; l=1; l=2, mr=2,3,4,5'		!random
			case (9)
				string = 'l=0; l=1; l=2'					!random
			case default
				string = ''
				write(*,*)	"[randomProjString]: only defined up to 9 states per atom at the moment, was called with n=",nStates
		end select 
		!
		return
	end subroutine


	subroutine run_w90setup(nntot_out, nnlist_out, nncell_out)
		integer,		intent(out)		::	nntot_out, nnlist_out(:,:), nncell_out(:,:,:)
		integer							:: 	at, n
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
			!if( mod(at,2) == 0)	atom_symbols(at)	= 'H '
			!if(	mod(at,2) == 1)	atom_symbols(at)	= 'He'
			!
			!atom_symbols(at)	= 'H'
 	
 			if(at<10)	then
 				write (atom_symbols(at), "(A1,I1)") "H", at
 	 		else if( at<100) then
 	 			write (atom_symbols(at), "(A1,I2)") "H", at
 	 		else
 	 			write(atom_symbols(at),"(A1,I19)") "H",at
 	 		end if

 	 		write(*,*)	'atom_symbols(',at,')=',atom_symbols(at)
 			atoms_cart(1:2,at)	= atPos(1:2,at)*aUtoAngstrm
			atoms_cart(3,at)	= 0.0_dp
		end do
		!kpt lattice (fractional co-ordinates relative to the reciprocal lattice vectors)
		kpt_latt(1,:)			= qpts(1,:) / ( recip_lattice(1,1) * aUtoAngstrm ) 
		kpt_latt(2,:)			= qpts(2,:)	/ ( recip_lattice(2,2) * aUtoAngstrm ) 
		kpt_latt(3,:)			= 0.0_dp
		!
		
		!WRITE INPUT FILE (wann setup)
		call write_W90setup_input()
		write(*,*)	"[w90Interf]: wrote w90 setup input file (.win)"
		!
		write(*,*)	"[run_w90setup]: try to call the wannier90 library mode (wannier_setup)"
		!	

		call wannier_setup(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice, &
								kpt_latt,num_bands_tot,num_atoms,atom_symbols,atoms_cart, &
								gamma_only,spinors,nntot,nnlist,nncell,num_bands,num_wann,proj_site, &
								proj_l,proj_m,proj_radial,proj_z,proj_x,proj_zona, &
								exclude_bands,proj_s,proj_s_qaxis)
		!
		write(*,'(a,i3,a)')	"[run_w90setup]: success, obtained fd sheme( nntot=",nntot,")"
		!
		nntot_out 	= nntot
		nnlist_out	= nnlist
		nncell_out	= nncell
		!DEBUG:
		if( nntot /= 4 )		write(*,*)	"[w90setup]: WARNING did not find exactly 4 nearest neighbours"
		do n = 1, num_bands_tot
			if( exclude_bands(n) /= 0) write(*,'(a,i4,a,i4)')"[w90setup]: WARNING exclude_bands(",n,")=",exclude_bands(n)
		end do
		!
		!			
		return
	end subroutine



	subroutine write_W90run_input()
		!input for wannierisation
		integer				:: qi, at
		!
		!create new
		open(unit=100,file=w90_Dir//seed_name//'.win',action='write',access='stream',form='formatted', status='replace')
		!
		!BASIC INFO
		write(100,*)	'num_wann  = ',nWfs
		write(100,*)	'mp_grid   = ', mp_grid(1) , ' ', mp_grid(2), ' ', mp_grid(3)
		write(100,*)	'num_iter  = ', nW90it
		if( nW90it >= 100 )		write(100,*)	'num_print_cycles	=', (nw90it / 10)
		 
		write(100,*)	
		write(100,*)	'shell_list =',shells
		write(100,*)	'skip_b1_tests = .true.'
		if( useBloch )		write(100,*)	'use_bloch_phases = true '
		if(.not. useBloch)	write(100,*)	'use_bloch_phases = false '
		write(100,*)	
		!
		!real lattice
		write(100,*)	'begin unit_cell_cart'
			write(100,*)	'ang'
			write(100,*)	real_lattice(1,1), ' ', real_lattice(1,2), ' ', real_lattice(1,3)	
			write(100,*)	real_lattice(2,1), ' ', real_lattice(2,2), ' ', real_lattice(2,3)	
			write(100,*)	real_lattice(3,1), ' ', real_lattice(3,2), ' ', real_lattice(3,3)		
		write(100,*)	'end unit_cell_cart'
		write(100,*)	
		!
		!atom cartesian
		write(100,*)	'begin atoms_cart'
			write(100,*)	'ang'
			do at = 1, num_atoms
				write(100,*) atom_symbols(at), atoms_cart(:,at)
			end do
		write(100,*)	'end atoms_cart'
		write(100,*)
		!
		!k points
		write(100,*)"begin kpoints"
        	do qi=1,num_kpts
        	    write(100,*)	kpt_latt(:,qi)
        	end do
        write(100,*)"end kpoints"
        write(100,*)
        !
        !PROJECTIONS
		!write(100,*)	'Begin Projections'
		!if( nWfs/ nAt == 3)	write(100,*)	'H: l=0;l=1,mr=2,3'	
		!if( nWfs/ nAt == 1)	write(100,*)	'H: l=0'
		!write(100,*)	'End Projections'
		!write(100,*)
		!
		!OUTPUT JOBS
		write(100,*)	'length_unit = Ang'
		write(100,*)	'write_u_matrices = .true.'
		write(100,*)	'write_tb = .true.'
		write(100,*)	'write_hr = .true.'
		write(100,*)	'write_XYZ = .true.'
		write(100,*)	
		!
		!
		!POST W90 PARAMETER
		write(100,*)	'kmesh = ', nKx, ' ', nKy, ' ', 1
  		write(100,*)	'geninterp = .true.'
  		write(100,*)	'geninterp_alsofirstder = .true.'
  		write(100,*)	'geninterp_single_file = .true.'
  		!

  		!PLOTTING JOBS
  		if( do_w90plot ) then
			write(100,*)	'wannier_plot = .true.'
			write(100,*)	'wvfn_formatted = .true.'
			write(100,*)	'wannier_plot_supercell =', nSCx, ' ', nSCy, ' ',1
		end if

		write(100,*)	
		close(100)
		!
		!
		return
	end subroutine



	subroutine writeW90KinterpMesh()
		!input file for post w90 module( energy and band derivative interpoltation)
		integer						:: ki
		!
		open(unit=115,file=w90_Dir//seed_name//'_geninterp.kpt',action='write',access='stream',form='formatted', status='replace')
		write(115,*)	'post w90 input file, gives info about interpolation mesh'
		write(115,*)	"frac"
		write(115,*)	nK
		do ki = 1, nK 
			write(115,*)	ki," ",kpts(1,ki)/recpLatt(1,1)," ",kpts(2,ki)/recpLatt(2,2), " ", 0.0_dp
		end do
		close(115)
		!
		return
	end subroutine





!READ HELPERS
	logical function B1condition(b_k, w_b)
		! test if
		!		sum_b{w_b * b_a * b_b}	= \delta_ab 
		! is true for all a,b
		real(dp),		intent(in)		:: 	b_k(:,:), w_b(:)
		logical							:: 	xx, xy, yy
		!
		xx	= .false.
		xy	= .false.
		yy	= .false.
		!
		if(		abs(	sum(w_b(:)*b_k(1,:)*b_k(1,:)) - 1.0_dp	)		< 0.1_dp				) 			xx = .true.
		if(		abs(	sum(w_b(:)*b_k(2,:)*b_k(2,:)) - 1.0_dp	)		< 0.1_dp				) 			yy = .true.
		if(		abs(	sum(w_b(:)*b_k(1,:)*b_k(2,:)) 			)		< 0.1_dp				)			xy = .true.
		!
		B1condition = xx .and. yy .and. xy
		!
		return
	end function




end module util_w90Interf


