module w90Interface	

	use mpi
	use mathematics,	only:	dp, PI_dp, i_dp, machineP, aUtoAngstrm, aUtoEv, myExp, isUnit, isIdentity
	use sysPara
	use basisIO,		only:	read_energies, read_Amn, read_Mmn
	use planeWave,		only:	calcMmat, calcAmatANA

	implicit none

	private
	public ::					setup_w90, write_w90_matrices, &
								read_FD_scheme, read_M_initial, read_U_matrix, readBandVelo, read_wann_centers, & 
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

		!WRITE INPUT FILE (wann setup)
		call write_W90setup_input()
		!
		!W90 
		write(*,*)	"[w90Interf]: wrote w90 setup input file (.win)"
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
	subroutine read_FD_scheme(nntot, nnlist, nncell, b_k, w_b)
		!call after w90 is finished
		!reads seedname.nnkp & seedname.wout
		integer, 					intent(out)		::	nntot
		integer,	allocatable, 	intent(out)		::	nnlist(:,:), nncell(:,:,:)
		real(dp),	allocatable,	intent(out)		::	b_k(:,:), w_b(:)
		integer										::	stat, qi, qnn, iLine(5), nn, start, end
		real(dp)									::	b_vec(3), weight
		character(len=*), 			parameter 		::	search_nnkp = "begin nnkpts", final_nnkp = "end nnkpts" 
		character(len=*), 			parameter 		::	search_wout=" |                  b_k Vectors (Ang^-1) and Weights (Ang^2)                  |"

		character(len=100)							::	line
		logical										::	finished, file_exists
		!
		!inquire file
		inquire(file=w90_dir//seedName//'.nnkp',exist=file_exists)
		if( .not. file_exists	) stop '[read_FD_scheme]: could not find .nnkp file'
		!
		!read file
		open(unit=330, iostat=stat, file=w90_dir//seedName//'.nnkp', form='formatted', status='old', action='read')
		finished = .false.
		do while( .not. finished .and. stat==0)
			read(330,"(a)",iostat=stat) line
			!find the block
			if(line == search_nnkp) then
				!GET nntot				
				read(330,*)	nntot
				allocate(	nnlist(		nQ, nntot)		)
				allocate(	nncell(	3, 	nQ, nntot)		)
				allocate(	b_k(	3,		nntot)		)
				allocate(	w_b(			nntot)		)
				!read data
				do qi = 1, nQ
					do qnn = 1, nntot
						read(330, *) iLine
						nnlist( iLine(1), qnn )		= iLine(2)
						nncell(	1:3, iLine(1), qnn)	= iLine(3:5)
					end do
				end do
				finished = .true.
			end if
		end do 
		close(330)
		if( finished ) then
			write(*,*)	"[readFDscheme]: read the .nnkp file"
			write(*,*)	"[readFDscheme]: nntot=",nntot
		else
			stop	"[readFDscheme]: did not find search string  in .nnkp file"
		end if
		!
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
				do qnn = 1, nntot
					!read line
					read(335,"(a)",iostat=stat) line
					!remove "|" from line  (this is neccessary to use pattern matching on string)
					start 	= scan(line,"|" ) 			+1
					end 	= scan(line,"|", .true.)	-1
					!extract values from substring
					read(line(start:end),*)	nn, b_vec, weight 
					!convert to a.u.
					b_k(1:3, nn) 	= b_vec(1:3)*aUtoAngstrm			![b_k] = (Ang^-1)
					w_b(nn)			= weight / aUtoAngstrm**2			![w_b] = (Ang^2)
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
		write(*,*)		" 	nn | 	b	| 	w_b "
		do nn = 1, nntot
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
		w_centers 	= w_centers 	/ aUtoAngstrm
		write(*,*)	"[read_wann_centers]: read the .xyz file"
		!
		!
		return
	end subroutine



	subroutine read_U_matrix(U_matrix)
		!READ U FROM W90
		complex(dp),	intent(out)	:: 	U_matrix(:,:,:)
		real(dp),		allocatable	:: 	krel(:,:)
		integer						:: 	stat, qi, n, m, dumI(3)
		real(dp)					:: 	val(2)
		logical						::	file_exists
		!
		U_matrix	= dcmplx(0.0_dp)

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
		if( dumI(1) /= size(U_matrix,3) ) 	stop	"[readUmatrix]:  ERROR - wrong qmesh size found in file (doesn't match allo size)"
		if(	dumI(2)	/= size(U_matrix,2)) 	stop	"[readUmatrix]:  ERROR - wrong number of wannier functions found in file (doesn't match allo size)"
		allocate(	krel( 3, dumI(1))	)
		!
		!body
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
					write(*,*)	"				x: k_w90= ",krel(1,qi)*2.0_dp*PI_dp/aX, " qpts=",qpts(1,qi)
					write(*,*)	"				y: k_w90= ",krel(2,qi)*2.0_dp*PI_dp/aY, " qpts=",qpts(2,qi)
					stop	"[readUmatrix]: ERROR k meshes are ordered diffferently"
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


	subroutine read_M_initial( M_init)
		!read the M_matrix written by w90prepMmat	
		complex(dp),		intent(out)		::	M_init(:,:,:,:)
		integer								::	qi,nn, n, m, &
												f_num_bands, f_num_kpts, f_nntot, f_qi, f_nnlist, f_nncell(1:3)
		real(dp)							::	realBuff(2)
		logical								:: 	foundFile
		!
		M_init = dcmplx(0.0_dp)
		!check for file
		inquire(file=w90_Dir//seedName//'.mmn', exist=foundFile)
		if( .not. foundFile ) stop 'the .mmn file could not be found'
		!
		!read file
		open(unit=120,file=w90_Dir//seedName//'.mmn',action='read',access='stream',form='formatted', status='old')
		read(120,*)
		read(120,*) f_num_bands, f_num_kpts, f_nntot		
		!
		do qi = 1, size(M_init,4)
			do nn = 1, size(M_init,3)
				!
				read(120,*)	f_qi, f_nnlist, f_nncell(1:3)
				if(	 qi /= 	f_qi )	stop	"[read_M_initial]: WARNING q mesh ordered differently"
				!
				do n = 1, size(M_init,2)
					do m = 1, size(M_init,1)
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
		open(unit=125,file=w90_Dir//'clone'//'.mmn',action='write',access='stream',form='formatted', status='replace')
		write(125,*)	"clone of the "//seedName//".mmn file"
		write(125,*)	f_num_bands, f_num_kpts, f_nntot
		do qi = 1, size(M_init,4)
			do nn = 1, size(M_init,3)
				!
				write(125,*)	qi, nn
				!
				do n = 1, size(M_init,2)
					do m = 1, size(M_init,1)
						write(125,*)	M_init(m,n,nn,qi)
					end do
				end do
				!
			end do
		end do
		close(125)

		!
		return
	end subroutine


	subroutine readBandVelo( v_vec )
		real(dp),		intent(out)		::	v_vec(:,:,:)
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
		do qi = 1, nQ
			do n = 1, nWfs
				read(320,*)	qInd, buffer
				v_vec(1,n,qInd)	= buffer(5)
				v_vec(2,n,qInd)	= buffer(6)
				v_vec(3,n,qInd)	= buffer(7) 
			end do
		end do
		close(320)
		!
		!ATOMIC UNITS CONVERSION:
		v_vec 	= v_vec  / (aUtoEv  * aUtoAngstrm)	 ! [v_vec] = eV / Angstroem
		write(*,*)	"[readBandVelo]: read the geninterp.dat file"
		!
		!	
		return
	end subroutine








































!prviate
!PREPARE & RUN W90
	subroutine write_W90setup_input()
		!write input file for wannier_setup call
		integer	:: i
		!
		open(unit=100,file=seed_name//'.win',action='write', status='replace')
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
		write(100,*)	'shell_list = ', (shells(i), i=1,size(shells) )
		write(100,*)	'skip_b1_tests = .true.'
		!

		! k point neighbours
		write(100,*)	'postproc_setup = .true.'
		write(100,*)
		!
		!PROJECTIONS
		write(100,*)	'Begin Projections'
		if(nWfs / nAt == 3) write(100,*)	'H: l=0;l=1,mr=2,3'
		if(nWfs / nAt == 1) write(100,*)	'H: l=0'
		write(100,*)	'End Projections'
		write(100,*)
		!
		!
		close(100)
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
			atom_symbols(at)	= 'H'
			atoms_cart(1:2,at)	= atPos(1:2,at)*aUtoAngstrm
			atoms_cart(3,at)	= 0.0_dp
		end do
		!kpt lattice (fractional co-ordinates relative to the reciprocal lattice vectors)
		kpt_latt(1,:)			= qpts(1,:) / ( recip_lattice(1,1) * aUtoAngstrm ) 
		kpt_latt(2,:)			= qpts(2,:)	/ ( recip_lattice(2,2) * aUtoAngstrm ) 
		kpt_latt(3,:)			= 0.0_dp
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
		write(100,*)	'Begin Projections'
		if( nWfs/ nAt == 3)	write(100,*)	'H: l=0;l=1,mr=2,3'	
		if( nWfs/ nAt == 1)	write(100,*)	'H: l=0'
		write(100,*)	'End Projections'
		write(100,*)
		!
		!OUTPUT JOBS
		write(100,*)	'write_u_matrices = .true.'
		write(100,*)	'write_tb = .true.'
		write(100,*)	'write_hr = .true.'
		write(100,*)	'write_XYZ = .true.'
		write(100,*)	
		!
		!PLOTTING JOBS
		!write(100,*)	'wannier_plot = .true.'
		!write(100,*)	'wannier_plot_supercell =', nSCx, ' ', nSCy, ' ',1
		!
		!POST W90 PARAMETER
		write(100,*)	'kmesh = ', nKx, ' ', nKy, ' ', 1
  		write(100,*)	'geninterp = .true.'
  		write(100,*)	'geninterp_alsofirstder = .true.'
  		write(100,*)	'geninterp_single_file = .true.'
  		!
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




end module


