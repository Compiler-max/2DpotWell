module w90Interface	

	use mpi
	use mathematics,	only:	dp, PI_dp, i_dp, machineP, aUtoAngstrm, aUtoEv, myExp
	use sysPara
	use planeWave,		only:	calcMmat
	use projection,		only:	calcAmatANA

	implicit none

	private
	public ::					run_w90, seed_name, U_matrix, wann_centres, wann_spreads, spread

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




!public:
	subroutine run_w90(ck, En)
		complex(dp),	intent(in)		:: 	ck(:,:,:)	
		real(dp),		intent(in)		:: 	En(:,:)
		complex(dp),	allocatable		::	M_matrix_orig(:,:,:,:), A_matrix(:,:,:), U_matrix_opt(:,:,:)
		real(dp),		allocatable		::	eigenvalues(:,:)
		logical,		allocatable		::	lwindow(:,:)
		logical							::	w90success
		integer							::	n, ierr
		!
		!PREP W90 INIT
		seed_name	= seedName
		call write_W90setup_input()
		!
		!W90 
		write(*,*)	"[w90Interf]: wrote w90 input file"
		call chdir(w90_Dir)
		call run_w90setup()
		call chdir("..")
		write(*,*)	"[w90Interf]: done with w90 setup"
		write(*,*)	"[w90Interf]: will use num_bands= ",num_bands, " bands"
		write(*,*)	"[w90Interf]: to gen   num_wann=  ",num_wann, " wnfs"
		write(*,*)	"[w90Interf]: start preparing wannierisation"
		!
		!ALLOCATE
		allocate(	M_matrix_orig(	num_bands	,	num_bands	,	nntot	,	num_kpts	)		)
		allocate(	 A_matrix( 		num_bands	,	num_wann				,	num_kpts	)		)
		allocate(	U_matrix_opt(	num_bands	,	num_wann				,	num_kpts	)		)
		allocate(	U_matrix(		num_wann	,	num_wann				,	num_kpts	)		)
		allocate(	eigenvalues(								num_bands	,	num_kpts	)		)
		allocate(	lwindow(									num_bands	,	num_kpts	)		)
		allocate(	wann_centres(	3,	num_wann)	)
		allocate(	wann_spreads(		num_wann)	)
		!
		!CALC W90 INPUT 
		call w90prepMmat(ck, M_matrix_orig)
		call w90prepAmat(ck, A_matrix)
		call w90prepEigVal(En, eigenvalues)
		
		!PREP W90
		call write_W90run_input()
		write(*,*)	"[w90Interf]: done preparing wannierization input matrices"
		!PREP PW90
		call writeW90KinterpMesh()
		write(*,*)	"[w90Interf]: wrote interpolation mesh file"
	
		!WANNIERISE
		!call wannier_run(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice, &
		!					kpt_latt,num_bands,num_wann,nntot,num_atoms,atom_symbols, &
		!					atoms_cart,gamma_only,M_matrix_orig,A_matrix,eigenvalues, &
		!					U_matrix,U_matrix_opt,lwindow,wann_centres,wann_spreads, &
		!					spread)

		!note: the library modus gives different results compared to executing w90 manually, why is currently unclear

		inquire(file="../thirdparty/wannier90/wannier90.x ", exist=w90success)
		if( w90success ) then
			call chdir(w90_Dir)
			call execute_command_line ("../../thirdparty/wannier90/wannier90.x "//seed_name, exitstat=ierr)
 	 		print *, "Exit status of wannier90 was ", ierr
			call chdir("..")
			write(*,*)	"[w90Interf]: wannierization finished"
		else
			write(*,*)	"[w90Interf]: please provide wannier90.x in /thirdpary/wannier90/ or execute wannier90 manually"
		end if
		!
		!
		end do
	end subroutine
















!prviate
	subroutine write_W90setup_input()
		integer				:: stat
		integer				:: qx, qy, qi, qxl, qxr, qyl, qyr, Gxl(3), Gxr(3), Gyl(3), Gyr(3)
		!
		seed_name			= 'wf1'
		!
		!Delete old input file
		open(unit=200, iostat=stat, file=w90_Dir//seed_name//'.win', status='old')
		if (stat == 0) close(200, status='delete')
		!Delete old wout file
		open(unit=210, iostat=stat, file=w90_Dir//seed_name//'.wout', status='old')
		if (stat == 0) close(210, status='delete')
		!		
		!
		!
		!create new
		open(unit=100,file=w90_Dir//seed_name//'.win',action='write', status='new')
		!
		!BASIC INFO
		write(100,*)	'num_wann  = ',nWfs
		write(100,*)	'num_bands = ',nBands
		write(100,*)	'mp_grid   = ', mp_grid(1) , ' ', mp_grid(2), ' ', mp_grid(3)
		if( useBloch )	write(100,*)	'use_bloch_phases = true '
		write(100,*)	
		!


		!write(100,*)	'shell_list = ',1
		!
		!k points
		!write(100,*)"begin kpoints"
        !do qi=1,nQ
        !    write(100,*)qpts(:,qi),0.0_dp
        !end do
        !write(100,*)"end kpoints"
        !write(100,*)

		! k point neighbours
		write(100,*)	'postproc_setup = .true.'
		write(100,*)	'begin nnkpts'
		!get neighbours and write to file
		do qx = 1, nQx
			do qy = 1, nQy
				qxl	= getLeft( qx,nQx)
				qxr	= getRight(qx,nQx)
				qyl	= getLeft(  qy,nQy)
				qyr = getRight( qy,nQy)
				!
				!GET GRID POSITION OF NEIGHBOURS
				qi	= getkindex(qx,qy)
				qxl	= getkindex(qxl,qy)
				qxr	= getkindex(qxr,qy)
				qyl	= getkindex(qx,qyl)
				qyr	= getkindex(qx,qyr)
				!
				!SHIFT NEIGHBOURS BACK TO FIRST BZ
				Gxl(:)	= 0
				Gxr(:)	= 0
				Gyl(:)	= 0
				Gyr(:)	= 0
				if( qx == 1 ) 	Gxl(1)	= - 1
				if( qx == nQx)	Gxr(1)	= + 1
				if( qy == 1 ) 	Gyl(2)	= - 1
				if( qy == nQy)	Gyr(2)	= + 1
				!
				write(100,*)	qi,' ',qxl,' ',Gxl(1),' ',Gxl(2),' ', Gxl(3)
				write(100,*)	qi,' ',qxr,' ',Gxr(1),' ',Gxr(2),' ', Gxr(3)
				write(100,*)	qi,' ',qyl,' ',Gyl(1),' ',Gyl(2),' ', Gyl(3)
				write(100,*)	qi,' ',qyr,' ',Gyr(1),' ',Gyr(2),' ', Gyr(3)
				!
			end do 
		end do
		write(100,*)	'end nnkpts'
		write(100,*)
		!
		!PROJECTIONS
		write(100,*)	'Begin Projections'
		!write(100,*)	'H: l=',0,';l=',1,',mr=',1,';l=',1,',mr=',2
		if(nWfs == 6) write(100,*)	'H: l=0;l=1,mr=2,3'
		if(nWfs == 2) write(100,*)	'H: l=0'
		
		!write(100,*)	'He: l=0;l=1,mr=1,2'
		write(100,*)	'End Projections'
		write(100,*)
		!
		!
		close(100)
		return
	end subroutine




	subroutine run_w90setup()
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
		real_lattice(3,3)	= aX * aY * aUtoAngstrm
		!reciprocal cell
		recip_lattice		= 0.0_dp
		recip_lattice(1,1)	= 2.0_dp * PI_dp / real_lattice(1,1) 	
		recip_lattice(2,2)	= 2.0_dp * PI_dp / real_lattice(2,2)
		recip_lattice(3,3)	= 2.0_dp * PI_dp / real_lattice(3,3) 		
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
			atoms_cart(1:2,at)	= atPos(1:2,at)*aUtoAngstrm
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
		!			
		return
	end subroutine



	subroutine write_W90run_input()
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
		write(100,*)	'shell_list =',shell
		write(100,*)	'skip_b1_tests = .true.'
		if( useBloch )	write(100,*)	'use_bloch_phases = true '
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
		!atom cart
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
		if(nWfs == 6)	write(100,*)	'H: l=0;l=1,mr=2,3'	
		if(nWfs == 2)	write(100,*)	'H: l=0'
		write(100,*)	'End Projections'
		write(100,*)
		!
		!JOBS for TB basis
		write(100,*)	'write_u_matrices = .true.'
		write(100,*)	'write_tb = .true.'
		write(100,*)	'write_hr = .true.'
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


	subroutine w90prepEigVal(au_units, eV_units)
		!convert eigenvalues from atomic units to eV
		real(dp),		intent(in)		:: 	au_units(:,:)
		real(dp),		intent(out)		::	eV_units(:,:)
		integer							:: 	qi, n		!
		!
		!WRITE FILE
		open(unit=110,file=w90_Dir//seed_name//'.eig',action='write',access='stream',form='formatted', status='replace')
		do qi = 1, num_kpts
			do n = 1, num_bands
				eV_units(n,qi)	= aUtoEv * au_units(n,qi)
				write(110,*)	n, ' ', qi, ' ', eV_units(n,qi)
			end do
		end do
		close(110)
		write(*,'(a,i3,a)')	"[#",myID,";w90prepEigVal]: wrote .eig file"
		!
		!DEBUG
		if(	size(au_units,2) /= nQ	)	write(*,'(a,i3,a)')		"[#",myID,";w90prepEigVal]: warning En has wrong numbers of kpts"
		!
		!
		return
	end subroutine


!M_MATRIX ROUTINES
	subroutine w90prepMmat(ck, M_mat)
		complex(dp),	intent(in)		::	ck(:,:,:) 	 !ck(			nG		,	nBands  	,	nQ	)	
		complex(dp),	intent(out)		::	M_mat(:,:,:,:)
		integer							::	qi, nn, n, m
		real(dp)						::	gShift(2)							
		!
		!FILL THE MATRIX
		M_mat	= dcmplx(0.0_dp)
		write(*,*)	"[w90prepMmat]: use ",nntot," nearest neighbours per k point"
		!!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(qi, nn, gShift)
		do qi = 1, nQ
			do nn = 1, nntot
				!calc overlap of unks
				if( nncell(3,qi,nn)/= 0 ) then
					M_mat(:,:,nn,qi)	= dcmplx(1.0_dp-5)
					write(*,*)	"[w90prepMmat]: WARNING nearest neighbours in z direction used!"
					if(qi==1  ) write(*,*)	"[w90prepMmat]: oLap set to zero for nn=",nn
				else
					gShift(1)			= nncell(1,qi,nn) * 2.0_dp * PI_dp / aX
					gShift(2)			= nncell(2,qi,nn) * 2.0_dp * PI_dp / aX
					!oLap				= UNKoverlap(n,m, qi, nnlist(qi,nn), gShift, ck)
					!call calcMmat(qi,nnlist(qi,nn), gShift, ck, M_loc(:,:,nn,qi))
					call calcMmat(qi, nnlist(qi,nn), gShift, nGq, Gvec, ck, M_mat(:,:,nn,qi))
				end if
			end do
		end do
		!!$OMP END PARALLEL DO	
		!
		!WRITE TO FILE
		open(unit=120,file=w90_Dir//seed_name//'.mmn',action='write',access='stream',form='formatted', status='replace')
		write(120,*)	'overlap matrix'
		write(120,*)	num_bands, ' ', num_kpts, ' ', nntot	
		!
		do qi = 1, num_kpts
			do nn = 1, nntot
				write(120,*)	qi,' ',nnlist(qi,nn),' ',nncell(1,qi,nn),' ',nncell(2,qi,nn),' ', nncell(3,qi,nn)
				do n = 1, num_bands
					do m = 1, num_bands
						write(120,*)	dreal(M_mat(m,n,nn,qi)), ' ', dimag(M_mat(m,n,nn,qi))
					end do
				end do		
			end do
		end do
		close(120)
		write(*,'(a,i3,a)')	"[#",myID,";w90prepMmat]: wrote .mmn file"
		!
		!
		return
	end subroutine

!A_MATRIX ROUTINES
	subroutine w90prepAmat(ck, A_mat)
		complex(dp),	intent(in)		:: ck(:,:,:) 	 !ck(			nG		,	nBands  	,	nQ	)
		complex(dp),	intent(out)		:: A_mat(:,:,:)
		integer							:: qi, n, m
		!
		if( .not. useBloch )	then
			!
			if(	size(ck,2) < num_bands )	write(*,'(a,i3,a)')"[#",myID,";w90prepAmat]: warning not enough abInitio coeff, try to increase nSolve in inpupt"
			if(	size(ck,3)/= qChunk )  	write(*,'(a,i3,a)')"[#",myID,";w90prepAmat]: warning ab initio exp. coefficients have wrong numbers of kpts"
			if(	num_wann/nAt > 3) 			write(*,'(a,i3,a)')"[#",myID,";w90prepAmat]: can not handle more then 3 wfs per atom at the moment"
			!
			!MATRIX SETUP
			!$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(qi)		
			do qi = 1, nQ
				call calcAmatANA(qi,ck(:,:,qi), A_mat(:,:,qi))
			end do
			!$OMP END PARALLEL DO
			!
			!WRITE TO FILE
			open(unit=120,file=w90_Dir//seed_name//'.amn',action='write',access='stream',form='formatted', status='replace')
			write(120,*)	'overlap matrix'
			write(120,*)	num_bands, ' ', num_kpts, ' ', num_wann
			do qi = 1, num_kpts
				do n = 1, num_wann
					do m = 1, num_bands	
						write(120,*)	m, ' ', n, ' ', ' ', qi, ' ', dreal(A_mat(m,n,qi)), ' ', dimag(A_mat(m,n,qi))
					end do
				end do
			end do
			write(*,'(a,i3,a)')	"[#",myID,";w90prepAmat]: wrote .amn file"
			!
		else
			if( myID==root )	write(*,*)	"[w90Interf]: projection disabled will use bloch phase"
		end if
		!
		!
		return 
	end subroutine










	integer function getLeft(i,N)
		!HELPER for calcConn
		!gets left (lower) neighbour, using the periodicity at boundary
		!
		integer,	intent(in)	:: i,N
		if(i.eq.1) then
			getLeft = N
		else
			getLeft = i-1
		end if
		!
		return
	end function


	integer function getRight(i,N)
		!HELPER for calcConn
		!gets right (upper) neighbour, using the periodicity at boundary
		!
		integer,	intent(in)	:: i,N
		if(i.eq.N) then
			getRight = 1
		else
			getRight = i+1
		end if
		!
		return
	end function





end module





!ALL GATHER M_MATRIX: 

!!!!M_MATRIX ROUTINES
!!!	subroutine w90prepMmat(ck_loc)
!!!		complex(dp),	intent(in)		:: ck_loc(:,:,:) 	 !ck(			nG		,	nBands  	,	nQ	)	
!!!		integer							:: qi, nn, n, m, mesgSize
!!!		integer,		allocatable		:: nGq_glob(:)
!!!		complex(dp),	allocatable		:: M_glob(:,:,:,:), M_loc(:,:,:,:), ck_glob(:,:,:)
!!!		real(dp),		allocatable		:: Gvec_glob(:,:,:), G_send(:,:,:)
!!!		real(dp)						:: gShift(2)
!!!		!
!!!							allocate(	nGq_glob(			nQ													)		)
!!!							allocate(	ck_glob(		GmaxGLOBAL	,	nSolve		,			nQ				)		)
!!!							allocate(	Gvec_glob(		dim			,	GmaxGLOBAL	, 			nQ				)		)
!!!							allocate(	G_send(			dim			,	GmaxGLOBAL	,			nQ				)		)
!!!							allocate(	M_loc(			num_bands	,	num_bands	,	nntot	,	qChunk		)		)
!!!		if(myID == root)	allocate(	M_glob(			num_bands	,	num_bands	,	nntot	,	num_kpts	)		)
!!!		if(myID /= root)	allocate(	M_glob(				0		,		0		,		0	,		0		)		)
!!!							
!!!		!
!!!		!GATHER _glob ARRAYS 
!!!		call MPI_ALLGATHER(	nGq							, 	qChunk					, 	MPI_INTEGER			, &
!!!							nGq_glob					, 	qChunk					, 	MPI_INTEGER			, 	MPI_COMM_WORLD, ierr)
!!!		if(myID==root) write(*,'(a,i3,a)')	"[#",myID,";w90prepMmat: bcasted nGq"
!!!
!!!		
!!!			G_send(1:dim, 1:GmaxGLOBAL,1:qChunk) = Gvec(1:dim, 1:GmaxGLOBAL, 1:qChunk)
!!!		call MPI_ALLGATHER(	G_send						,	dim*GmaxGLOBAL*qChunk	,	MPI_DOUBLE_PRECISION, &
!!!							Gvec_glob					,  	dim*GmaxGLOBAL*qChunk	,	MPI_DOUBLE_PRECISION, 	MPI_COMM_WORLD, ierr)	
!!!		if(myID==root) write(*,'(a,i3,a)')	"[#",myID,";w90prepMmat: bcasted Gvec"
!!!
!!!
!!!		call MPI_ALLGATHER(	ck_loc						,	GmaxGLOBAL*nSolve*qChunk, 	MPI_DOUBLE_COMPLEX	, &	
!!!							ck_glob(1:GmaxGLOBAL, :,:)	,  	GmaxGLOBAL*nSolve*qChunk, 	MPI_DOUBLE_COMPLEX	, 	MPI_COMM_WORLD, ierr)	
!!!		if(myID==root) write(*,'(a,i3,a)')	"[#",myID,";w90prepMmat: bcasted basis coeff"
!!!	
!!!
!!!
!!!		if(myID==root) write(*,'(a,i3,a)')	"[#",myID,";w90prepMmat: basis broadcasted"
!!!		!
!!!		!FILL THE MATRIX
!!!		M_loc	= dcmplx(0.0_dp)
!!!		write(*,*)	"[w90prepMmat]: use ",nntot," nearest neighbours per k point"
!!!		do qi = myID*qChunk+1, myID*qChunk+qChunk
!!!			do nn = 1, nntot
!!!				!calc overlap of unks
!!!				if( nncell(3,qi,nn)/= 0 ) then
!!!					M_loc(:,:,nn,qi)	= dcmplx(1.0_dp-5)
!!!					write(*,*)	"[w90prepMmat]: WARNING nearest neighbours in z direction used!"
!!!					if(qi==1  ) write(*,*)	"[w90prepMmat]: oLap set to zero for nn=",nn
!!!				else
!!!					gShift(1)			= nncell(1,qi,nn) * 2.0_dp * PI_dp / aX
!!!					gShift(2)			= nncell(2,qi,nn) * 2.0_dp * PI_dp / aX
!!!					!oLap				= UNKoverlap(n,m, qi, nnlist(qi,nn), gShift, ck)
!!!					!call calcMmat(qi,nnlist(qi,nn), gShift, ck, M_loc(:,:,nn,qi))
!!!					call calcMmat(qi, nnlist(qi,nn), gShift, nGq_glob, Gvec_glob, ck_glob, M_loc(:,:,nn,qi))
!!!				end if
!!!			end do
!!!		end do
!!!		!
!!!		!COLLECT
!!!		mesgSize = num_bands**2 * nntot * qChunk
!!!		call MPI_GATHER(M_loc, mesgSize, MPI_DOUBLE_COMPLEX, M_glob, mesgSize, MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD, ierr )
!!!		if(myID==root) write(*,'(a,i3,a)')	"[#",myID,";w90prepMmat]: overlap matrix collected"
!!!		!
!!!		!
!!!		!WRITE TO FILE
!!!		if( myID == root ) then
!!!			!W90 input
!!!			open(unit=120,file=w90_Dir//seed_name//'.mmn',action='write',access='stream',form='formatted', status='replace')
!!!			write(120,*)	'overlap matrix'
!!!			write(120,*)	num_bands, ' ', num_kpts, ' ', nntot	
!!!			!
!!!			do qi = 1, num_kpts
!!!				do nn = 1, nntot
!!!					write(120,*)	qi,' ',nnlist(qi,nn),' ',nncell(1,qi,nn),' ',nncell(2,qi,nn),' ', nncell(3,qi,nn)
!!!					do n = 1, num_bands
!!!						do m = 1, num_bands
!!!							write(120,*)	dreal(M_glob(m,n,nn,qi)), ' ', dimag(M_glob(m,n,nn,qi))
!!!						end do
!!!					end do		
!!!				end do
!!!			end do
!!!			close(120)
!!!			write(*,'(a,i3,a)')	"[#",myID,";w90prepMmat]: wrote .mmn file"
!!!			
!!!			!POST W90 input
!!!			call writeABiN_basCoeff(ck_glob)
!!!			write(*,'(a,i3,a)')	"[#",myID,";w90prepMmat]: wrote  basis coeff"
!!!			call writeABiN_basis(nGq_glob, Gvec_glob)
!!!			write(*,'(a,i3,a)')	"[#",myID,";w90prepMmat]: wrote  basis"
!!!		end if	
!!!		!
!!!		!
!!!		return
!!!	end subroutine