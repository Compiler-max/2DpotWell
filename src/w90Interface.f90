module w90Interface	

	use mathematics,	only:	dp, PI_dp, i_dp, machineP, aUtoAngstrm, aUtoEv, myExp
	use sysPara
	use blochWf,		only:	UNKoverlap
	use projection,		only:	calcAmatANA

	implicit none

	private
	public ::					w90Interf, seed_name

	character(len=3)					:: 	seed_name
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
	subroutine w90Interf(ck, En)
		complex(dp),	intent(in)		:: ck(:,:,:)	
		real(dp),		intent(in)		:: En(:,:)
		!
		seed_name	= seedName
		call writeW90input()

		write(*,*)	"[w90Interf]: wrote w90 input file"
		call w90setup()
		write(*,*)	"[w90Interf]: done with w90 setup"
		write(*,*)	"[w90Interf]: will use num_bands= ",num_bands, " bands"
		write(*,*)	"[w90Interf]: to gen   num_wann=  ",num_wann, " wnfs"
		write(*,*)	"[w90Interf]: start preparing wannierisation"
		

		call w90prepEigVal(En)
		call w90prepMmat(ck)
		if( .not. useBloch )	then
			call w90prepAmat(ck)
			write(*,*)	"[w90Interf]: projection overlap matrix calculated"
		else
			write(*,*)	"[w90Interf]: projection disabled will use bloch phase"
		end if
		call writeW90inputPost()
		write(*,*)	"[w90Interf]: done preparing wannierization input matrices"
		call writeW90KinterpMesh()
		write(*,*)	"[w90Interf]: wrote interpolation mesh file"


		!call wannier_run(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice, &
		!					kpt_latt,num_bands,num_wann,nntot,num_atoms,atom_symbols, &
		!					atoms_cart,gamma_only,M_matrix_orig,A_matrix,eigenvalues, &
		!					U_matrix,U_matrix_opt,lwindow,wann_centres,wann_spreads, &
		!					spread)


	end subroutine







!prviate
	subroutine writeW90input()
		integer				:: stat
		integer				:: qx, qy, qi, qxl, qxr, qyl, qyr, Gxl(3), Gxr(3), Gyl(3), Gyr(3)
		!
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
				Gxl(:)	= 0.0_dp
				Gxr(:)	= 0.0_dp
				Gyl(:)	= 0.0_dp
				Gyr(:)	= 0.0_dp
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
		write(100,*)	'H: l=0;l=1,mr=2,3'
		
		!write(100,*)	'He: l=0;l=1,mr=1,2'
		write(100,*)	'End Projections'
		write(100,*)
		!
		!
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



	subroutine writeW90inputPost()
		integer				:: qi, at
		!
		!create new
		open(unit=100,file=seed_name//'.win',action='write',access='stream',form='formatted', status='replace')
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
		write(100,*)	'H: l=0;l=1,mr=2,3'	
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
		return
	end subroutine





	subroutine w90prepEigVal(En)
		!convert eigenvalues from atomic units to eV
		real(dp),		intent(in)		:: En(:,:)
		integer							:: qi, n
		real(dp),		allocatable		:: eigenvalues(:,:)
		!
		allocate(	eigenvalues( num_bands,	num_kpts	)	)
		!
		if(	size(En,1) /= num_bands	) write(*,*)"[w90prepEigVal]: warning En has wrong numbers of bands"
		if(	size(En,2) /= num_kpts 	) write(*,*)"[w90prepEigVal]: warning En has wrong numbers of kpts"
		!
		!fill new file
		open(unit=110,file=seed_name//'.eig',action='write',access='stream',form='formatted', status='replace')
		do qi = 1, num_kpts
			do n = 1, num_bands
				eigenvalues(n,qi)	= En(n,qi) * aUtoEv
				write(110,*)	n, ' ', qi, ' ', eigenvalues(n,qi)
			end do
		end do
		close(110)
		!
		!
		return
	end subroutine


	subroutine writeW90KinterpMesh()
		!input file for post w90 module( energy and band derivative interpoltation)
		integer						:: ki
		!
		open(unit=115,file=seed_name//'_geninterp.dat',action='write',access='stream',form='formatted', status='replace')
		write(115,*)
		write(115,*)	"frac"
		write(115,*)	nK
		do ki = 1, nK 
			write(115,*)	ki," ",kpts(1,ki)/recpLatt(1,1)," ",kpts(2,ki)/recpLatt(2,ki), " ", 0.0_dp
		end do
		close(115)
		!
		return
	end subroutine



!M_MATRIX ROUTINES
	subroutine w90prepMmat(ck)
		complex(dp),	intent(in)		:: ck(:,:,:) 	 !ck(			nG		,	nBands  	,	nQ	)	
		integer							:: qi, nn, n, m
		complex(dp)						:: oLap
		complex(dp),	allocatable		:: M_matrix(:,:,:,:)
		real(dp)						:: gShift(2)
		!
		allocate(	M_matrix(		num_bands	,	num_bands	,	nntot	,	num_kpts	)		)
		!
		if(	size(ck,2)/= num_bands ) write(*,*)"[w90prepEigVal]: waring ab initio exp. coefficients have wrong number of bands"
		if(	size(ck,3)/= num_kpts )  write(*,*)"[w90prepEigVal]: waring ab initio exp. coefficients have wrong numbers of kpts"
		!
		!SETUP THE MATRIX
		M_matrix	= dcmplx(0.0_dp)
		!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(qi, nn, m, n, oLap, gShift)
		do qi = 1, num_kpts
			do nn = 1, nntot
				do m = 1 , num_bands
					do n = 1, num_bands
						!calc overlap of unks
						if( nncell(3,qi,nn)/= 0 ) then
							oLap	= dcmplx(1.0_dp-5)
							write(*,*)	"[w90prepMmat]: WARNING nearest neighbours in z direction used!"
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
		!$OMP END PARALLEL DO
		!
		!WRITE TO FILE
		open(unit=120,file=seed_name//'.mmn',action='write',access='stream',form='formatted', status='replace')
		write(120,*)	'overlap matrix'
		write(120,*)	num_bands, ' ', num_kpts, ' ', nntot	
		!
		do qi = 1, num_kpts
			do nn = 1, nntot
				write(120,*)	qi,' ',nnlist(qi,nn),' ',nncell(1,qi,nn),' ',nncell(2,qi,nn),' ', nncell(3,qi,nn)
				do n = 1, num_bands
					do m = 1, num_bands
						write(120,*)	dreal(M_matrix(m,n,nn,qi)), ' ', dimag(M_matrix(m,n,nn,qi))
					end do
				end do
				
			end do
		end do
		close(120)	
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



	












!A_MATRIX ROUTINES
	subroutine w90prepAmat(ck)
		complex(dp),	intent(in)		:: ck(:,:,:) 	 !ck(			nG		,	nBands  	,	nQ	)
		integer							:: qi, n, m
		complex(dp),	allocatable		:: A_matrix(:,:,:)
		!
		allocate(	A_matrix(	num_bands	,		num_wann	,	num_kpts	)	)
		!
		if(	size(ck,2)/= num_bands )	write(*,*)"[w90prepAmat]: waring ab initio exp. coefficients have wrong number of bands"
		if(	size(ck,3)/= num_kpts )  	write(*,*)"[w90prepAmat]: waring ab initio exp. coefficients have wrong numbers of kpts"
		if(	num_wann/nAt > 3) 			write(*,*)"[w90prepAmat]: can not handle more then 3 wfs per atom at the moment"
		!
		!MATRIX SETUP
		A_matrix	= dcmplx(0.0_dp)
		do qi = 1, nQ
			call calcAmatANA(qi,ck(:,:,qi), A_matrix(:,:,qi))
		end do
		!
		!WRITE TO FILE
		open(unit=120,file=seed_name//'.amn',action='write',access='stream',form='formatted', status='replace')
		write(120,*)	'overlap matrix'
		write(120,*)	num_bands, ' ', num_kpts, ' ', num_wann
		do qi = 1, num_kpts
			do n = 1, num_wann
				do m = 1, num_bands	
					write(120,*)	m, ' ', n, ' ', ' ', qi, ' ', dreal(A_matrix(m,n,qi)), ' ', dimag(A_matrix(m,n,qi))
				end do
			end do
		end do

		return 
	end subroutine






end module
