module basisIO
	use sysPara
	use mpi
	use omp_lib
	use mathematics,	only:	dp, i_dp, PI_dp, machineP, aUtoEv, aUtoAngstrm

	implicit none

	private
	public :: 		writeABiN_basVect, writeABiN_energy, writeABiN_basis, writeABiN_basCoeff, writeABiN_velo, writeABiN_Amn,  writeABiN_Mmn,&
					read_coeff, read_gVec, read_energies, read_velo, readBasis, read_Amn, read_Mmn


	character(len=1024)				::	format='(a,i7.7)'


	contains





!WRITE
	subroutine writeABiN_basVect(qi, Gvec)
		integer,		intent(in)		::	qi
		real(dp),		intent(in)		::	Gvec(:,:)
		character(len=20)				::	filename
		
		!REAL
		write(filename, format) raw_dir//'gVec.',qi
		open(unit=210, file=filename		, form='unformatted', access='stream', action='write',status='replace') 
		write(210)	Gvec
		close(210)
		!
		return
	end subroutine



	subroutine writeABiN_energy(qi, En_loc)
		integer,		intent(in)		::	qi
		real(dp),		intent(in)		:: 	En_loc(:)
		character(len=20)				::	filename
		!
		!WRITE TO FILE
		write(filename, format) raw_dir//'enK.',qi
		open(unit=200, file=filename, form='unformatted', access='stream', action='write', status='unknown')
		write(200)	En_loc
		close(200)
		!
		!
		return
	end subroutine


	subroutine writeABiN_basCoeff(qi, ck)
		integer,		intent(in)		::	qi
		complex(dp),	intent(in)		::	ck(:,:)
		character(len=20)				::	filename
		!
		!REAL
		write(filename, format) raw_dir//'ckR.',qi
		open(unit=210, file=filename		, form='unformatted', access='stream', action='write',status='replace') 
		write(210)	dreal(ck)
		close(210)
		!
		!IMAG
		write(filename, format) raw_dir//'ckI.',qi
		open(unit=215, file=filename		, form='unformatted', access='stream', action='write',status='replace') 
		write(215)	dimag(ck)
		close(215)
		!
		return
	end subroutine


	subroutine writeABiN_Amn(qi, Amn)
		integer,		intent(in)		::	qi
		complex(dp),	intent(in)		::	Amn(:,:)
		character(len=20)				::	filename		
		!
		!REAL
		write(filename, format) raw_dir//'AmnR.',qi
		open(unit=210, file=filename		, form='unformatted', access='stream', action='write',status='replace') 
		write(210)	dreal(Amn)
		close(210)
		!
		!IMAG
		write(filename, format) raw_dir//'AmnI.',qi
		open(unit=215, file=filename		, form='unformatted', access='stream', action='write',status='replace') 
		write(215)	dimag(Amn)
		close(215)
		!
		return
	end subroutine


	subroutine writeABiN_Mmn(qi, Mmn)
		integer,		intent(in)		::	qi
		complex(dp),	intent(in)		::	Mmn(:,:,:)
		character(len=20)				::	filename		
		!
		!
		!REAL
		write(filename, format) raw_dir//'MmnR.',qi
		open(unit=210, file=filename		, form='unformatted', access='stream', action='write',status='replace') 
		write(210)	dreal(Mmn)
		close(210)
		!
		!IMAG
		write(filename, format) raw_dir//'MmnI.',qi
		open(unit=215, file=filename		, form='unformatted', access='stream', action='write',status='replace') 
		write(215)	dimag(Mmn)
		close(215)
		!
		return
	end subroutine



	subroutine writeABiN_velo(qi, velo)
		integer,		intent(in)		::	qi
		complex(dp),	intent(in)		::	velo(:,:,:)
		character(len=20)				::	filename		
		!
		!
		!REAL
		write(filename, format) raw_dir//'velR.',qi
		open(unit=210, file=filename		, form='unformatted', access='stream', action='write',status='replace') 
		write(210)	dreal(velo)
		close(210)
		!
		!IMAG
		write(filename, format) raw_dir//'velI.',qi
		open(unit=215, file=filename		, form='unformatted', access='stream', action='write',status='replace') 
		write(215)	dimag(velo)
		close(215)
		!
		return
	end subroutine





	subroutine writeABiN_basis(nGq_loc, Gvec_loc)
		integer,		intent(in)		::	nGq_loc(:)
		real(dp),		intent(in)		::	Gvec_loc(:,:,:)
		integer,		allocatable		::	nGq_glob(:)
		real(dp),		allocatable		::	Gvec_glob(:,:,:)

		integer							::	qi
		!
		!ALLOCATE TARGET
		if( myID == root ) then
			allocate(	nGq_glob(	nQ	) 							)	
			allocate( 	Gvec_glob( size(Gvec,1), size(Gvec,2), nQ)	)			
		else
			allocate(	nGq_glob( 0 	)							)			
			allocate(	Gvec_glob( 0, 0, 0)							)
		end if
		!
		!GATHER
		call MPI_GATHER( nGq_loc	, qChunk, MPI_INTEGER, nGq_glob		, qChunk, MPI_INTEGER, MPI_COMM_WORLD, ierr)	
		call MPI_GATHER( Gvec_loc	, qChunk, MPI_DOUBLE_PRECISION, Gvec_glob	, qChunk, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)	
		!
		!WRITE TO FILE
		if( myID == root ) then
			!NGQ
			open(unit=215, file=raw_dir//'nGq.dat'		, form='unformatted', access='stream', action='write', status='replace')
			write(215)	nGq_glob
			close(215)
			!REAL GVEC
			open(unit=220, file=raw_dir//'Gvec.dat'		, form='unformatted', access='stream', action='write',status='replace') 
			do qi = 1, size(Gvec_glob,3)
				write(220)	Gvec_glob(:,:,qi)
			end do
			close(220)
			write(*,*)	"[writeABiN_basis]: wrote nGq and Gvec to binary files"
		end if
		!
		return
	end subroutine








!READ
	subroutine	read_coeff(qi,ck)
		integer,		intent(in)		:: 	qi
		complex(dp),	intent(out)		:: 	ck(:,:)
		real(dp),		allocatable		:: 	buffer(:,:)
		character(len=20)				:: 	filename
		logical							:: 	file_exists
		!
		allocate(	buffer( size(ck,1), size(ck,2) 	)		)
		!
		!
		!UNK REAL PART
		write(filename, format) raw_dir//'ckR.',qi
		inquire(file=filename,exist=file_exists)
		if(.not. file_exists)	then
			write(*,'(a,i3,a,i5,a)')	"[#",myID,";read_coeff]: could not find ckR.",qi," file in rawData"
			stop 'could not find basis coefficients (ckR files)'
		end if
		open(unit=700, file=filename ,form='unformatted',access='stream',action='read')
		read(700) buffer
		ck(:,:)	= dcmplx(buffer)
		close(700)
		!
		!UNK IMAG PART
		write(filename, format) raw_dir//'ckI.',qi
		inquire(file=filename,exist=file_exists)
		if(.not. file_exists)	then
			write(*,'(a,i3,a,i5,a)')	"[#",myID,";read_coeff]: could not find ckI.",qi," file in rawData"
			stop 'could not find basis coefficients (ckI files)'
		end if
		open(unit=710, file=filename,form='unformatted',access='stream',action='read')
		read(710) buffer
		ck(:,:)	= ck(:,:) + i_dp * dcmplx(buffer)
		close(710)
		!
		!
		return
	end subroutine


	subroutine read_gVec(qi, Gvec)
		integer,		intent(in)		::	qi
		real(dp),		intent(out)		::	Gvec(:,:)
		character(len=20)				::	filename
		logical							:: 	file_exists
		!
		!real part
		write(filename, format) raw_dir//'gVec.',qi
		!
		inquire(file=filename,exist=file_exists)
		if(.not. file_exists)	stop 'could not find basis coefficients (gVec files)'
		!
		open(unit=700, file=filename ,form='unformatted',access='stream',action='read')
		read(700) Gvec
		close(700)
		!
		!
		return
	end subroutine


	subroutine read_energies(En)
		real(dp),		intent(out)		::	En(:,:)
		integer							::	qi
		character(len=20)				::	filename
		logical							::	file_exists
		!
		!BAND ENERGIES
		do qi = 1, size(En,2)	
				write(filename, format) raw_dir//'enK.',qi
				inquire(file=filename,exist=file_exists)
				if(.not. file_exists)	stop 'could not find all eigenvalues (enK files)'
				open(unit=720, file=filename, form='unformatted', access='stream', action='read')
				read(720)	En(:,qi)
				close(720)
		end do
		write(*,*)	"[read_energies]: read the enK.qi files (abinit files)"
		!
		return
	end subroutine

	subroutine readBasis()
		!
		open(unit=730, file=raw_dir//"nGq.dat" ,form='unformatted',access='stream',action='read')
		read(730)	nGq
		close(730)
		!
		open(unit=740, file=raw_dir//"Gvec.dat", form='unformatted', access='stream', action='read')
		read(740)	Gvec
		close(740)
		!
		!
		return
	end subroutine


	subroutine read_Amn(Amn)
		complex(dp),		intent(out)		::	Amn(:,:,:)
		integer								::	qi
		character(len=20)					::	filename
		real(dp),		allocatable			::	buffer(:,:)
		!
		allocate(	buffer( size(Amn,1), size(Amn,2) )		)
		!
		do qi = 1, size(Amn,3)
			!REAL
			write(filename, format) raw_dir//'AmnR.',qi
			open(unit=210, file=filename		, form='unformatted', access='stream', action='read',status='old') 
			read(210)	buffer
			close(210)
			Amn(:,:,qi)	= dcmplx(buffer(:,:))
			!
			!IMAG	
			write(filename, format) raw_dir//'AmnI.',qi
			open(unit=215, file=filename		, form='unformatted', access='stream', action='read',status='old') 
			read(215)	buffer
			close(215)
			Amn(:,:,qi)	= Amn(:,:,qi) + i_dp * dcmplx(buffer(:,:))
			!
		end do
		write(*,*)	"[read_Amn]: finished"
		!
		return
	end subroutine


	subroutine read_Mmn(Mmn)
		complex(dp),		intent(out)		::	Mmn(:,:,:,:)
		integer								::	qi
		character(len=20)					::	filename
		real(dp),		allocatable			::	buffer(:,:,:)
		!
		allocate(	buffer( size(Mmn,1), size(Mmn,2), size(Mmn,3) )		)
		!
		do qi = 1, size(Mmn,4)
			!REAL
			write(filename, format) raw_dir//'MmnR.',qi
			open(unit=210, file=filename		, form='unformatted', access='stream', action='read',status='old') 
			read(210)	buffer
			close(210)
			Mmn(:,:,:,qi)	= dcmplx(buffer(:,:,:))
			!
			!IMAG	
			write(filename, format) raw_dir//'MmnI.',qi
			open(unit=215, file=filename		, form='unformatted', access='stream', action='read',status='old') 
			read(215)	buffer
			close(215)
			Mmn(:,:,:,qi)	= Mmn(:,:,:,qi) + i_dp * dcmplx(buffer(:,:,:))
			!
		end do
		write(*,*)	"[read_Mmn]: finished"
		!
		return
	end subroutine


	subroutine read_velo(velo)
		complex(dp),	intent(out)			::	velo(:,:,:,:)
		integer								::	qi
		character(len=20)					::	filename
		real(dp),		allocatable			::	buffer(:,:,:)
		!
		allocate(	buffer( size(velo,1), size(velo,2), size(velo,3) )		)
		!
		velo = dcmplx(0.0_dp)
		do qi = 1, size(velo,4)
			!REAL
			write(filename, format) raw_dir//'velR.',qi
			open(unit=210, file=filename		, form='unformatted', access='stream', action='read',status='old') 
			read(210)	buffer
			close(210)
			velo(:,:,:,qi)	= dcmplx(buffer(:,:,:))
			!
			!IMAG
			write(filename, format) raw_dir//'velI.',qi
			open(unit=210, file=filename		, form='unformatted', access='stream', action='read',status='old') 
			read(210)	buffer
			close(210)
			velo(:,:,:,qi)	= velo(:,:,:,qi) + i_dp * dcmplx(buffer(:,:,:))
		end do
		write(*,*)	"[read_velo]: finished reading velR.qi & velI.qi files (abinit files)"
		!
		return
	end subroutine




end module basisIO












	!subroutine writeABiN_basCoeff(ck_loc)
	!	complex(dp),	intent(in)		::	ck_loc(:,:,:)
	!	real(dp),		allocatable		::	buffer(:,:)
	!	complex(dp),	allocatable		:: 	ck_glob(:,:,:)
	!	integer							:: 	qi, mesgSize
	!	!
	!	!ALLOCATE TARGET
	!	if( myID == root ) then
	!		allocate(	ck_glob( size(ck_loc,1), size(ck_loc,2), nQ		)		)
	!		allocate(	buffer(	size(ck_glob,1), size(ck_glob,2)		)		)
	!	else
	!		allocate(	ck_glob(0,0,0)	)
	!	end if
	!	!
	!	!GATHER
	!	mesgSize = size(ck_loc,1) * size(ck_loc,2) * size(ck_loc,3)
	!	call MPI_GATHER(	ck_loc, mesgSize, MPI_DOUBLE_COMPLEX, ck_glob, mesgSize, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr )
	!	!
	!	!WRITE TO FILE
	!	if(myID == root ) then
	!		!REAL PART
	!		open(unit=210, file=raw_dir//'ckR.dat'		, form='unformatted', access='stream', action='write',status='replace') 
	!		do qi = 1, nQ
	!			buffer	= dreal(ck_glob(:,:,qi)) 
	!			write(210)	buffer
	!		end do
	!		close(210)
	!		!IMAG PART
	!		open(unit=211, file=raw_dir//'ckI.dat'		, form='unformatted', access='stream', action='write', status='replace')
	!		do qi = 1, nQ
	!			buffer	= dimag(ck_glob(:,:,qi)) 
	!			write(211)	buffer
	!		end do
	!		close(211)
	!		write(*,*)	"[writeABiN_basCoeff]: wrote basis coefficients to ckR (real part) and ckI (imag part)"
	!		!
	!	end if
	!	!
	!	!DEBUG
	!	if( size(ck_loc,3) /= qChunk )	write(*,'(a,i3,a)')	"[#",myID,";writeABiN_basCoeff]: wrong qpt size of ck detected"
	!	!
	!	!
	!	return
	!end subroutine




