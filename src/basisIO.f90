module basisIO
	use sysPara
	use mpi
	use omp_lib
	use mathematics,	only:	dp, i_dp, PI_dp, machineP, aUtoEv, aUtoAngstrm

	implicit none

	private
	public :: 		writeABiN_energy, writeABiN_basis, writeABiN_basCoeff, &
					readAbIn, readBasis


	character(len=1024)				::	format='(a,i7.7)'


	contains





!WRITE
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
	subroutine	readAbIn(ck, En)
		complex(dp),	intent(out)		:: ck(:,:,:)
		real(dp),		intent(out)		:: En(:,:)
		real(dp),		allocatable		:: buffer(:,:)
		integer							:: qi
		character(len=20)				:: filename
		!
		allocate(	buffer( size(ck,1), size(ck,2) 	)		)
		!
		!UNK REAL PART
		do qi = 1 , size(ck,3)
			write(filename, format) raw_dir//'ckR.',qi
			open(unit=700, file=filename ,form='unformatted',access='stream',action='read')
			read(700) buffer
			ck(:,:,qi)	= dcmplx(buffer)
			close(700)
		end do
		!
		!UNK IMAG PART
		do qi = 1 , size(ck,3)
			write(filename, format) raw_dir//'ckI.',qi
			open(unit=710, file=filename,form='unformatted',access='stream',action='read')
			read(710) buffer
			ck(:,:,qi)	= ck(:,:,qi) + i_dp * dcmplx(buffer)
			close(710)
		end do
		write(*,*)	"[readAbIn]: read the basis coefficients"
		!
		!BAND ENERGIES
		do qi = 1, size(En,2)	
				write(filename, format) raw_dir//'enK.',qi
				open(unit=720, file=filename, form='unformatted', access='stream', action='read')
				read(720)	En(:,qi)
				close(720)
		end do
		write(*,*)	"[readAbIn]: read eigenvalues"
		!
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





end module basisIO