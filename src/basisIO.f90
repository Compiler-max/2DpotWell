module basisIO
	use sysPara
	use mpi
	use omp_lib
	use mathematics,	only:	dp, i_dp, PI_dp, machineP, aUtoEv, aUtoAngstrm

	implicit none

	private
	public :: 		writeABiN_energy, writeABiN_basis, writeABiN_basCoeff, &
					readAbIn, readBasis





	contains





!WRITE
	subroutine writeABiN_energy(En_loc)
		real(dp),		intent(in)		:: 	En_loc(:,:)
		real(dp),		allocatable		:: 	En_glob(:,:)
		integer							::	mesgSize
		!
		!ALLOCATE TARGET
		if(myID == root )	allocate(	En_glob( size(En_loc,1), nQ	)		)
		if(myID /= root )	allocate(	En_glob(		0,		0	)		)
		!
		!GATHER
		mesgSize = size(En_loc,1)*size(En_loc,2)
		call MPI_GATHER( En_loc, mesgSize, MPI_DOUBLE_PRECISION, En_glob, mesgSize, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
		!
		!WRITE TO FILE
		if(myID == root ) then
			open(unit=200, file=raw_dir//'bandStruct.dat', form='unformatted', access='stream', action='write', status='unknown')
			write(200)	En_glob
			close(200)
			write(*,'(a,i3,a)')		"[#",myID,";writeABiN_energy]: wrote energy bands to file"
		end if
		!
		!DEBUG
		if( size(En_loc,2) /= qChunk )	write(*,'(a,i3,a)')		"[#",myID,";writeABiN_energy]: bands have wrong qpt size"
		!
		return
	end subroutine




	subroutine writeABiN_basCoeff(ck_loc)
		complex(dp),	intent(in)		::	ck_loc(:,:,:)
		real(dp),		allocatable		::	buffer(:,:)
		complex(dp),	allocatable		:: 	ck_glob(:,:,:)
		integer							:: 	qi, mesgSize
		!
		!ALLOCATE TARGET
		if( myID == root ) then
			allocate(	ck_glob( size(ck_loc,1), size(ck_loc,2), nQ		)		)
			allocate(	buffer(	size(ck_glob,1), size(ck_glob,2)		)		)
		else
			allocate(	ck_glob(0,0,0)	)
		end if
		!
		!GATHER
		mesgSize = size(ck_loc,1) * size(ck_loc,2) * size(ck_loc,3)
		call MPI_GATHER(	ck_loc, mesgSize, MPI_DOUBLE_COMPLEX, ck_glob, mesgSize, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr )
		!
		!WRITE TO FILE
		if(myID == root ) then
			!REAL PART
			open(unit=210, file=raw_dir//'ckR.dat'		, form='unformatted', access='stream', action='write',status='replace') 
			do qi = 1, size(ck_glob,3)
				buffer	= dreal(ck_glob(:,:,qi)) 
				write(210)	buffer
			end do
			close(210)
			!IMAG PART
			open(unit=211, file=raw_dir//'ckI.dat'		, form='unformatted', access='stream', action='write', status='replace')
			do qi = 1, size(ck_glob,3)
				buffer	= dimag(ck_glob(:,:,qi)) 
				write(211)	buffer
			end do
			close(211)
			write(*,*)	"[writeABiN_basCoeff]: wrote basis coefficients to ckR (real part) and ckI (imag part)"
			!
		end if
		!
		!DEBUG
		if( size(ck_loc,3) /= qChunk )	write(*,'(a,i3,a)')	"[#",myID,";writeABiN_basCoeff]: wrong qpt size of ck detected"
		!
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
		call MPI_GATHER( Gvec_loc	, qChunk, MPI_INTEGER, Gvec_glob	, qChunk, MPI_INTEGER, MPI_COMM_WORLD, ierr)	
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
		real(dp),		allocatable		:: buffer(:,:), eBuff(:)
		integer							:: qi
		!
		allocate(	buffer( size(ck,1), size(ck,2) 	)		)
		allocate(	eBuff(size(En,1)				)		)
		!
		!
		!call readGvec()
		!
		!UNK REAL PART
		open(unit=700, file=raw_dir//"ckR.dat",form='unformatted',access='stream',action='read')
		do qi = 1 , size(ck,3)
			read(700) buffer
			ck(:,:,qi)	= dcmplx(buffer)
		end do
		close(700)
		!
		!UNK IMAG PART
		open(unit=710, file=raw_dir//"ckI.dat",form='unformatted',access='stream',action='read')
		do qi = 1 , size(ck,3)
			read(710) buffer
			ck(:,:,qi)	= ck(:,:,qi) + i_dp * dcmplx(buffer)
		end do
		close(710)
		!
		!BAND ENERGIES
		open(unit=720, file=raw_dir//"bandStruct.dat",form='unformatted',access='stream',action='read')
		do qi = 1, size(En,2)	
				read(720) eBuff
				En(1:nSolve,qi)	= eBuff(1:nSolve)
		end do
		close(720)
		!
		!
	
		!
		return
	end subroutine


	subroutine readBasis()
		real(dp),		allocatable		:: buffer(:,:)
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