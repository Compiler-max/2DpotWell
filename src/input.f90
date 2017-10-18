module input
	use mathematics,	only:	dp, i_dp
	use sysPara
	use output,		 	only:		writeEnAndUNK 


	implicit none
	private

	public ::	filesExist, readHam

	contains





!public:
	logical function filesExist()
		logical				:: unkR, unkI, En
		!
		inquire( file="./rawData/unkR.dat"			, exist=unkR	)
		inquire( file="./rawData/unkI.dat"			, exist=unkI	)
		inquire( file="./rawData/bandStruct.dat"	, exist=En		)
		!
		if( .not. unkR ) write(*,*)"[filesExist]: real part of unk not found"
		if( .not. unkI ) write(*,*)"[filesExist]: imag part of unk not found"
		if( .not. En   ) write(*,*)"[filesExist]: energies not found"
		!
		filesExist = ( unkR .and. unkI )  .and. En
		!
		return
	end function


	subroutine	readHam(unk, En)
		complex(dp),	intent(out)		:: unk(:,:,:)
		real(dp),		intent(out)		:: En(:,:)
		real(dp),		allocatable		:: buffer(:,:,:), eBuff(:,:)
		!
		allocate(	buffer( size(unk,1), size(unk,2), size(unk,3) )		)
		allocate(	eBuff(nG,nQ)										)
		!
		!
		if( filesExist() ) then
			!
			!UNK REAL PART
			buffer	= 0.0_dp
			open(unit=700, file="./rawData/unkR.dat",form='unformatted',access='stream',action='read')
				read(700) buffer
			close(700,status='delete')
			unk	= buffer
			!
			!UNK IMAG PART
			buffer	= 0.0_dp
			open(unit=710, file="./rawData/unkI.dat",form='unformatted',access='stream',action='read')
				read(710) buffer
			close(710,status='delete')
			unk	= unk + i_dp * buffer
			!
			!BAND ENERGIES
			open(unit=720, file="./rawData/bandStruct.dat",form='unformatted',access='stream',action='read')
				read(720) eBuff
			close(720,status='delete')
			En(1:nBands,:)	= eBuff(1:nBands,:)
			!
			!
			call writeEnAndUNK(eBuff, unk)
		else
			write(*,*)	"[readHam]: could not find all necessary files"
			unk	= dcmplx(0.0_dp)
			En	= 0.0_dp
		end if
		!
		!
		return
	end subroutine






end module input