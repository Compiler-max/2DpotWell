module input
	use mathematics,	only:	dp, i_dp
	use sysPara
	use output,		 	only:		writeEnAndCK


	implicit none
	private

	public ::	filesExist, readHam

	contains





!public:
	logical function filesExist()
		logical				:: ckR, ckI, En
		!
		inquire( file="./rawData/ckR.dat"			, exist=ckR	)
		inquire( file="./rawData/ckI.dat"			, exist=ckI	)
		inquire( file="./rawData/bandStruct.dat"	, exist=En		)
		!
		if( .not. ckR ) write(*,*)"[filesExist]: real part of basCoeff not found"
		if( .not. ckI ) write(*,*)"[filesExist]: imag part of basCoeff not found"
		if( .not. En   ) write(*,*)"[filesExist]: energies not found"
		!
		filesExist = ( ckR .and. ckI )  .and. En
		!
		return
	end function


	subroutine	readHam(ck, En)
		complex(dp),	intent(out)		:: ck(:,:,:)
		real(dp),		intent(out)		:: En(:,:)
		real(dp),		allocatable		:: buffer(:,:,:), eBuff(:,:)
		!
		allocate(	buffer( size(ck,1), size(ck,2), size(ck,3) )		)
		allocate(	eBuff(nG,nQ)										)
		!
		!
		if( filesExist() ) then
			!
			!UNK REAL PART
			buffer	= 0.0_dp
			open(unit=700, file="./rawData/ckR.dat",form='unformatted',access='stream',action='read')
				read(700) buffer
			close(700)
			ck	= buffer
			!
			!UNK IMAG PART
			buffer	= 0.0_dp
			open(unit=710, file="./rawData/ckI.dat",form='unformatted',access='stream',action='read')
				read(710) buffer
			close(710)
			ck	= ck + i_dp * buffer
			!
			!BAND ENERGIES
			open(unit=720, file="./rawData/bandStruct.dat",form='unformatted',access='stream',action='read')
				read(720) eBuff
			close(720)
			En(1:nBands,:)	= eBuff(1:nBands,:)
			!
			!
			
			!call writeEnAndCK(eBuff, ck)
		else
			write(*,*)	"[readHam]: could not find all necessary files"
			ck	= dcmplx(0.0_dp)
			En	= 0.0_dp
		end if
		!
		!
		return
	end subroutine






end module input