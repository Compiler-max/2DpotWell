module output
	!module contains several routines for printing and writing data
	use mathematics,	only:	dp, PI_dp
	use blochWf,		only:	calcBasis
	use sysPara 


	implicit none
	private

	public ::	writeMeshInfo, writeMeshBin, writeEnAndCK, writeCkASunk, writeConnCurv, writeWannFiles, writePolFile, &
				printMat, printTiming , writePeierls,  writeInterpBands, writeEnH, printBasisInfo


	interface printMat
		module procedure printCmplxMat
	end interface printMat

	contains















!public:




	subroutine writeMeshInfo()
		!writes the generated meshes to readable txt file for debuggin purpose
		!
		!
		integer		:: i
		open(unit=100,file='meshInfo.txt',action='write')
		!
		!R MESH
		write(100,*)"*******************R POINT MESH******************************"
		do i = 1, nR
			write(100,'(a,i4,a,f15.10,a,f15.10,a)')				"r(",i,") = (", rpts(1,i) ,", ",rpts(2,i)," )"
		end do		
		!
		!R CELL	
		write(100,*)"*******************RCELL VECTOR******************************"
		do i = 1, nSC
			write(100,'(a,i4,a,f15.10,a,f15.10,a)')				"R(",i,") = (", Rcell(1,i) ,", ", Rcell(2,i), " )"
		end do

		!
		!G VECTOR
		write(100,*)"*******************G VECTOR MESH******************************"
		write(100,*)"nG0=",nG0
		do i = 1, nG
			write(100,'(a,i4,a,f15.12,a,f15.12,a)')				"G(",i,") = (", Gvec(1,i) , ", ", Gvec(2,i) , ")"
		end do
		!!
		!K MESH
		write(100,*)"*******************K POINT MESH******************************"
		do i = 1, nQ
			write(100,'(a,i4,a,f15.12,a,f15.12,a)')				"q(",i,") = (", qpts(1,i) , "," , qpts(2,i), " )"
		end do
		!!
		!K 	INTERPOLATION
		write(100,*)"*******************K INTERPOLATION MESH******************************"
		do i = 1, nK
			write(100,'(a,i4,a,f15.12,a,f15.12,a)')				"k(",i,") = (", kpts(1,i) , "," , kpts(2,i), " )"
		end do
		!!
		!ATOM POSITION AND RADII
		write(100,*)"*******************ATOMS******************************"
		do i = 1, nAt
			write(100,'(a,i3,a,f8.5,a,f8.5,a,a,f8.5,a,f8.5,a)')	"atom(",i,") at (",atPos(1,i),", ",atPos(2,i),") ",&
																"w. rad = (", atR(1,i), ", ", atR(2,i), " )"	
		end do
		!
		!TRIAL ORBITAL HEIGHT
		write(100,*)"*******************TIRAL ORBITALS******************************"
		do i = 1, nAt
			write(100,'(a,i2,a,f15.12)')				"g_trial(",i,") = ", trialOrbVAL(i) 
		end do
		!
		close(100)
		return
	end subroutine


	subroutine writeMeshBin()
		!SYS PARA
		open(unit=300,file='rawData/sysPara.dat',form='unformatted',access='stream',action='write')
		write(300) nAt
		write(300) nG
		write(300) nQ
		write(300) nQx
		write(300) nQy
		write(300) nR
		write(300) nRx
		write(300) nRy
		write(300) nWfs
		write(300) nSC
		write(300) nSCx
		write(300) nSCy
		write(300) nK
		write(300) nKx
		write(300) nKy
		close(300)
		!
		!CELL INFO
		open(unit=305,file='rawData/cellInfo.dat',form='unformatted',access='stream',action='write')
		write(305) aX
		write(305) aY
		close(305)
		!
		!
		!ATOM INFORMATION
		open(unit=306,file='rawData/atPos.dat',form='unformatted',access='stream',action='write')
		write(306) atPos
		close(306)
		!
		open(unit=307,file='rawData/atR.dat',form='unformatted',access='stream',action='write')
		write(307) atR
		close(307)
		!
		!
		!R MESH
		open(unit=310,file='rawData/rpts.dat',form='unformatted',access='stream',action='write')
		write(310) rpts
		close(310)
		!
		!K MESH
		open(unit=320,file='rawData/qpts.dat',form='unformatted',access='stream',action='write')
		write(320) qpts
		close(320)
		!
		!K INTERPOLATION
		open(unit=325,file='rawData/kpts.dat',form='unformatted',access='stream',action='write')
		write(325) kpts
		close(325)
		!
		!
		return
	end subroutine


	subroutine writeEnAndCK(EnT, ck)
		real(dp),		intent(in)		:: EnT(:,:)
		complex(dp),	intent(in)		:: ck(:,:,:)
		real(dp),		allocatable		:: buffer(:,:,:)
		!
		allocate(	buffer(	size(ck,1), size(ck,2), size(ck,3)	)		)
		!
		open(unit=200, file='rawData/bandStruct.dat', form='unformatted', access='stream', action='write')
		write(200)	EnT
		close(200)
		!
		buffer	= dreal(ck) 
		open(unit=210, file='rawData/ckR.dat'		, form='unformatted', access='stream', action='write',status='replace') 
		write(210)	buffer
		close(210)
		!
		buffer	= dimag(ck)
		open(unit=211, file='rawData/ckI.dat'		, form='unformatted', access='stream', action='write')
		write(211)	buffer
		write(*,*)	"[writeEnAndCK]: wrote eigenvalues and eigencoefficients"

		return
	end subroutine


	subroutine writeCkASunk(ck, ckW)
		complex(dp),	intent(in)		:: ck(:,:,:), ckW(:,:,:)
		complex(dp),	allocatable		:: basVec(:), unk(:,:,:), unkW(:,:,:)
		real(dp),		allocatable		:: buffer(:,:,:), bufferW(:,:,:)
		integer							:: qi, ri
		!
		allocate( basVec(nG)		)
		allocate( unk(nR,nBands,nQ)	)
		allocate( unkW(nR,nWfs,nQ)	)
		allocate( buffer(nR,nBands,nQ))
		allocate( bufferW(nR,nBands,nQ))

		!
		!CALCULATE UNKs on real space grid
		do qi = 1, nQ
			do ri = 1, nR
				call calcBasis(qi, ri, basVec)
				unk(ri,:,qi) = matmul(basVec, ck(:,:,qi)	)
				unkW(ri,:,qi)= matmul(basVec, ckW(:,:,qi)	)
			end do
		end do
		!
		!WRITE RESULTS
		!UNK
		buffer	= dreal(unk)
		open(unit=400,file='rawData/unkR.dat',form='unformatted',access='stream',action='write')
		write(400)	buffer
		close(400)
		!
		buffer	= dimag(unk)
		open(unit=405,file='rawData/unkI.dat',form='unformatted',access='stream',action='write')
		write(405)	buffer
		close(405)
		!
		!UNK (W)
		bufferW	= dreal(unkW)
		open(unit=410,file='rawData/ROTunkR.dat',form='unformatted',access='stream',action='write')
		write(410)	bufferW
		close(410)
		!
		bufferW	= dimag(unkW)
		open(unit=415,file='rawData/ROTunkI.dat',form='unformatted',access='stream',action='write')
		write(415)	bufferW
		close(415)
		!
		!
		return
	end subroutine


	subroutine writeUNKs(unk)
		complex(dp),	intent(in)		:: unk(:,:,:)
		real(dp),		allocatable		:: buffer(:,:,:)
		!
		allocate(	buffer(size(unk,1), size(unk,2), size(unk,3))	)
		!REAL PART
		buffer	= dreal(unk)
		open(unit=400,file='rawData/ROTunkR.dat',form='unformatted',access='stream',action='write')
		write(400)	buffer
		close(400)
		!IMAG PART
		buffer	= dimag(unk)
		open(unit=405,file='rawData/ROTunkI.dat',form='unformatted',access='stream',action='write')
		write(405)	buffer
		close(405)
		!!
		return
	end subroutine


	subroutine writeConnCurv(Aconn, Fcurv)
		complex(dp),		intent(in)		:: Aconn(:,:,:,:), Fcurv(:,:,:,:)
		real(dp),			allocatable		:: buffer(:,:,:)
		integer								:: ki
		!
		allocate(	buffer( size(Aconn,1),size(Aconn,2),size(Aconn,3) )		)
		!
		!CONNECTION
		open(unit=410,file='rawData/AconnR.dat',form='unformatted',access='stream',action='write')
		do ki = 1, size(Aconn,4)
			buffer	= dreal(Aconn(:,:,:,ki))
			write(410)	buffer
		end do
		close(410)
		!
		!
		!CURVATURE
		open(unit=420,file='rawData/FcurvR.dat',form='unformatted',access='stream',action='write')
		do ki = 1, size(Fcurv,4)
			buffer	= dreal(Fcurv(:,:,:,ki))
			write(420) buffer
		end do	
		close(420)
		!
		!
		!
		return
	end subroutine



	subroutine writeWannFiles(wnF, wCent, wSprd)
		complex(dp),	intent(in)		:: wnF(:,:,:)
		real(dp),		intent(in)		:: wCent(:,:), wSprd(:,:)
		real(dp),		allocatable		:: wnfR(:,:,:), wnfI(:,:,:)
		integer							:: n
		!
		allocate(	wnfR(	size(wnF,1), size(wnF,2), size(wnF,3)		)				)
		allocate(	wnfI(	size(wnF,1), size(wnF,2), size(wnF,3)		)				)
		!
		wnfR	= dreal(wnF)
		wnfI	= dimag(wnF)
		!
		!WANNIER FUNCTIONS:
		open(unit=500,file='rawData/wnfR.dat',form='unformatted',access='stream',action='write')
		write(500)	wnfR
		close(500)
		!
		open(unit=505,file='rawData/wnfI.dat',form='unformatted',access='stream',action='write')
		write(505)	wnfI
		close(505)
		!
		!
		!
		!CENTERS AND SPREADS:
		open(unit=510,file='rawData/wCent.dat',form='unformatted',access='stream',action='write')
		write(510)	wCent
		close(510)
		!
		open(unit=515,file='rawData/wSprd.dat',form='unformatted',access='stream',action='write')
		write(515)	wCent
		close(515)
		!
		!TEXT FILE
		open(unit=516,file='wannier.txt',action='write')
		write(516,'(a)')	"****************atom positions****************************"
		do n = 1, nAt
			write(516,'(a,i3,a,f6.3,a,f6.3,a)')	"atom=,",n,	"centered at (",atPos(1,n),", ",atPos(2,n),")."
		end do
		!
		write(516,'(a)')	"****************Wannier functions****************************"
		do n = 1, nWfs
			write(516,'(a,i3,a,f10.5,a,f10.5,a,f10.5,a,f10.5)') &	
					"n=",n	,"wCent= (",wCent(1,n),", ",wCent(2,n), &
					").wSprd=(",wSprd(1,n),", ",wSprd(2,n)
			write(516,'(a,i3,a,f10.5,a,f10.5,a,f10.5,a,f10.5)')	&
					"n=",n	,"wCent= (",dmod(wCent(1,n),aX),", ",dmod(wCent(2,n),aY), & 
					").wSprd=(",wSprd(1,n),", ",wSprd(2,n)
		!															
		end do
		close(516)
		!
		!
		return
	end subroutine


	subroutine writeInterpBands(Ew)
		real(dp),		intent(in)		:: Ew(:,:)
		open(unit=520,file='rawData/Ewann.dat',form='unformatted',access='stream',action='write')
		write(520)	Ew
		close(520)
		!
		!
		return
	end subroutine


	subroutine writePeierls(ckP, EnP)
		complex(dp),	intent(in)		:: ckP(:,:,:)
		real(dp),		intent(in)		:: EnP(:,:)
		real(dp),		allocatable		:: buffer(:,:,:)
		!
		allocate(	buffer( size(ckP,1)	,	size(ckP,2)	,	size(ckP,3) 					)			)
		
		!real(UNK)
		buffer	= dreal(ckP)
		open(unit=700,file='rawData/ckPeiR.dat',form='unformatted',access='stream',action='write')
		write(700)	buffer
		close(700)
		!imag(UNK)
		buffer	= dimag(ckP)
		open(unit=705,file='rawData/ckPeiI.dat',form='unformatted',access='stream',action='write')
		write(705)	buffer
		close(705)
		!
		!Conn
		open(unit=710,file='rawData/EnPei.dat',form='unformatted',access='stream',action='write')
		write(710)	EnP
		close(710)
		!
		
		!
		return
	end subroutine


	subroutine writeEnH(EnH)
		real(dp),		intent(in)		:: EnH(:,:)
		!
		open(unit=800,file='rawData/EnInterP.dat',form='unformatted',access='stream',action='write')
		write(800)	EnH
		close(800)
		!
		return
	end subroutine


	subroutine writePolFile(pWann, pBerry, pNiu, pPei )
		real(dp),		intent(in)		:: pWann(2), pBerry(2), pNiu(3), pPei(3)
		!	
		!	
		open(unit=600,file='polOutput.txt',action='write')
		write(600,*)"**************POLARIZATION OUTPUT FILE**********************"
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"**************ZION:"
		write(600,*) Zion
		write(600,*)"*"
		write(600,*)"*"
		!
		write(600,*)"**************POL:"
	

	
		write(600,'(a,f8.4,a,f8.4,a,a,f8.4,a,f8.4,a)')	"pWann =  (",  pWann(1)	,	", ",	pWann(2),		"),",& 
												" moded=(",dmod(pWann(1),aX/vol),", ",dmod(pWann(2),aY/vol),")."
		!
		write(600,'(a,f8.4,a,f8.4,a,a,f8.4,a,f8.4,a)')	"pBerry=  (",		pBerry(1)	,	", ",	pBerry(2),		"),",& 
												" moded=(",dmod(pBerry(1),aX/vol),", ",dmod(pBerry(2),aY/vol),")."
		!write(600,'(a,e16.9,a,f16.12,a,f16.12,a)')	"pInt= ",norm2(pInt)," * (", &	
		!														pInt(1)/norm2(pInt),	", ",	pInt(2)/norm2(pInt),	")"
		!write(600,'(a,e16.9,a,f16.12,a,f16.12,a)')	"pIon= ",norm2(pIon)," * (", &	
		!														pIon(1)/norm2(pIon)	,	", ",	pIon(2)/norm2(pIon),	")"
		!write(600,'(a,e16.9,a,f16.12,a,f16.12,a)')	"pTot= ",norm2(pTot)," * (", &	
		!														pTot(1)/norm2(pTot),	", ",	pTot(2)/norm2(pTot),	")"
		!
		!
		write(600,*)"**************PERTURBATION:"
		write(600,'(a,f16.12,a,f16.12,a,f16.12,a,f16.12,a)')	"Bext= ",norm2(Bext) ," * (", &
											Bext(1)/norm2(Bext),	", ",	Bext(2)/norm2(Bext),", ", Bext(3)/norm2(Bext),	")"
		write(600,*)"*"
		!write(600,'(a,e16.9,a,f16.12,a,f16.12,a,f16.12,a)')	"pNiu= ",norm2(pNiu)," * (", &	
		!									pNiu(1)/norm2(pNiu),	", ",	pNiu(2)/norm2(pNiu),", ", pNiu(3)/norm2(pNiu),	")"

		write(600,'(a,f16.8,a,f16.8,a,f16.8,a)')	"pNiu= (", 	pNiu(1),	", ",	pNiu(2),", ", pNiu(3),	")"

		write(600,'(a,f16.8,a,f16.8,a,f16.8,a)')	"pPei= (", 	pPei(1),	", ",	pPei(2),", ", pPei(3),	")"
	
		close(600)
		!
		!
		return
	end subroutine



	subroutine printBasisInfo()
		!nG to Gcut RATIO
		if(			nG		< 		vol * Gcut * dsqrt(Gcut) /	(2.0_dp*PI_dp**2)		) then
			write(*,*)	"[main]: increase nG or decrease Gcut"
		else
			write(*,*)	"[main]: it seem to be enough basis functions for the choosen Gcut parameter."
		end if
		!nR to Gcut RATIO
		if(		nR	<		vol * dsqrt(Gcut) 	/ PI_dp) then
			write(*,*)	"[main]: need more real space points or smaller Gcut"
		else
			write(*,*)	"[main]: real space mesh seems fine enough"
		end if
		!real space grid TO supercells
		if(		nRx / nSCx <  nQx) write(*,*)	"[main]: need more reals space points or less k points"

		return 
	end subroutine



	subroutine printTiming(aT,kT,pT,wT,bT,peiT,oT,mastT)
		real,	intent(in)	:: aT,kT,pT,wT,bT,peiT,oT,mastT
		!
		print '    ("r&alloc  time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									aT 				, 100*aT		   	/mastT
		print '    ("k-solver time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									kT 				, 100*kT		   	/mastT
		print '    ("project. time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									pT 				, 100*pT		   	/mastT
		print '    ("wannier method          = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									wT				, 100*wT		   	/mastT
		print '    ("berry method            = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									bT				, 100*bT		   	/mastT							
		print '    ("peierls method          = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									peiT			, 100*peiT		   	/mastT													
		print '    ("writing  time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									oT 				, 100*oT		   	/mastT
		print '    ("other    time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									(mastT-aT-kT-pT-wT-bT-oT) 	, 100*(mastT-aT-kT-pT-wT-bT-oT)/mastT
		print '    ("overall  time spend     = ",f15.7," seconds.")', mastT
		!
		return
	end subroutine














!privat:	
	subroutine printCmplxMat(n, M)
		integer    		, intent(in)    :: n
		complex(dp)		, intent(in)    :: M(:,:)
		character(len=3)				:: imag_unit(n,n) 
		integer 						:: i1,i2
		!
		imag_unit = '+i*'
		WHERE(AIMAG(M)<0.)imag_unit = '-i*'
		!
		do i2 = 1,n
			write(*,'(100(a,f7.3,a,f7.3,a,Xxxx))') 		 ('(',real(M(i2,i1)),imag_unit(i2,i1),abs(aimag(M(i2,i1))),')',i1=1,n )
		end do
		!
		return
	end subroutine












end module output
