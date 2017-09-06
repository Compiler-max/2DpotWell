module output
	!module contains several routines for printing and writing data
	use mathematics,	only:	dp, PI_dp
	use sysPara 


	implicit none
	private

	public ::	writeMeshInfo, writeMeshBin, writeWaveFunc, writeWannFiles, writePolFile, &
				printMat, printTiming 


	interface printMat
		module procedure printCmplxMat
	end interface printMat

	contains















!public:
	subroutine writeMeshInfo()
		!writes the generated meshes to readable txt file for debuggin purpose
		!
		!
		integer		:: i, j
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
			write(100,'(a,i3,a,f7.5,a,f7.5,a,a,f7.5,a,f7.5,a)')	"atom(",i,") at (",atPos(1,i),", ",atPos(2,i),") ",&
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


	subroutine writeWaveFunc(unk, Aconn, Fcurv)
		complex(dp),	intent(in)		:: unk(:,:,:)
		real(dp),		intent(in)		:: Aconn(:,:,:,:), Fcurv(:,:,:)
		real(dp),		allocatable		:: unkR(:,:,:), unkI(:,:,:)
		!
		allocate(	unkR(	size(unk,1)		, size(unk,2)	, size(unk,3)		)			)
		allocate(	unkI(	size(unk,1)		, size(unk,2)	, size(unk,3)		)			)
		!
		unkR 	= dreal(unk)
		unkI 	= dimag(unk)
		!
		!LATTICE PERIODIC FUNCTIONS
		open(unit=400,file='rawData/unkR.dat',form='unformatted',access='stream',action='write')
		write(400)	unkR
		close(400)
		!
		open(unit=405,file='rawData/unkI.dat',form='unformatted',access='stream',action='write')
		write(405)	unkI
		close(405)
		!
		!
		!CONNECTION
		open(unit=410,file='rawData/AconnR.dat',form='unformatted',access='stream',action='write')
		write(410)	Aconn
		close(410)
		!
		!
		!CURVATURE	
		open(unit=420,file='rawData/FcurvR.dat',form='unformatted',access='stream',action='write')
		write(420) Fcurv		!420!!!!!!open(unit=410,file='rawData/AconnR.dat',form='unformatted',access='stream',action='write')
		close(420)
		!
		!
		return
	end subroutine



	subroutine writeWannFiles(wnF, wCent, wSprd)
		complex(dp),	intent(in)		:: wnF(:,:,:)
		real(dp),		intent(in)		:: wCent(:,:), wSprd(:,:)
		real(dp),		allocatable		:: wnfR(:,:,:), wnfI(:,:,:)
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
		!
		return
	end subroutine


	subroutine writeInterpolFiles(Aconn)
		complex(dp),	intent(in)		::	Aconn(:,:,:)
		real(dp),		allocatable		::	AconnR(:,:,:), AconnI(:,:,:)
		allocate(	AconnR(	size(Aconn,1)	, size(Aconn,2)	, size(Aconn,3)		)			)	
		allocate(	AconnI(	size(Aconn,1)	, size(Aconn,2)	, size(Aconn,3)		)			)
		!
		AconnR	= dreal(Aconn)
		AconnI	= dimag(Aconn)
		!
		open(unit=550,file='rawData/AintR.dat',form='unformatted',access='stream',action='write')
		write(550)	AconnR
		close(550)
		!
		open(unit=555,file='rawData/AintI.dat',form='unformatted',access='stream',action='write')
		write(555)	AconnI
		close(555)
		!
		return
	end subroutine


	subroutine writePolFile(pEl, pIon, pTot, pElA, pInt, pNiu, pPei )
		real(dp),		intent(in)		:: pEl(2), pIon(2), pTot(2), pElA(2), pInt(2), pNiu(3), pPei(3)
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
		write(600,'(a,f16.12,a,f16.12,a,f16.12,a)')	"pEl = ",norm2(pEl)	," * (", &	
																pEl(1)/norm2(pEl) 	,	", ",	pEl(2)/norm2(pEl),		")"
		write(600,'(a,f16.12,a,f16.12,a,f16.12,a)')	"pElA= ",norm2(pElA)," * (", &	
																pElA(1)/norm2(pElA)	,	", ",	pElA(2)/norm2(pelA),	")"
		write(600,'(a,f16.12,a,f16.12,a,f16.12,a)')	"pInt= ",norm2(pInt)," * (", &	
																pInt(1)/norm2(pInt),	", ",	pInt(2)/norm2(pInt),	")"
		write(600,'(a,f16.12,a,f16.12,a,f16.12,a)')	"pIon= ",norm2(pIon)," * (", &	
																pIon(1)/norm2(pIon)	,	", ",	pIon(2)/norm2(pIon),	")"
		write(600,'(a,f16.12,a,f16.12,a,f16.12,a)')	"pTot= ",norm2(pTot)," * (", &	
																pTot(1)/norm2(pTot),	", ",	pTot(2)/norm2(pTot),	")"
		!
		!
		write(600,*)"**************PERTURBATION:"
		write(600,'(a,f16.12,a,f16.12,a,f16.12,a,f16.12,a)')	"Bext= ",norm2(Bext) ," * (", &
											Bext(1)/norm2(Bext),	", ",	Bext(2)/norm2(Bext),", ", Bext(3)/norm2(Bext),	")"
		write(600,*)"*"
		write(600,'(a,f16.12,a,f16.12,a,f16.12,a,f16.12,a)')	"pNiu= ",norm2(pNiu)," * (", &	
											pNiu(1)/norm2(pNiu),	", ",	pNiu(2)/norm2(pNiu),", ", pNiu(3)/norm2(pNiu),	")"
		write(600,'(a,f16.12,a,f16.12,a,f16.12,a,f16.12,a)')	"pPei= ",norm2(pPei)," * (", &	
											pPei(1)/norm2(pPei),	", ",	pPei(2)/norm2(pPei),", ", pPei(3)/norm2(pPei),	")"
		close(600)
		!
		!
		return
	end subroutine




	subroutine printTiming(aT,kT,wI,scT,peiT,wT,oT,mastT)
		real,	intent(in)	:: aT, kT, wT, oT, wI, scT,peiT, mastT
		!
		print '    ("r&alloc  time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									aT 				, 100*aT		   	/mastT
		print '    ("k-solver time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									kT 				, 100*kT		   	/mastT
		print '    ("0 order pol. time spend = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									wT 				, 100*wT		   	/mastT
		print '    ("wannier interpolation   = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									wI 				, 100*wI		   	/mastT
		print '    ("semiclassics            = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									scT 			, 100*scT		   	/mastT							
		print '    ("peierls substitution    = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									peiT 			, 100*peiT		   	/mastT														
		print '    ("writing  time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									oT 				, 100*oT		   	/mastT
		print '    ("other    time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									(mastT-wT-kT-oT-wI-scT) 	, 100*(mastT-aT-wT-kT-oT-wI-scT)/mastT
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
	end












end module output
