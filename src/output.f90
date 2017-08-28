module output
	!module contains several routines for printing and writing data
	use mathematics,	only:	dp, PI_dp
	use sysPara    , 	only: 	aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot,&
								nG,nG0, Gcut, nKx, nKy,nK, nSC, nR, nRx, nRy, dx, dy, dkx, dky, nWfs, &
								Gvec, atPos, atR, kpts, rpts, Rcell, trialOrbVAL, Zion


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
	end


	subroutine writeMeshBin()
		!SYS PARA
		open(unit=300,file='rawData/sysPara.dat',form='unformatted',access='stream',action='write')
		write(300) nAt
		write(300) nG
		write(300) nK
		write(300) nKx
		write(300) nKy
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
		open(unit=320,file='rawData/kpts.dat',form='unformatted',access='stream',action='write')
		write(320) kpts
		close(320)
		!
		!
		return
	end


	subroutine writeWaveFunc(unk, Aconn)
		complex(dp),	intent(in)		:: unk(:,:,:), Aconn(:,:,:)
		!
		!
		!
		!LATTICE PERIODIC FUNCTIONS
		open(unit=400,file='rawData/unkR.dat',form='unformatted',access='stream',action='write')
		write(400)	dreal(unk)
		close(400)
		!
		open(unit=405,file='rawData/unkI.dat',form='unformatted',access='stream',action='write')
		write(405)	dimag(unk)
		close(405)
		!
		!
		!CONNECTION
		open(unit=410,file='rawData/AconnR.dat',form='unformatted',access='stream',action='write')
		write(410)	dreal(Aconn)
		close(410)
		!
		open(unit=415,file='rawData/AconnI.dat',form='unformatted',access='stream',action='write')
		write(415)	dimag(Aconn)
		close(415)
		!
		!
		return
	end



	subroutine writeWannFiles(wnF, wCent, wSprd)
		complex(dp),	intent(in)		:: wnF(:,:,:)
		real(dp),		intent(in)		:: wCent(:,:), wSprd(:,:)
		!
		!
		!
		!WANNIER FUNCTIONS:
		open(unit=500,file='rawData/wnfR.dat',form='unformatted',access='stream',action='write')
		write(500)	dreal(wnF)
		close(500)
		!
		open(unit=505,file='rawData/wnfI.dat',form='unformatted',access='stream',action='write')
		write(505)	dimag(wnF)
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
	end





	subroutine writePolFile(pEl, pIon, pTot, pElA )
		real(dp),		intent(in)		:: pEl(2), pIon(2), pTot(2), pElA(2)
		!	
		!	
		open(unit=600,file='polOutput.txt',action='write')
		write(600,*)"**************POLARIZATION OUTPUT FILE**********************"
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"ZION:"
		write(600,*) Zion
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"POL:"
		write(600,'(a,f16.12,a,f16.12,a)')	"pEl = (",	pEl(1) ,	", ",	pEl(2),		")"
		write(600,'(a,f16.12,a,f16.12,a)')	"pElA= (",	pElA(1),	", ",	pElA(2),	")"
		write(600,'(a,f16.12,a,f16.12,a)')	"pIon= (",	pIon(1),	", ",	pIon(2),	")"
		write(600,'(a,f16.12,a,f16.12,a)')	"pTot= (",	pTot(1),	", ",	pTot(2),	")"
		close(600)
		!
		return
	end




	subroutine printTiming(aT,kT,wT,oT,mastT)
		real,	intent(in)	:: aT, kT, wT, oT, mastT
		!
		print '    ("r&alloc  time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									aT 				, 100*aT		   	/mastT
		print '    ("k-solver time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									kT 				, 100*kT		   	/mastT
		print '    ("0 order pol. time spend = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									wT 				, 100*wT		   	/mastT
		print '    ("writing  time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									oT 				, 100*oT		   	/mastT
		print '    ("other    time spend     = ",f15.7," seconds = ",f15.7,"% of overall time")',& 
									(mastT-wT-kT) 	, 100*(mastT-aT-wT-kT-oT)/mastT
		print '    ("overall  time spend     = ",f15.7," seconds.")', mastT
		!
		return
	end














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
