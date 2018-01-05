module output
	!module contains several routines for printing and writing data
	use omp_lib
	use mathematics,	only:	dp, PI_dp, machineP, aUtoEv, aUtoAngstrm
	use blochWf,		only:	calcBasis
	use sysPara 


	implicit none
	private

	public ::	writeMeshInfo, writeMeshBin, writeEnAndCK, writeCkASunk, writeConnCurv, writeWannFiles, writePolFile, &
				printMat, printTiming , writePeierls,  writeInterpBands, writeEnH, printBasisInfo, & 
				writeVeloHtxt, writeVeloEffTB, writeUmat, writeInterpU, writeBerryInterpU,& 
				writeHtb, writeHtbBerry, writeRtbBerry, writeEnEffTB, writeEnBerry, writeEnAbInitio


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
			write(100,'(a,i6,a,f15.7,a,f15.7,a)')				"r(",i,") = (", rpts(1,i) ,", ",rpts(2,i)," )"
		end do		
		!
		!R CELL	
		write(100,*)"*******************RCELL VECTOR******************************"
		do i = 1, nSC
			write(100,'(a,i6,a,f15.7,a,f15.7,a)')				"R(",i,") = (", Rcell(1,i) ,", ", Rcell(2,i), " )"
		end do

		!
		!G VECTOR
		write(100,*)"*******************G VECTOR  test MESH******************************"
		write(100,*)"nG0=",nG0
		do i = 1, nG
			write(100,'(a,i6,a,f12.4,a,f12.4,a)')				"G(",i,") = (", Gtest(1,i) , ", ", Gtest(2,i) , ")"
		end do
		!!
		!K MESH
		write(100,*)"*******************K POINT MESH******************************"
		do i = 1, nQ
			write(100,'(a,i6,a,f15.7,a,f15.7,a)')				"q(",i,") = (", qpts(1,i) , "," , qpts(2,i), " )"
		end do
		!!
		!K 	INTERPOLATION
		write(100,*)"*******************K INTERPOLATION MESH******************************"
		do i = 1, nK
			write(100,'(a,i6,a,f15.7,a,f15.7,a)')				"k(",i,") = (", kpts(1,i) , "," , kpts(2,i), " )"
		end do
		!!
		!ATOM POSITION AND RADII
		write(100,*)"*******************ATOMS******************************"
		do i = 1, nAt
			write(100,'(a,i3,a,f8.3,a,f8.3,a,a,f8.3,a,f8.3,a)')	"atom(",i,") at (",atPos(1,i),", ",atPos(2,i),") ",&
																"w. rad = (", atR(1,i), ", ", atR(2,i), " )"	
		end do
		!
		!
		close(100)
		return
	end subroutine


	subroutine writeMeshBin()
		!SYS PARA
		integer				:: stat
		open(unit=300, iostat=stat, file='rawData/sysPara.dat', status='old')
		if (stat == 0) close(300, status='delete')
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
		write(300) nSolve
		close(300)
		!
		!CELL INFO
		open(unit=305, iostat=stat, file='rawData/cellInfo.dat', status='old')
		if (stat == 0) close(305, status='delete')
		open(unit=305,file='rawData/cellInfo.dat',form='unformatted',access='stream',action='write')
		write(305) aX
		write(305) aY
		close(305)
		!
		!
		!ATOM INFORMATION
		open(unit=306, iostat=stat, file='rawData/atPos.dat', status='old')
		if (stat == 0) close(306, status='delete')
		open(unit=306,file='rawData/atPos.dat',form='unformatted',access='stream',action='write')
		write(306) atPos
		close(306)
		!
		open(unit=307, iostat=stat, file='rawData/atR.dat', status='old')
		if (stat == 0) close(307, status='delete')
		open(unit=307,file='rawData/atR.dat',form='unformatted',access='stream',action='write')
		write(307) atR
		close(307)
		!
		!
		!R MESH
		open(unit=310, iostat=stat, file='rawData/rpts.dat', status='old')
		if (stat == 0) close(310, status='delete')
		open(unit=310,file='rawData/rpts.dat',form='unformatted',access='stream',action='write')
		write(310) rpts
		close(310)
		!
		!K MESH
		open(unit=320, iostat=stat, file='rawData/qpts.dat', status='old')
		if (stat == 0) close(320, status='delete')
		open(unit=320,file='rawData/qpts.dat',form='unformatted',access='stream',action='write')
		write(320) qpts
		close(320)
		!
		!K INTERPOLATION
		open(unit=325, iostat=stat, file='rawData/kpts.dat', status='old')
		if (stat == 0) close(325, status='delete')
		open(unit=325,file='rawData/kpts.dat',form='unformatted',access='stream',action='write')
		write(325) kpts
		close(325)
		!
		!
		return
	end subroutine


	subroutine writeEnAndCK(EnT, ck, nGq)
		real(dp),		intent(in)		:: EnT(:,:)
		complex(dp),	intent(in)		:: ck(:,:,:)
		integer,		intent(in)		:: nGq(:)
		real(dp),		allocatable		:: buffer(:,:)
		integer							:: qi
		!
		!
		allocate(	buffer(	size(ck,1), size(ck,2)	)		)
		!
		open(unit=200, file='rawData/bandStruct.dat', form='unformatted', access='stream', action='write', status='replace')
		write(200)	EnT
		close(200)
		!
	
		open(unit=210, file='rawData/ckR.dat'		, form='unformatted', access='stream', action='write',status='replace') 
		do qi = 1, size(ck,3)
			buffer	= dreal(ck(:,:,qi)) 
			write(210)	buffer
		end do
		close(210)
		!
		open(unit=211, file='rawData/ckI.dat'		, form='unformatted', access='stream', action='write', status='replace')
		do qi = 1, size(ck,3)
			buffer	= dimag(ck(:,:,qi)) 
			write(211)	buffer
		end do
		close(211)
		!
		open(unit=215, file='rawData/nGq.dat'		, form='unformatted', access='stream', action='write', status='replace')
		write(215)	nGq
		close(215)

		write(*,*)	"[writeEnAndCK]: wrote eigenvalues and eigencoefficients and nGq info"
		!
		return
	end subroutine


	subroutine writeCkASunk(ck, ckW)
		complex(dp),	intent(in)		:: ck(:,:,:), ckW(:,:,:)
		complex(dp),	allocatable		:: basVec(:), unk(:,:,:), unkW(:,:,:)
		integer							:: qi, ri, stat, Gmax
		!
		allocate( unk(nR,nBands,nQ)	)
		allocate( unkW(nR,nWfs,nQ)	)
		write(*,*)"[writeCkASunk]: allocated arrays "
		!
		!CALCULATE UNKs on real space grid
		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(qi, ri, basVec, Gmax ) 
		allocate(	basVec(nG)		)
		!$OMP DO SCHEDULE(STATIC) COLLAPSE(1) 
		do qi = 1, nQ
			Gmax	= nGq(qi)
			do ri = 1, nR
				call calcBasis(qi, ri, basVec)
				unk(ri,:,qi) = matmul(basVec(1:Gmax),  ck(1:Gmax,:,qi)	)
				unkW(ri,:,qi)= matmul(basVec(1:Gmax), ckW(1:Gmax,:,qi)	)
			end do
		end do
		!$OMP END DO
		!$OMP END PARALLEL
		write(*,*)"[writeCkASunk]: calculated unks on real space grid"
		!
		!WRITE RESULTS
		!UNK
		open(unit=400, iostat=stat, file='rawData/unkR.dat', status='old')
		if (stat == 0) close(400, status='delete')
		open(unit=405, iostat=stat, file='rawData/unkI.dat', status='old')
		if (stat == 0) close(405, status='delete')
		open(unit=410, iostat=stat, file='rawData/ROTunkR.dat', status='old')
		if (stat == 0) close(410, status='delete')
		open(unit=415, iostat=stat, file='rawData/ROTunkI.dat', status='old')
		if (stat == 0) close(415, status='delete')


		open(unit=400,file='rawData/unkR.dat',form='unformatted',access='stream',action='write')
		open(unit=405,file='rawData/unkI.dat',form='unformatted',access='stream',action='write')
		open(unit=410,file='rawData/ROTunkR.dat',form='unformatted',access='stream',action='write')
		open(unit=415,file='rawData/ROTunkI.dat',form='unformatted',access='stream',action='write')
		do qi = 1, nQ
			write(400) dreal(unk(:,:,qi))
			write(405) dimag(unk(:,:,qi))
			!
			write(410) dreal(unkW(:,:,qi))
			write(415) dimag(unkW(:,:,qi))
			!buffer(:,:,qi)	= dreal(unk(:,:,qi))
		end do
		close(400)
		close(405)
		close(410)
		close(415)
		!
		write(*,*)"[writeCkASunk]: wrote all files"
		!
		!
		return
	end subroutine



	subroutine writeUmat(U_mat)	
		complex(dp),		intent(in)		:: U_mat(:,:,:)
		integer								:: qi, n, m
		!
		open(unit=800,file='U_mat.txt',action='write')
		write(800,*)	"U_mat after read in from berry"
		do qi = 1, size(U_mat,3)
			write(800,*)	qpts(1,qi)/recpLatt(1,1)," ",qpts(2,qi)/recpLatt(2,2)
			do n = 1, size(U_mat,2)
				do m = 1, size(U_mat,1)
					write(800,'(a,i3,a,i3,a,f14.10,a,f14.10)')	" ",m," ",n," ",dreal(U_mat(m,n,qi))," ",dimag(U_mat(m,n,qi))
				end do
			end do
			write(800,*)
		end do 
		close(800)
		!
		return
	end subroutine


	subroutine writeInterpU(U_mat)
		complex(dp),	intent(in)		:: U_mat(:,:,:)
		integer								:: ki, n, m
		!
		open(unit=805,file='U_matTB.txt',action='write')
		write(805,*)	"U_mat interpolated by effTB"
		do ki = 1, size(U_mat,3)
			write(805,*)	kpts(1,ki)/recpLatt(1,1)," ",kpts(2,ki)/recpLatt(2,2)
			do n = 1, size(U_mat,2)
				do m = 1, size(U_mat,1)
					write(805,'(a,i3,a,i3,a,f14.10,a,f14.10)')	" ",m," ",n," ",dreal(U_mat(m,n,ki))," ",dimag(U_mat(m,n,ki))
				end do
			end do
			write(805,*)
		end do 
		close(805)
		!
		return
	end subroutine


	subroutine writeBerryInterpU( U_mat)
		complex(dp),	intent(in)		:: U_mat(:,:,:)
		integer								:: ki, n, m
		!
		open(unit=805,file='U_matBerryInterp.txt',action='write')
		write(805,*)	"U_mat interpolated by berry"
		do ki = 1, size(U_mat,3)
			write(805,*)	kpts(1,ki)/recpLatt(1,1)," ",kpts(2,ki)/recpLatt(2,2)
			do n = 1, size(U_mat,2)
				do m = 1, size(U_mat,1)
					write(805,'(a,i3,a,i3,a,f14.10,a,f14.10)')	" ",m," ",n," ",dreal(U_mat(m,n,ki))," ",dimag(U_mat(m,n,ki))
				end do
			end do
			write(805,*)
		end do 
		close(805)
		!
		return
	end subroutine

	subroutine writeEnAbInitio(En)
		real(dp),		intent(in)			:: En(:,:)
		integer								:: ki, n
		!
		open(unit=805,file='En_AbInitio.txt',action='write')
		write(805,*)	"energies from solver"
		write(805,*)	size(En,2)
		do ki = 1, size(En,2)
			write(805,*)	"#ki=",ki
			do n = 1, size(En,1)
				write(805,*)	En(n,ki)
			end do
		end do
		close(805)
		!
		return
	end subroutine

	subroutine writeEnEffTB( En)
		real(dp),		intent(in)			:: En(:,:)
		integer								:: ki, n
		!
		open(unit=805,file='En_InterpTB.txt',action='write')
		write(805,*)	"interpolated energies by pW90 in a.u."
		write(805,*)	size(En,2)
		do ki = 1, size(En,2)
			write(805,*)	"#ki=",ki
			do n = 1, size(En,1)
				write(805,*)	En(n,ki)
			end do
		end do
		close(805)
		!
		return
	end subroutine

	subroutine writeEnBerry(En)
		real(dp),		intent(in)			:: En(:,:)
		integer								:: ki, n
		!
		open(unit=805,file='En_InterpBerry.txt',action='write')
		write(805,*)	"interpolated energies by berry in a.u."
		write(805,*)	size(En,2)
		do ki = 1, size(En,2)
			write(805,*)	"#ki=",ki
			do n = 1, size(En,1)
				write(805,*)	En(n,ki)
			end do
		end do
		close(805)
		!	
		return
	end subroutine


	subroutine writeHtb( H_tb )
		complex(dp),		intent(in)		:: 	H_tb(:,:,:)
		integer								::	R, n, m

		open(unit=805,file='H_tb.txt',action='write')
		write(805,*)	"H_tb (eV) read in by effTB"
		do R = 1, size(H_tb,3)
			write(805,*)	Rcell(1,R), " ", Rcell(2,R)
			do n = 1, size(H_tb,2)
				do m = 1, size(H_tb,1)
					write(805,'(a,i3,a,i3,a,e16.8,a,e16.8)')	" ",m," ",n," ",dreal(H_tb(m,n,R)*aUtoEv)," ",dimag(H_tb(m,n,R)*aUtoEv)
				end do
			end do
			write(805,*)
		end do 
		close(805)

		return
	end subroutine

	subroutine writeHtbBerry( H_tb )
		complex(dp),		intent(in)		:: 	H_tb(:,:,:)
		integer								::	R, n, m

		open(unit=807,file='H_tbBerry.txt',action='write')
		write(807,*)	"H_tb (eV) calculated by Berry"
		do R = 1, size(H_tb,3)
			write(807,*)	Rcell(1,R), " ", Rcell(2,R)
			do n = 1, size(H_tb,2)
				do m = 1, size(H_tb,1)
					write(807,'(a,i3,a,i3,a,e16.8,a,e16.8)')	" ",m," ",n," ",dreal(H_tb(m,n,R)*aUtoEv)," ",dimag(H_tb(m,n,R)*aUtoEv)
				end do
			end do
			write(807,*)
		end do 
		close(805)

		return
	end subroutine


	subroutine writeRtbBerry( r_tb )
		complex(dp),		intent(in)		:: 	r_tb(:,:,:,:)
		integer								::	R, n, m

		open(unit=807,file='r_tbBerry.txt',action='write')
		write(807,*)	"r_tb (angstroem) calculated by Berry"
		do R = 1, size(r_tb,4)
			write(807,*)	Rcell(1,R), " ", Rcell(2,R)
			do n = 1, size(r_tb,3)
				do m = 1, size(r_tb,2)
					write(807,'(a,i3,a,i3,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8,a,e16.8)')	&
							" ",m," ",n," ",dreal(r_tb(1,m,n,R)*aUtoAngstrm)," ",dimag(r_tb(1,m,n,R)*aUtoAngstrm),&
										" ",dreal(r_tb(2,m,n,R)*aUtoAngstrm)," ",dimag(r_tb(2,m,n,R)*aUtoAngstrm),&
										" ",dreal(r_tb(3,m,n,R)*aUtoAngstrm)," ",dimag(r_tb(3,m,n,R)*aUtoAngstrm)
				end do
			end do
			write(807,*)
		end do 
		close(805)

		return
	end subroutine



	subroutine writeVeloHtxt(velo)
		complex(dp),	intent(in)		:: velo(:,:,:,:) !veloK(3		, 	nWfs	,	nWfs	,	nK		)
		integer		qi, n, m
		complex(dp)						:: vSum
		real(dp)						:: Ederiv

		Ederiv	= aUtoEv * aUtoAngstrm

		open(unit=800,file='veloBerry.txt',action='write')
		write(800,*)"*******************berry velocities******************************"
		write(800,*)"nQ  = ",size(velo,4)
		write(800,*)"nWfs= ",size(velo,2)
		do qi = 1, size(velo,4)
			write(800,*)"*"
			write(800,*)"*"
			write(800,*)" qi=",qi
			!
			!x component
			vSum = dcmplx(0.0_dp)
			do m = 1, size(velo,3)
				do n = 1, size(velo,2)
					write(800,'(a,i3,a,i3,a,e17.10,a,e17.10,a,e17.10)')  "n=",n," m=",m,&
										" v_x=",dreal(velo(1,n,m,qi)),"+i*",dimag(velo(1,n,m,qi)),"; abs=",cdabs(velo(1,n,m,qi))
					vSum = vSum + velo(1,n,m,qi)
				end do
			end do
			write(800,'(a,e17.10,a,e17.10)') "sum v_x = ",dreal(vSum),"+i*",dimag(vSum)
			!
			!y component
			vSum = dcmplx(0.0_dp)
			do m = 1, size(velo,3)
				do n = 1, size(velo,2)
					write(800,'(a,i3,a,i3,a,e17.10,a,e17.10,a,e17.10)')  "n=",n," m=",m,&
										" v_y=",dreal(velo(2,n,m,qi)),"+i*",dimag(velo(2,n,m,qi)),"; abs=",cdabs(velo(2,n,m,qi))
					vSum = vSum + velo(2,n,m,qi)
				end do
			end do
			write(800,'(a,e17.10,a,e17.10)') "sum v_y = ",dreal(vSum),"+i*",dimag(vSum)
			!
			!z component
			vSum = dcmplx(0.0_dp)
			do m = 1, size(velo,3)
				do n = 1, size(velo,2)
					write(800,'(a,i3,a,i3,a,e17.10,a,e17.10,a,e17.10)')  "n=",n," m=",m,&
										" v_z=",dreal(velo(3,n,m,qi)),"+i*",dimag(velo(3,n,m,qi)),"; abs=",cdabs(velo(3,n,m,qi))
					vSum = vSum + velo(3,n,m,qi)
				end do
			end do
			write(800,'(a,e17.10,a,e17.10)') "sum v_z = ",dreal(vSum),"+i*",dimag(vSum)
		end do

		return
	end subroutine


	subroutine writeVeloEffTB(velo)
		complex(dp),	intent(in)		:: 	velo(:,:,:,:) !veloK(3		, 	nWfs	,	nWfs	,	nK		)
		integer							::	qi, n, m
		complex(dp)						:: vSum

		open(unit=800,file='veloTB.txt',action='write')
		write(800,*)"*******************eff TB velocities******************************"
		write(800,*)"nQ  = ",size(velo,4)
		write(800,*)"nWfs= ",size(velo,2)
		do qi = 1, size(velo,4)
			write(800,*)"*"
			write(800,*)"*"
			write(800,*)" qi=",qi
			!
			!x component
			vSum = dcmplx(0.0_dp)
			do m = 1, size(velo,3)
				do n = 1, size(velo,2)
					write(800,'(a,i3,a,i3,a,e17.10,a,e17.10,a,e17.10)')  "n=",n," m=",m,&
										" v_x=",dreal(velo(1,n,m,qi)),"+i*",dimag(velo(1,n,m,qi)),"; abs=",cdabs(velo(1,n,m,qi))
					vSum = vSum + velo(1,n,m,qi)
				end do
			end do
			write(800,'(a,e17.10,a,e17.10)') "sum v_x = ",dreal(vSum),"+i*",dimag(vSum)
			!
			!y component
			vSum = dcmplx(0.0_dp)
			do m = 1, size(velo,3)
				do n = 1, size(velo,2)
					write(800,'(a,i3,a,i3,a,e17.10,a,e17.10,a,e17.10)')  "n=",n," m=",m,&
										" v_y=",dreal(velo(2,n,m,qi)),"+i*",dimag(velo(2,n,m,qi)),"; abs=",cdabs(velo(2,n,m,qi))
					vSum = vSum + velo(2,n,m,qi)
				end do
			end do
			write(800,'(a,e17.10,a,e17.10)') "sum v_y = ",dreal(vSum),"+i*",dimag(vSum)
			!
			!z component
			vSum = dcmplx(0.0_dp)
			do m = 1, size(velo,3)
				do n = 1, size(velo,2)
					write(800,'(a,i3,a,i3,a,e17.10,a,e17.10,a,e17.10)')  "n=",n," m=",m,&
										" v_z=",dreal(velo(3,n,m,qi)),"+i*",dimag(velo(3,n,m,qi)),"; abs=",cdabs(velo(3,n,m,qi))
					vSum = vSum + velo(3,n,m,qi)
				end do
			end do
			write(800,'(a,e17.10,a,e17.10)') "sum v_z = ",dreal(vSum),"+i*",dimag(vSum)
		end do

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
		open(unit=520,file='rawData/EnTB.dat',form='unformatted',access='stream',action='write')
		write(520)	Ew
		close(520)
		!
		!
		return
	end subroutine


	subroutine writePeierls(ckP, EnP)
		complex(dp),	intent(in)		:: ckP(:,:,:)
		real(dp),		intent(in)		:: EnP(:,:)
		real(dp),		allocatable		:: buffer(:,:)
		integer							:: qi
		!
		allocate(	buffer( size(ckP,1)	,	size(ckP,2)		)	)
		
		!real(UNK)
		
		open(unit=700,file='rawData/ckPeiR.dat',form='unformatted',access='stream',action='write')
		do qi = 1, size(ckP,3)
			buffer	= dreal(ckP(:,:,qi))
			write(700)	buffer
		end do
		close(700)
		!imag(UNK)
		open(unit=705,file='rawData/ckPeiI.dat',form='unformatted',access='stream',action='write')
		do qi = 1, size(ckP,3)
			buffer	= dimag(ckP(:,:,qi))
			write(700)	buffer
		end do
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
		open(unit=800,file='rawData/EnBerry.dat',form='unformatted',access='stream',action='write')
		write(800)	EnH
		close(800)
		!
		return
	end subroutine


	subroutine writePolFile(pWann, pBerry, pNiuF2, pNiuF3, pPei )	!writePolFile(pWann, pBerry, pNiu, pPei )
		real(dp),		intent(in)		:: pWann(2), pBerry(3), pNiuF2(3), pNiuF3(3), pPei(3)
		real(dp)						:: pNiu(3), aUtoConv
		!	
		aUtoConv = 1.602176565_dp / 5.2917721092_dp * 1e-4_dp  ! converts from [a.u.] to  [mücro C / cm]
		!	
		open(unit=600,file='polOutput.txt',action='write')
		write(600,*)"**************POLARIZATION OUTPUT FILE**********************"
		write(600,*)" via wavefunction method"
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
		write(600,*) "aX/vol=",aX/vol,"aY/vol=",aY/vol

	
		write(600,'(a,f16.7,a,f16.7,a,a,f16.7,a,f16.7,a)')	"pWann =  (",  pWann(1)	,	", ",	pWann(2),		") [a.u.],",& 
												" moded=(",dmod(pWann(1),aX/vol)*aUtoConv,", ",dmod(pWann(2),aY/vol)*aUtoConv,") [muC/cm]."
		!
		write(600,'(a,f16.7,a,f16.7,a,f16.7,a,a,f16.7,a,f16.7,a)')	"pBerry=  (",		pBerry(1)	,	", ",&	
																					pBerry(2), " ,", pBerry(3)	,	") [a.u.],",& 
												" moded=(",dmod(pBerry(1),aX/vol)*aUtoConv,", ",dmod(pBerry(2),aY/vol)*aUtoConv,") [muC/cm]."
		!write(600,'(a,e16.9,a,f16.12,a,f16.12,a)')	"pInt= ",norm2(pInt)," * (", &	
		!														pInt(1)/norm2(pInt),	", ",	pInt(2)/norm2(pInt),	")"
		!write(600,'(a,e16.9,a,f16.12,a,f16.12,a)')	"pIon= ",norm2(pIon)," * (", &	
		!														pIon(1)/norm2(pIon)	,	", ",	pIon(2)/norm2(pIon),	")"
		!write(600,'(a,e16.9,a,f16.12,a,f16.12,a)')	"pTot= ",norm2(pTot)," * (", &	
		!														pTot(1)/norm2(pTot),	", ",	pTot(2)/norm2(pTot),	")"
		!
		!
		write(600,*)"**************PERTURBATION:"
		write(600,*) "states considered for perturbation nStates=",nSolve
		if( norm2(Bext) > machineP ) then
			write(600,'(a,f16.12,a,f16.12,a,f16.12,a,f16.12,a)')	"Bext= ",norm2(Bext) ," * (", &
											Bext(1)/norm2(Bext),	", ",	Bext(2)/norm2(Bext),", ", Bext(3)/norm2(Bext),	")"
		else
			write(600,'(a,f16.8,a,f16.8,a,f16.8,a)')	"Bext= (", 	Bext(1),	", ",	Bext(2),", ", Bext(3),	")"
		end if
		write(600,*)"*"
		

		write(600,*)"**************FIRST ORDER POL:"
		!NIU
		write(600,'(a,f16.8)') "F3 prefactor = ",prefactF3
		write(600,'(a,e16.7,a,e16.7,a,e16.7,a,a,e16.7,a,e16.7,a)')	"pNiuF2= (", 	pNiuF2(1),	", ",	pNiuF2(2),", ", pNiuF2(3),	")[a.u.]",&
														" moded=(",dmod(pNiuF2(1),aX/vol)*aUtoConv,", ",dmod(pNiuF2(2),aY/vol)*aUtoConv,") [muC/cm]."
		!
		write(600,'(a,e16.7,a,e16.7,a,e16.7,a,a,e16.7,a,e16.7,a)')	"pNiuF3= (", 	pNiuF3(1),	", ",	pNiuF3(2),", ", pNiuF3(3),	")[a.u.]",&
														" moded=(",dmod(pNiuF3(1),aX/vol)*aUtoConv,", ",dmod(pNiuF3(2),aY/vol)*aUtoConv,") [muC/cm]."
		!
		pNiu(:)	= pNiuF2(:) + pNiuF3(:)
		!												
		write(600,'(a,e16.7,a,e16.7,a,e16.7,a,a,e16.7,a,e16.7,a)')	"pNiu  = (", 	pNiu(1),	", ",	pNiu(2),", ", pNiu(3),	")[a.u.]",&
														" moded=(",dmod(pNiu(1),aX/vol)*aUtoConv,", ",dmod(pNiu(2),aY/vol)*aUtoConv,") [muC/cm]."
		
		!PEIERLS													
		write(600,'(a,e16.7,a,e16.7,a,e16.7,a,a,e16.7,a,e16.7,a)')	"pPei  = (", 	pPei(1),	", ",	pPei(2),", ", pPei(3),	")[a.u.]",&
																" moded=(",dmod(pPei(1),aX/vol)*aUtoConv,", ",dmod(pPei(2),aY/vol)*aUtoConv,") [muC/cm]."
		close(600)
		!
		!
		return
	end subroutine



	subroutine printBasisInfo()
		!nG to Gcut RATIO
		if(			nG		< 		vol * Gcut * dsqrt(Gcut) /	(2.0_dp*PI_dp**2)		) then
			write(*,*)	"[BasisInfo]:	vol: increase Gcut"
		else
			write(*,*)	"[BasisInfo]: 	vol: o.k. Gcut"
		end if
		!nR to Gcut RATIO
		if(		nR	<		vol * dsqrt(Gcut) 	/ PI_dp) then
			write(*,*)	"[BasisInfo]: nR: increase real space grid (or decrease Gcut)"
		else
			write(*,*)	"[BasisInfo]: nR: o.k."
		end if
		!real space grid TO supercells
		if(		nRx / nSCx <  10)  then
			write(*,*)	"[BasisInfo]: nR per cell: need more grid points in each cell"
		else	
			write(*,*)	"[BasisInfo]: nR per cell: o.k."
		end if
		!
		!
		return 
	end subroutine



	subroutine printTiming(aT,kT,wT,pwT,bT,oT,mastT)
		real,	intent(in)	:: aT,kT,wT,pwT,bT,oT,mastT
		!
		print '    ("r&alloc  time spend     = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									aT 				, 100.0_dp*aT		   	/mastT
		print '    ("k-solver time spend     = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									kT 				, 100.0_dp*kT		   	/mastT
		print '    ("prep w90 time spend     = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									wT				, 100.0_dp*wT		   	/mastT
		print '    ("post w90 - eff TB       = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									pwT				, 100.0_dp*pwT		   	/mastT
		print '    ("berry method            = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									bT				, 100.0_dp*bT		   	/mastT														
		print '    ("writing  time spend     = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									oT 				, 100.0_dp*oT		   	/mastT
		print '    ("other    time spend     = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									(mastT-aT-kT-wT-bT-oT) 	, 100.0_dp*(mastT-aT-kT-wT-bT-oT)/mastT
		print '    ("overall  time spend     = ",f16.4," seconds.")', mastT
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
