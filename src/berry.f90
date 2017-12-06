module berry
	!module contains methods related to the Berry phase theory of modern polarization
	!	i.e. calculation of connection, velocities, curvatures and polarization
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP,  myExp, myLeviCivita, nIntegrate
	use sysPara
	use effTB,			only:	TBviaKspace
	use wannInterp,		only:	DoWannInterpol
	use polarization,	only:	calcPolViaA
	use semiClassics,	only:	calcFirstOrdP
	use peierls,		only:	peierlsMethod
	use wannier,		only:	wannMethod
	use output,			only:	writeCkASunk, writePolFile, writeVeloHtxt, writeEnH, writeUmat, writeBerryInterpU, &
								writeHtbBerry, writeRtbBerry, writeEnBerry

	implicit none

	private
	public	::	berryMethod



	integer								::	num_wann, num_kpts
	character(len=3)					::	seed_name
	complex(dp),	allocatable			::	ck(:,:,:), ckW(:,:,:), Uq(:,:,:)
	real(dp),		allocatable			:: 	EnQ(:,:), krel(:,:)








	contains




!public



	subroutine berryMethod()
		!todo
		real(dp)						:: 	pBerry(3), pNiuF2(3), pNiuF3(3), pPei(3)
		real(dp),		allocatable		:: 	EnK(:,:), R_real(:,:)
		complex(dp),	allocatable		:: 	AconnK(:,:,:,:), FcurvK(:,:,:,:), veloK(:,:,:,:), &
											tHopp(:,:,:), rHopp(:,:,:,:) , Uk(:,:,:)
		real(dp)						::	pWann(3)
		integer							::	gi, qi, n, m, R
		!					
		!
		allocate(			tHopp(					nWfs	, 	nWfs	,	nSc		)			)
		allocate(			rHopp(		3		,	nWfs	, 	nWfs	, 	nSC		)			)			
		allocate(			EnK(					nWfs	, 				nK		)			)
		allocate(			AconnK(		3		, 	nWfs	,	nWfs	,	nK		)			)
		allocate(			FcurvK(		3		, 	nWfs	, 	nWfs	,	nK		)			)
		allocate(			veloK(		3		, 	nWfs	,	nWfs	,	nK		)			)
		allocate(			Uk(						nWfs	,	nWfs	,	nK		)			)

		allocate(			R_real(		3		,							nSC		)			)
		allocate(	ck(		nG,		nBands,	nQ		)	)
		allocate(	ckW(	nG, 	nWfs,  	nQ		)	)
		allocate(	EnQ(	nG,				nQ		)	)

		R_real(3,:)	= 0.0_dp
		do R = 1, nSC
			R_real(1:2,R)	= Rcell(1:2,R)
		end do 

		!read in U matrix (yields nQ, nWfs) 
		if( useRot )  then
			call readUmatrix()
			write(*,*)	"[berryMethod]: read U matrix "
		else
			allocate( Uq( nWfs, nWfs, nQ )	)
			Uq = dcmplx(0.0_dp)
			num_wann	= nWfs
			num_kpts	= nQ
			do qi = 1, nQ
				do n = 1, nWfs
					Uq(n,n,qi)	= dcmplx(1.0_dp)
				end do
			end do
			write(*,*)	"[berryMethod]: U matrix set as Identity"
		end if

		!read in ck, En
		call readHam(ck, EnQ)

		!rotate ck, get ckW
		ckW	= dcmplx(0.0_dp)
		if(	useRot ) then
			do qi = 1, nQ
				do gi = 1, nGq(qi)			! u^(H) = u^(W) U -> u^(W) = u^(H) U^dagger
					do n = 1, num_wann
						!SUM OVER m
						do m = 1, num_wann
							ckW(gi,n,qi)	=  ckW(gi,n,qi) + Uq(m,n,qi)   * ck(gi,m,qi)	
						end do		
					end do	
				end do
			end do
			write(*,*)	"[berryMethod]: applied U matrix to basis coefficients"
		else 
			if( nWfs <= nBands) then
				ckW(:,:,:)	= ck(:,1:nWfs,:)
				write(*,*)	"[berryMethod]: rotations disabled. Will use initial electronic structure coeff"
			else 
				ckW	= dcmplx(0.0_dp)
				write(*,*)	"[berryMethod]: critical error, less nBands then nWfs, coeff set to zero..."
			end if
		end if

		!SET UP EFFECTIVE TIGHT BINDING MODELL
		call TBviaKspace(ckW, EnQ, Uq, tHopp, rHopp)
		call writeHtbBerry(tHopp)
		call writeRtbBerry(rHopp)
		write(*,*)	"[berryMethod]: set up effective tight binding model (k-Space method)"

		!INTERPOLATE CONN,CURV, VELO
		call DoWannInterpol(ckW, rHopp, tHopp, R_real, EnK, Uk, AconnK, FcurvK, veloK)
		call writeEnBerry(EnK)
		call writeBerryInterpU(Uk)
		write(*,*)	"[berrryMethod]: interpolation done"


		
		
		!INTEGRATE CONNECTION
		call calcPolViaA(AconnK,pBerry)
		write(*,'(a,f12.8,a,f12.8,a)')	"[berrryMethod]: calculated zero order pol=(",pBerry(1),", ",pBerry(2),")."
		



		!1st ORDER SEMICLASSICS
		if(doNiu) then
			write(*,*)	"[berrryMethod]: now calc first order pol"
			call calcFirstOrdP(FcurvK, AconnK, veloK, EnK, pNiuF2, pNiuF3)
			write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[berryMethod]: pNiuF2=(",pNiuF2(1),", ",pNiuF2(2),", ",pNiuF2(3),")."
			write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[berryMethod]: pNiuF3=(",pNiuF3(1),", ",pNiuF3(2),", ",pNiuF3(3),")."
		end if


		!PEIERLS SUBSTITUTION
		if(doPei) then
			write(*,*)	"[berrryMethod]: now calc first order pol via peierls sub."
			call peierlsMethod(ckW, tHopp, pPei)
			write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[berryMethod]: pPei=(",pPei(1),", ",pPei(2),", ",pPei(3),")."
		end if

		!WANNIER
		if(doWanni) then
			write(*,*) "[berryMethod]: now calc the Wannier functions by hand"
			call wannMethod(ckW, pWann)
			write(*,'(a,f8.4,a,f8.4,a,f8.4,a)')	"[berryMethod]: pWann=(",pWann(1),", ",pWann(2),", ",pWann(3),")."
		end if

		!OUTPUT
		pWann = 0.0_dp
		call writePolFile(pWann, pBerry, pNiuF2, pNiuF3, pPei )
		call writeVeloHtxt(veloK)
		call writeUmat(Uq)
		if( writeBin )	call writeCkASunk(ck, ckW)
		if( writeBin )	call writeEnH(EnK)

		write(*,*)	"[berrryMethod]: all done"
		!
		!
		!
		return
	end subroutine


!private
	subroutine readUmatrix()
		integer						:: stat, qi, n, m, dumI(3)
		real(dp)					:: val(2)
		!todo
		seed_name	= seedName
		!READ U MATRIX
		open(unit=300,iostat=stat, file=seed_name//'_u.mat', status='old',action='read')
		if( stat /= 0)	write(*,*)	"[readUmatrix]: warning did not file _u.mat file"
		read(300,*)
		read(300,*) dumI(1:3)
		if( dumI(1) /= nQ ) write(*,*)	"[readUmatrix]: warning num_kpts=",dumI(1)," nQ=",nQ
		if(	dumI(2)	/= nWfs) write(*,*)	"[readUmatrix]: warning num_wann=",dumI(2)," nWfs=",nWfs
		num_kpts	= dumI(1)
		num_wann	= dumI(2)
		!
		allocate(	krel(		3,							num_kpts)			)
		allocate(	Uq(		num_wann,	num_wann,	num_kpts)			)	
		!
		do qi = 1,  num_kpts
			read(300,*)
			read(300,*) krel(1:3,qi)
			do n = 1, num_wann
				do m = 1, num_wann
					read(300,*)	val(1:2)
					Uq(m,n,qi)	= dcmplx(val(1))	+	i_dp	*	dcmplx(val(2))
				end do
			end do
			!DEBUG
			if( abs( krel(1,qi)*2.0_dp*PI_dp/aX - qpts(1,qi)) > machineP) then
				write(*,*)	"[readUmatrix]: warning k meshes are ordered diffferently"
				write(*,*)	"				x: k_w90= ",krel(1,qi)*2.0_dp*PI_dp/aX, " qpts=",qpts(1,qi)
				write(*,*)	"				y: k_w90= ",krel(2,qi)*2.0_dp*PI_dp/aY, " qpts=",qpts(2,qi)
			end if
		end do
		close(300)

		return
	end subroutine


	subroutine	readHam(ck, En)
		complex(dp),	intent(out)		:: ck(:,:,:)
		real(dp),		intent(out)		:: En(:,:)
		real(dp),		allocatable		:: buffer(:,:), eBuff(:)
		integer							:: qi
		!
		allocate(	buffer( size(ck,1), size(ck,2) 	)		)
		allocate(	eBuff(size(En,1)				)		)
		!
		!
		!
		!UNK REAL PART
		open(unit=700, file="./rawData/ckR.dat",form='unformatted',access='stream',action='read')
		do qi = 1 , size(ck,3)
			read(700) buffer
			ck(:,:,qi)	= dcmplx(buffer)
		end do
		close(700)
		!
		!UNK IMAG PART
		open(unit=710, file="./rawData/ckI.dat",form='unformatted',access='stream',action='read')
		do qi = 1 , size(ck,3)
			read(710) buffer
			ck(:,:,qi)	= ck(:,:,qi) + i_dp * dcmplx(buffer)
		end do
		close(710)
		!
		!BAND ENERGIES
		open(unit=720, file="./rawData/bandStruct.dat",form='unformatted',access='stream',action='read')
		do qi = 1, size(En,2)	
				read(720) eBuff
				En(1:nBands,qi)	= eBuff(1:nBands)
		end do
		close(720)
		!
		!
	
		!
		return
	end subroutine



end module berry