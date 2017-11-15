module berry
	!module contains methods related to the Berry phase theory of modern polarization
	!	i.e. calculation of connection, velocities, curvatures and polarization
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, myExp, myLeviCivita, nIntegrate
	use sysPara
	use input,			only:	readHam
	use effTB,			only:	TBviaKspace
	use wannInterp,		only:	DoWannInterpol
	use polarization,	only:	calcPolViaA
	use semiClassics,	only:	calcFirstOrdP
	use peierls,		only:	peierlsMethod
	use wannier,		only:	wannMethod
	use output,			only:	writeCkASunk, writePolFile
	implicit none

	private
	public	::	berryMethod



	integer								::	num_wann, num_kpts
	character(len=3)					::	seed_name
	complex(dp),	allocatable			::	ck(:,:,:), ckW(:,:,:), Uq(:,:,:)
	real(dp),		allocatable			:: 	EnQ(:,:), krel(:,:)








	contains




!public



	subroutine berryMethod(pBerry, pNiuF2, pNiuF3, pPei)
		!todo
		real(dp),		intent(out)		:: 	pBerry(2), pNiuF2(3), pNiuF3(3), pPei(3)
		real(dp),		allocatable		:: 	EnK(:,:)
		complex(dp),	allocatable		:: 	AconnK(:,:,:,:), FcurvK(:,:,:,:), veloK(:,:,:,:), &
											tHopp(:,:,:), rHopp(:,:,:,:) 
		real(dp)						::	pWann(3)
		integer							::	gi, qi
		!					
		!
		write(*,*)	"[berrryMethod]: hello from Berry"
		if( setBasis() /= 0 ) write(*,*)	"[berryMethod]: could not read nGq file"

		
		write(*,*)	"[berrryMethod]: wrote the U matrix from file"

		allocate(			tHopp(					nWfs	, 	nWfs	,	nSc		)			)
		allocate(			rHopp(		2		,	nWfs	, 	nWfs	, 	nSC		)			)			
		allocate(			EnK(					nWfs	, 				nK		)			)
		allocate(			AconnK(		3		, 	nWfs	,	nWfs	,	nK		)			)
		allocate(			FcurvK(		3		, 	nWfs	, 	nWfs	,	nK		)			)
		allocate(			veloK(		3		, 	nWfs	,	nWfs	,	nK		)			)

		allocate(	ck(		nG,		nBands,	nQ		)	)
		allocate(	ckW(	nG, 	nWfs,  	nQ		)	)
		allocate(	EnQ(	nG,				nQ		)	)

		!read in U matrix (yields nQ, nWfs) 
		call readUmatrix()

		!read in ck, En
		call readHam(ck, EnQ)

		!rotate ck, get ckW
		do qi = 1, nQ
			do gi = 1, nGq(qi)
				ckW(gi,:,qi)	= matmul( ck(gi,:,qi), Uq(:,:,qi))
			end do
		end do

		!SET UP EFFECTIVE TIGHT BINDING MODELL
		call TBviaKspace(ckW, EnQ, Uq, tHopp, rHopp)
		write(*,*)	"[berryMethod]: set up effective tight binding model (k-Space method)"

		!INTERPOLATE CONN,CURV, VELO
		call DoWannInterpol(ckW, rHopp, tHopp, EnK, AconnK, FcurvK, veloK)
		write(*,*)	"[berrryMethod]: interpolation done"

		
		
		!INTEGRATE CONNECTION
		call calcPolViaA(AconnK,pBerry)
		write(*,'(a,f12.8,a,f12.8,a)')	"[berrryMethod]: calculated zero order pol=(",pBerry(1),", ",pBerry(2),")."
		



		!1st ORDER SEMICLASSICS
		if(doNiu) then
			write(*,*)	"[berrryMethod]: now calc first order pol"
			call calcFirstOrdP(FcurvK, AconnK, veloK, EnK, pNiuF2, pNiuF3)
			write(*,'(a,f9.4,a,f9.4,a,f9.4,a)')	"[berryMethod]: pNiuF2=(",pNiuF2(1),", ",pNiuF2(2),", ",pNiuF2(3),")."
			write(*,'(a,f9.4,a,f9.4,a,f9.4,a)')	"[berryMethod]: pNiuF3=(",pNiuF3(1),", ",pNiuF3(2),", ",pNiuF3(3),")."
		end if


		!PEIERLS SUBSTITUTION
		if(doPei) then
			write(*,*)	"[berrryMethod]: now calc first order pol via peierls sub."
			call peierlsMethod(ckW, tHopp, pPei)
			write(*,'(a,f8.4,a,f8.4,a,f8.4,a)')	"[berryMethod]: pPei=(",pPei(1),", ",pPei(2),", ",pPei(3),")."
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

		if( writeBin )	call writeCkASunk(ck, ckW)
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
		end do
		close(300)

		return
	end subroutine



end module berry