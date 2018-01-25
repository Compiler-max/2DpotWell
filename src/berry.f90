module berry
	!module contains methods related to the Berry phase theory of modern polarization
	!	i.e. calculation of connection, velocities, curvatures and polarization
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP,  myExp, myLeviCivita, aUtoAngstrm, aUtoEv
	use sysPara
	use planeWave,		only:	calcVeloGrad, calcConnOnCoarse
	use polarization,	only:	calcPolViaA
	use semiClassics,	only:	calcFirstOrdP
	use peierls,		only:	peierlsMethod
	use output,			only:	writePolFile, writeVeloHtxt

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
		real(dp),		allocatable		:: 	R_real(:,:)
		complex(dp),	allocatable		:: 	AconnQ(:,:,:,:), FcurvQ(:,:,:,:),veloQ(:,:,:,:)					
		real(dp)						::	pWann(3)
		real(dp),		allocatable		::	v_Band(:,:,:)
		!					
		!COARSE
		allocate(			ck(			Gmax	,	nSolve				,	nQ		)			)
		allocate(			ckW(		Gmax	, 	nWfs				,  	nQ		)			)
		allocate(			EnQ(		nSolve							,	nQ		)			)
		allocate(			AconnQ(		3		, 	nWfs	,	nWfs	,	nQ		)			)
		allocate(			FcurvQ(		3		,	nWfs	,	nWfs	,	nQ		)			)
		allocate(			veloQ(		3		, 	nSolve	,	nSolve	,	nQ		)			)
		allocate(			R_real(		3		,							nSC		)			)
		allocate(			v_Band(		3		,			nSolve		,	nQ		)			)
		!
		!READ IN QUANTITIES
		call initRead(R_real, ck, EnQ)

		!
		!ROTATE
		call applyRot(ck, Uq, ckW)
		!
		!CONNECTION
		call calcConnOnCoarse(ckW, AconnQ)
		call calcPolViaA(AconnQ,pBerry)
		write(*,*)"[berryMethod]: coarse rotated pol =(",pBerry(1),", ",pBerry(2),", ", pBerry(3),")."
		!
		!1st ORDER SEMICLASSICS
		if(doNiu) then
			write(*,*)	"[berrryMethod]: now calc first order pol"
			veloQ	= dcmplx(0.0_dp)
			FcurvQ	= dcmplx(0.0_dp)	!does not matter since <FcurvQ,AconnQ> is always zero in 2D
			call calcVelo(ck , Uq , AconnQ, EnQ ,  veloQ)
			!
			call calcFirstOrdP(FcurvQ, AconnQ, veloQ, EnQ, pNiuF2, pNiuF3)
			write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[berryMethod]: pNiuF2=(",pNiuF2(1),", ",pNiuF2(2),", ",pNiuF2(3),")."
			write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[berryMethod]: pNiuF3=(",pNiuF3(1),", ",pNiuF3(2),", ",pNiuF3(3),")."
		else
			pNiuF2 = 0.0_dp
			pNiuF3 = 0.0_dp
		end if


		!PEIERLS SUBSTITUTION
		if(doPei) then
			write(*,*)	"[berrryMethod]: now calc first order pol via peierls sub."
			call peierlsMethod(ckW, pPei)
			write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[berryMethod]: pPei=(",pPei(1),", ",pPei(2),", ",pPei(3),")."
		else
			pPei = 0.0_dp
		end if

		!WANNIER
		!if(doWanni) then
		!	write(*,*) "[berryMethod]: now calc the Wannier functions by hand"
		!	!call wannMethod(ckW, pWann)
		!	write(*,'(a,f8.4,a,f8.4,a,f8.4,a)')	"[berryMethod]: pWann=(",pWann(1),", ",pWann(2),", ",pWann(3),")."
		!else
		!	pWann	= 0.0_dp
		!end if
		pWann = 0.0_dp

		!OUTPUT
		call writePolFile(pWann, pBerry, pNiuF2, pNiuF3, pPei )
		call writeConnTxt( AconnQ )
		call writeVeloHtxt( veloQ)!*aUtoEv*aUtoAngstrm )				
		!if( writeBin )	call write(ck, ckW)

		write(*,*)	"[berrryMethod]: all done"
		!
		!
		!
		return
	end subroutine




















!private
	subroutine initRead(R_real, ck, EnQ)
		real(dp),		intent(out)		:: 	R_real(:,:), EnQ(:,:)
		complex(dp),	intent(out)		:: 	ck(:,:,:)
		integer							::	R, qi, n 
		!
		!INIT SUPERCELL
		R_real(3,:)	= 0.0_dp
		do R = 1, nSC
			R_real(1:2,R)	= Rcell(1:2,R)
		end do 
		!
		!READ U FROM W90
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
		!
		!READ ABINITIO
		call readHam(ck, EnQ)
		!
		!
		return
	end subroutine



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
				En(1:nSolve,qi)	= eBuff(1:nSolve)
		end do
		close(720)
		!
		!
	
		!
		return
	end subroutine


	subroutine applyRot(ck, Uq, ckW)
		complex(dp),	intent(in)		::	ck(:,:,:), Uq(:,:,:)
		complex(dp),	intent(out)		::	ckw(:,:,:)
		integer							::	qi, gi, n, m
		!
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
				write(*,'(a,a)')	"[berryMethod]: rotations disabled.",&
																" Will use initial electronic structure coeff"
			else 
				ckW	= dcmplx(0.0_dp)
				write(*,*)	"[berryMethod]: critical error, less nBands then nWfs, coeff set to zero..."
			end if
		end if
		!
		!
		return
	end subroutine



	subroutine calcVelo(ck, U_mat , A_mat, En_vec ,  v_mat)
		complex(dp),	intent(in)		::	ck(:,:,:), U_mat(:,:,:), A_mat(:,:,:,:)
		real(dp),		intent(in)		::	En_vec(:,:)
		complex(dp),	intent(out)		::	v_mat(:,:,:,:)
		complex(dp),	allocatable		:: 	Abar(:,:,:), U(:,:), Ucjg(:,:), tmp(:,:)
		real(dp),		allocatable		::	v_Band(:,:,:)
		integer							::	n, m, qi, a

		!
		if( doVeloNum ) then
			!BLOUNT
			write(*,*)	"[beryMethod/calcVelo]: velo via blount formula - warning this is deprecated please set doVeloNum = f"
			allocate(			Abar(		3		,	nWfs	,	nWfs				)			)
			allocate(			tmp(					nWfs	,	nWfs				)			)
			allocate(			U(						nWfs	,	nWfs				)			)
			allocate(			Ucjg(					nWfs	,	nWfs				)			)
			allocate(			v_Band(		3		,			nWfs		,	nQ		)			)
			!
			call readBandVelo( v_Band )
			do qi = 1, nQ
				!(H) GAUGE
				if( doGaugBack ) then
					!GET ROTATED QUANTITIES
					U	 = U_mat(:,:,qi)
					Ucjg = transpose( dconjg(U)	)
					do a = 1, 3
						tmp(:,:)	= matmul(	A_mat(a,:,:,qi) 	,	U	)
						Abar(a,:,:)	= matmul(	Ucjg				,	tmp	)
					end do
					!
					!APPLY
					do m = 1, nWfs
						do n = 1, nWfs
							if(n==m)	v_mat(1:3,n,n,qi)	= v_Band(1:3,n,qi)
							if(n/=m) 	v_mat(1:3,n,m,qi)	= - i_dp * dcmplx( En_vec(m,qi)-En_vec(n,qi) ) * Abar(1:3,n,m)
						end do
					end do
				!(W) GAUGE
				else
					do m = 1, nWfs
						do n = 1, nWfs
							if(n==m)	v_mat(1:3,n,n,qi)	= v_Band(1:3,n,qi)
							!if(n/=m) 	v_mat(1:3,n,m,qi)	= - i_dp * dcmplx( En_vec(m,qi)-En_vec(n,qi) ) * A_mat(1:3,n,m,qi)
							if(n/=m) 	v_mat(1:3,n,m,qi)	= - i_dp * dcmplx( En_vec(m,qi)-En_vec(n,qi) ) * A_mat(1:3,n,m,qi)
							!if(n/=m)	v_mat(1:3,n,m,qi)	= - i_dp * dcmplx( En_vec(m,qi)-En_vec(n,qi) ) * Abar(1:3,n,m)
						end do
					end do
				end if
			end do
		else
			!PLANE WAVE GRADIENT
			write(*,*)	"[beryMethod/calcVelo]: velo via plane wave gradients"
			call calcVeloGrad( ck, v_mat)
		end if


		return
	end subroutine

	subroutine readBandVelo( v_vec )
		real(dp),		intent(out)		::	v_vec(:,:,:)
		real(dp)						::	buffer(7)
		integer							::	stat, qi, n, qInd
		!
		write(*,*)  seed_name//'_geninterp.dat'
		open(unit=320,iostat=stat, file=seedName//'_geninterp.dat',form='formatted', status='old',action='read')
		if( stat/= 0 ) then 
			v_vec = 0.0_dp
			write(*,*)	"[readBandVelo]: could not find _geninterp.dat file.. velocities set to zero"
		else 
			read(320,*)
			read(320,*)
			read(320,*)
			do qi = 1, nQ
				do n = 1, nWfs
					read(320,*)	qInd, buffer
					v_vec(1,n,qInd)	= buffer(5)
					v_vec(2,n,qInd)	= buffer(6)
					v_vec(3,n,qInd)	= buffer(7) 
				end do
			end do
			close(320)
		end if
		!
		!ATOMIC UNITS CONVERSION:
		v_vec 	= v_vec  / (aUtoEv  * aUtoAngstrm)	 ! [v_vec] = eV / Angstroem
		!
		!	
		return
	end subroutine


	subroutine writeConnTxt( A_mat )
		complex(dp),	intent(in)		::	A_mat(:,:,:,:)
		integer							::	qi, n, m
		!
		open(unit=350,file='AconnBerry.txt',action='write', status='replace')
		write(350,*)	"connection calculated via berryMethod"
		write(350,*)	"n m real(A_x) imag(A_x) real(A_y) imag(A_y) real(A_z) imag(A_z)"
		do qi = 1, size(A_mat,4)
			write(350,*)	"qi=",	qi
			do m = 1, size(A_mat,3)
				do n = 1, size( A_mat,2)
					write(350,'(i3,a,i3,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4)')	n," ",m,&	
													"   ",dreal(A_mat(1,n,m,qi))," ",dimag(A_mat(1,n,m,qi)),&
													"   ",dreal(A_mat(2,n,m,qi))," ",dimag(A_mat(2,n,m,qi)),&
													"   ",dreal(A_mat(3,n,m,qi))," ",dimag(A_mat(3,n,m,qi))
				end do
			end do
		end do
		close(350)
		!
		return
	end subroutine


end module berry