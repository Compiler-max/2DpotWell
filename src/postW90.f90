module postW90
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP,  & 
								aUtoAngstrm, aUtoEv, & 
								myExp, myLeviCivita, eigSolver
	use sysPara
	use wannInterp,		only:	DoWannInterpol
	use polarization,	only:	calcPolWannCent, calcPolViaA
	use semiclassics,	only:	calcFirstOrdP
	use output,			only:	writeInterpBands, writeVeloEffTB

	implicit none

	private
	public	::					effTBmodel
	

	character(len=3)				:: 	seed_name
	integer							:: 	num_kpts, num_wann, nrpts, R_null
	integer,		allocatable		:: 	wigStzDegn(:), R_vect(:,:)
	real(dp)						::	recip_latt(3,3)
	real(dp),		allocatable		:: 	krel(:,:), En_vec(:,:), R_real(:,:), wCent(:,:)
	complex(dp),	allocatable		::	H_tb(:,:,:), r_tb(:,:,:,:), &
										A_mat(:,:,:,:), Om_tens(:,:,:,:,:),	H_mat(:,:,:), Ha_mat(:,:,:,:), &
										U_mat(:,:,:), Om_mat(:,:,:,:), v_mat(:,:,:,:)


	contains

	

!public
	subroutine effTBmodel()
		real(dp)					:: pWann(3), pConn(3), pNiuF2(3), pNiuF3(3), pPei(3)
		complex(dp), allocatable	:: dummy(:,:,:)



		pWann 	= 0.0_dp
		pConn	= 0.0_dp
		pNiuF2	= 0.0_dp
		pNiuF3	= 0.0_dp 
		pPei	= 0.0_dp
		!read seed_name from eStructure input file
		seed_name= seedName
		!If TB file found do calc
		if( readTBsingle() ) then
			write(*,*)	"[effTBmodel]: done reading eff tb matrices"
			!OWN IMPL
			call wannInterpolator()
			!BERRY IMPL
			!allocate(	dummy(				num_wann, 	num_wann ,  nQ	)	)
			!allocate(	A_mat(		3,		num_wann,	num_wann,	nK	)	)		
			!allocate(	En_vec(						num_wann	,	nK	)	)
			!allocate(	v_mat(		3,		num_wann,	num_wann,	nK	)	)
			!allocate(	Om_mat(		3,		num_wann,	num_wann,	nK	)	)
			!dummy = dcmplx(0.0_dp)	!need to set doVeloNum = true
			!call DoWannInterpol( dummy, r_tb, H_tb, R_real, En_vec, A_mat, Om_mat, v_mat)
			write(*,*)	"[effTBmodel]: done interpolating to k mesh with nK=",nK
			if( pw90GaugeB ) then
				call gaugeTrafo()
				write(*,*)	"[effTBmodel]: done gauging back to (H) gauge"
			end if
			!calc all desired polarizations
			call polWrapper(pWann, pConn, pNiuF2, pNiuF3, pPei)
		else
			write(*,*)	"[effTBmodel]: did not find input file, no calculations performed"
		end if
		!
		!output file
		call writePw90pol( pWann, pConn, pNiuF2, pNiuF3, pPei)
		call writeVeloEffTB(v_mat)
		if( writeBin ) call writeInterpBands(En_vec)
		!
		!
		return
	end subroutine










!privat
	logical function readTBsingle( )
		integer						:: 	stat, cnt, offset, R, n, m, i, mn(2), dumI(3), line15(15)
		real(dp)					::	real2(2), real6(6), real3(3)
		!try opening file
		open(unit=310, iostat=stat, file=seed_name//'_tb.dat', status='old', action='read' )
		if( stat /= 0)  then
			write(*,*) "[readTBsingle]: warning, file _tb.dat not found"
			readTBsingle 	= .false.
			recip_latt		= 0.0_dp
		else
			readTBsingle	= .true.
			!
			read(310,*)
			!recip lattice (read into buffer, avoids compiler warning)
			read(310,*) 		real3(:)
			recip_latt(1,:)	= 	real3(:)
			read(310,*)			real3(:)
			recip_latt(2,:)	= 	real3(:)
			read(310,*)			real3(:)
			recip_latt(3,:)	= 	real3(:)
			!sys info
			read(310,*) num_wann
			read(310,*)	nrpts
			!
			allocate( 	wigStzDegn(									nrpts	)		)
			allocate(	R_vect(			3,							nrpts	)		)
			allocate(	R_real(			3,							nrpts	)		)
			allocate(	H_tb(				num_wann,	num_wann, 	nrpts	)		)
			allocate(	r_tb(			3,	num_wann,	num_wann,	nrpts	)		)
			allocate(	wCent(			3,		num_wann					)		)	
			!
			!read degeneracy of each wigner seitz grid point
			cnt 	= 0
			offset 	= 0
			if(	nrpts <= 15 ) then
				read(310,*) wigStzDegn
			else
				do while ( cnt < nrpts )		!read 15 entries per line till nrpts real2ues are read 
					if( nrpts - cnt >= 15	) then
						read(310,*)		line15
						cnt	= cnt + 15
						do i = 1 , 15
							if(offset+i <= size(wigStzDegn))	wigStzDegn(offset+i)	= line15(i) 
						end do
						offset=	offset + 15
					else						!the last line might contain less then 15 entries
						read(310,*)	wigStzDegn( (offset+1):(offset+(nrpts-cnt)) )
						cnt = cnt + (nrpts-cnt)
						offset = offset + (nrpts-cnt)
					end if
				end do
			end if
			!
			!READ HOPPINGS
			do R = 1, nrpts
				!skip first line
				read(310,*)
				!second line fractional real6 of R
				read(310,*)	dumI(1:3)
				R_vect(1:3,R)	= dumI(1:3)
				if( R_vect(1,R)==0 .and. R_vect(2,R)==0 .and. R_vect(3,R)==0 ) R_null = R
				!get indices & hoppings
				do n = 1, num_wann
					do m = 1, num_wann
						read(310,*) mn(1:2), real2(1:2)
						H_tb(mn(1),mn(2),R)	= dcmplx(real2(1)) + i_dp * dcmplx(real2(2))
					end do
				end do
			end do
			!
			!READ POSITIONS
			do R = 1, nrpts
				!skip first line
				read(310,*)
				!second line fractional real6 of R
				read(310,*)	dumI(1:3)
				if( dumI(1) /= R_vect(1,R) .or. dumI(2) /= R_vect(2,R) .or. dumI(3) /= R_vect(3,R) ) then
					write(*,*)	"[readTBsingle]: warning, detected R vector ordering issue while reading positions"
				end if
				!
				do n = 1, num_wann
					do m = 1, num_wann
						read(310,*) mn(1:2), real6(1:6)
						r_tb(1,mn(1),mn(2),R)	= dcmplx(real6(1)) + i_dp * dcmplx(real6(2))
						r_tb(2,mn(1),mn(2),R)	= dcmplx(real6(3)) + i_dp * dcmplx(real6(4))
						r_tb(3,mn(1),mn(2),R)	= dcmplx(real6(5)) + i_dp * dcmplx(real6(6))
					end do
				end do
			end do
			!
			!CONVERT BACK TO [a.u.]
			H_tb	= H_tb /	aUtoEv 
			r_tb	= r_tb /	aUtoAngstrm
			!
			!get centers
			do n = 1, num_wann
				wCent(1:3,n)	= r_tb(1:3,n,n,R_null)
			end do
			!
			!calculate real R vector
			do R = 1 , nrpts
				R_real(1,R)	= R_vect(1,R)	* aX
				R_real(2,R)	= R_vect(2,R)	* aY
				R_real(3,R)	= R_vect(3,R)	* 0
				if( abs(R_real(1,R)-Rcell(1,R)) > machineP ) then
					write(*,*) "[readTB]: warning  Rcell and R_real dont match (x comp)" 
					write(*,*) "			R_real(x)=",R_real(1,R)," Rcell=",Rcell(1,R)
					write(*,*) "			R_real(y)=",R_real(2,R)," Rcell=",Rcell(2,R)
				end if
				!if( abs(R_real(2,R)-Rcell(2,R)) > machineP ) write(*,*) "[readTB]: warning Rcell and R_real dont match(y comp)"
			end do
			!
		end if
		!
		!
		return
	end function


	subroutine wannInterpolator()
		integer						:: ki, R, a, b, c
		complex(dp)					:: phase
		!
		allocate(	A_mat(		3,		num_wann,	num_wann,	nK	)	)
		allocate(	Om_tens(	3,	3,	num_wann,	num_wann,	nK	)	)
		allocate(	H_mat(				num_wann,	num_wann,	nK	)	)
		allocate(	Ha_mat(		3,		num_wann,	num_wann,	nK	)	)
		!
		allocate(	En_vec(						num_wann	,	nK	)	)
		allocate(	U_mat(				num_wann,	num_wann,	nK	)	)
		allocate(	v_mat(		3,		num_wann,	num_wann,	nK	)	)
		allocate(	Om_mat(		3,		num_wann,	num_wann,	nK	)	)
		!
		A_mat	= dcmplx(0.0_dp)
		Om_tens = dcmplx(0.0_dp)
		H_mat	= dcmplx(0.0_dp)
		Ha_mat	= dcmplx(0.0_dp)
		En_vec	= 0.0_dp
		v_mat	= dcmplx(0.0_dp)
		!
		!SET UP K SPACE MATRICES
		do ki = 1 , nK
			do R = 1, nrpts
				phase		= myExp( 	dot_product(kpts(1:2,ki),R_real(1:2,R))		) !/ dcmplx(real(nrpts,dp))
				!
				H_mat(:,:,ki)		= H_mat(:,:,ki)	 		+ phase 								* H_tb(:,:,R)
				do a = 1, 3
					Ha_mat(a,:,:,ki)	= Ha_mat(a,:,:,ki) 	+ phase * i_dp * dcmplx(R_real(a,R))	* H_tb(:,:,R)
					A_mat(a,:,:,ki)		= A_mat(a,:,:,ki)	+ phase									* r_tb(a,:,:,R)
					!
					do b = 1, 3
						Om_tens(a,b,:,:,ki)	= Om_tens(a,b,:,:,ki)	+  phase * i_dp * dcmplx(R_real(a,R)) * r_tb(b,:,:,R)
						Om_tens(a,b,:,:,ki)	= Om_tens(a,b,:,:,ki)	-  phase * i_dp * dcmplx(R_real(b,R)) * r_tb(a,:,:,R)
					end do
				end do
			end do
		end do
		!ENERGY INTERPOLATION
		do ki = 1, nK
			U_mat(:,:,ki)	= H_mat(:,:,ki)
			call eigSolver(U_mat(:,:,ki),	En_vec(:,ki))
		end do
	
		!
		!CURVATURE TO MATRIX
		do ki = 1, nK
			do c = 1, 3
				do b = 1, 3
					do a = 1,3
						Om_mat(c,:,:,ki)	= myLeviCivita(a,b,c) * Om_tens(a,b,:,:,ki)
					end do
				end do
			end do
		end do
		!
		!VELOCITIES
		call calcVelo()
		!
		!
		return
	end subroutine


	subroutine calcVelo()
		integer							:: ki, m, n, i
		complex(dp),	allocatable		:: Hbar(:,:,:), Abar(:,:,:), Ucjg(:,:), U(:,:), tmp(:,:)
		!
		allocate(		Hbar(	size(Ha_mat,1),	size(Ha_mat,2),	size(Ha_mat,3)		)	)		
		allocate(		Abar(	size(A_mat ,1),	size(A_mat ,2),	size(A_mat ,3)		)	)
		allocate(		U(	size(U_mat,1),	size(U_mat,2)							)	)
		allocate(		Ucjg(	size(U_mat,1),	size(U_mat,2)						)	)
		allocate(		tmp(	size(U_mat,1),	size(U_mat,2)						)	)
		!
		do ki = 1, nK
			!GAUGE BACK
			U				= U_mat(:,:,ki)
			Ucjg			= dconjg(	transpose(U)	)
			do i = 1, 3
				!ROTATE TO HAM GAUGE
				tmp			= matmul(	Ha_mat(i,:,:,ki)	, Ucjg			)	
				Hbar(i,:,:)	= matmul(	U				, tmp				)	
				!
				tmp			= matmul(	A_mat(i,:,:,ki)		, Ucjg			)	
				Abar(i,:,:)	= matmul(	U				, tmp				)
				!APPLY ROTATION
				do m = 1, num_wann
					do n = 1, num_wann
						if( n==m )	v_mat(i,n,n,ki) = Hbar(i,n,n)
						if( n/=m )	v_mat(i,n,m,ki) = - i_dp * dcmplx( En_vec(m,ki) - En_vec(n,ki) ) * Abar(i,n,m) 
						!v_mat(1:3,n,m,ki)	=  Ha_mat(1:3,n,m,ki)	- i_dp * dcmplx( En_vec(m,ki) - En_vec(n,ki) ) * A_mat(1:3,n,m,ki) 
						!DEBUG
						if( n/=m .and. abs(Hbar(i,n,m)) > 0.1_dp ) then
							write(*,'(a,i1,a,i3,a,i3,a,f8.4,a,f8.4,a,f8.4)')"[calcVelo]: found off diag band deriv i=",i,&
									" n=",n," m=",m, "v_nm=",dreal(Hbar(i,n,m)), "+i*",dimag(Hbar(i,n,m))," abs=",abs(Abar(i,n,n))
						end if
					end do
				end do
			end do	
			
			!NO GAUGE BACK
			!do m = 1, num_wann
			!	do n = 1, num_wann
			!		if( n==m )	v_mat(1:3,n,n,ki) = Ha_mat(1:3,n,n,ki)
			!		if( n/=m )	v_mat(1:3,n,m,ki) =  - i_dp * dcmplx( En_vec(m,ki) - En_vec(n,ki) ) * A_mat(1:3,n,m,ki) 
			!		!v_mat(1:3,n,m,ki)	=  Ha_mat(1:3,n,m,ki)	- i_dp * dcmplx( En_vec(m,ki) - En_vec(n,ki) ) * A_mat(1:3,n,m,ki) 
			!		!DEBUG
			!		if( n/=m .and. abs(Hbar(i,n,m)) > 0.1_dp ) then
			!				write(*,'(a,i1,a,i3,a,i3,a,f8.4,a,f8.4,a,f8.4)')"[calcVelo]: found off diag band deriv i=",i,&
			!						" n=",n," m=",m, "v_nm=",dreal(Hbar(i,n,m)), "+i*",dimag(Hbar(i,n,m))," abs=",abs(Abar(i,n,n))
			!		end if
			!	end do
			!end do
		end do
		!
		return
	end subroutine




























	subroutine gaugeTrafo()
		integer								:: 	ki, a, b, n, m
		complex(dp),		allocatable		:: 	Ucjg(:,:), tmp(:,:), &
												HaBar(:,:,:), DaH(:,:,:), AaBar(:,:,:), OmBar(:,:,:,:)
		!
		allocate(	Ucjg(			num_wann,	num_wann		)	)
		allocate(	tmp(			num_wann,	num_wann		)	)
		allocate(	HaBar(	3,		num_wann,	num_wann		)	)
		allocate(	DaH(	3,		num_wann,	num_wann		)	)
		allocate(	AaBar(	3,		num_wann,	num_wann		)	)
		allocate(	OmBar(	3,	3,	num_wann,	num_wann		)	)
		!
		
		do ki = 1, nK
			!SET UP THE WORK ARRAYS
			Ucjg				= dconjg(	transpose( U_mat(:,:,ki) )		)
			do a = 1, 3
				!conn
				tmp(:,:)		= matmul(	A_mat(a,:,:,ki) 	, 		Ucjg		)
				AaBar(a,:,:)	= matmul(	U_mat(:,:,ki)		,	tmp(:,:)		)
				!HaBar
				tmp(:,:)		= matmul(	Ha_mat(a,:,:,ki) 	, 		Ucjg		)
				HaBar(a,:,:)	= matmul(	U_mat(:,:,ki)		,	tmp(:,:)		)
				!OmBar
				do b = 1, 3
					tmp(:,:)		= matmul(	Om_tens(a,b,:,:,ki)	,	Ucjg	)		
					OmBar(a,b,:,:)	= matmul(	U_mat(:,:,ki)				,	tmp(:,:)		)
				end do
				!DaH
				do n = 1, num_wann
					do m = 1, num_wann
						if(n/=m)	DaH(a,n,m)	= HaBar(a,n,m) /	(	En_vec(m,ki) - En_vec(n,ki)	)
						if(n==m)	DaH(a,n,m)	= dcmplx(0.0_dp)
					end do
				end do
				
			end do
			!
			!
			!YIELD DESIRED QUANTITIES
			!conn
			A_mat(:,:,:,ki) 	= AaBar + i_dp * DaH
			!velo
			do n = 1, num_wann
				do m = 1, num_wann
					v_mat(1:3,n,m,ki)	= HaBar(1:3,n,m) - i_dp * ( En_vec(m,ki) - En_vec(n,ki) ) * AaBar(1:3,n,m)
					
				end do 
			end do
			!curv
			do a = 1, 3
				do b = 1, 3
					tmp	= dcmplx(0.0_dp)
					tmp	= tmp + OmBar(a,b,:,:)
					tmp	= tmp - 		matmul( DaH(a,:,:), AaBar(b,:,:) ) 	+ 		matmul( AaBar(b,:,:) , DaH(a,:,:)	)
					tmp = tmp + 		matmul( DaH(b,:,:), AaBar(a,:,:) )	-		matmul( AaBar(a,:,:) , DaH(b,:,:)	)
					tmp = tmp - i_dp * 	matmul(	DaH(a,:,:), DaH(b,:,:)	 )  + i_dp *matmul( 	DaH(b,:,:), DaH(a,:,:)	)
					!
					Om_tens(a,b,:,:,ki)	= tmp(:,:)
				end do
			end do
		end do
		!
		return
	end subroutine


	subroutine polWrapper(pWann, pConn, pNiuF2, pNiuF3, pPei)
		real(dp),		intent(out)		:: pWann(3), pConn(3), pNiuF2(3), pNiuF3(3), pPei(3)
		!
		!
		!POLARIZATION CALC
		write(*,*)	"[effTBmodel]: start calculating electric polarization:"
		write(*,'(a,f6.3,a,f6.3,a,f6.3,a)')	"R_real(:,R_null)= (",&
												R_real(1,R_null),", ",R_real(2,R_null),", ",R_real(3,R_null),")."
		!0th order pol
		pWann	= 0.0_dp
		call calcPolWannCent( wCent, pWann )
		write(*,'(a,f10.4,a,f10.4,a,f10.4,a)')	"[effTBmodel]: pWann=(",pWann(1),", ",pWann(2),", ",pWann(3),")."
		call calcPolViaA(A_mat, pConn)
		write(*,'(a,f10.5,a,f10.5,a,f10.5,a)')	"[effTBmodel]: pConn=(",pConn(1),", ",pConn(2),", ",pConn(3),")."
		!1st order pol
		call calcFirstOrdP(Om_mat,A_mat,v_mat,En_vec,pNiuF2, pNiuF3) !calcFirstOrdP(Fcurv, Aconn, Velo, En, p1F2, p1F3)
		write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[effTBmodel]: pNiuF2=(",pNiuF2(1),", ",pNiuF2(2),", ",pNiuF2(3),")."
		write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[effTBmodel]: pNiuF3=(",pNiuF3(1),", ",pNiuF3(2),", ",pNiuF3(3),")."
		
		!1st order peierls
		!if(	doPei ) then
		!	!todo
		!	pPei	= 0.0_dp
		!	write(*,'(a,f10.5,a,f10.5,a,f10.5,a)')	"[effTBmodel]: pPei=(",pPei(1),", ",pPei(2),", ",pPei(3),")."
		!end if
		!
		!
		return
	end subroutine


	subroutine writePw90pol( pWann, pConn, pNiuF2, pNiuF3, pPei)
		real(dp),		intent(in)		::	pWann(3), pConn(3), pNiuF2(3), pNiuF3(3), pPei(3)
		real(dp)						:: 	pNiu(3)
		!
		open(unit=600,file='pW90pol.txt',action='write')
		write(600,*)"**************POST W90 POLARIZATION**********************"
		write(600,*)"via effective tight binding model"
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
		!
		!
		write(600,'(a,f16.7,a,f16.7,a,f16.7,a,a,f16.7,a,f16.7,a)')	"pWann =  (",  pWann(1)	,	", ",	pWann(2),	", ",pWann(3),	"),",& 
												" moded=(",dmod(pWann(1),aX/vol),", ",dmod(pWann(2),aY/vol),")."
		!
		write(600,'(a,f16.7,a,f16.7,a,f16.7,a,a,f16.7,a,f16.7,a)')	"pBerry=  (",		pConn(1)	,	", ",	pConn(2),	", ",pConn(3),"),",& 
												" moded=(",dmod(pConn(1),aX/vol),", ",dmod(pConn(2),aY/vol),")."
		!write(600,'(a,e16.9,a,f16.12,a,f16.12,a)')	"pInt= ",norm2(pInt)," * (", &	
		!														pInt(1)/norm2(pInt),	", ",	pInt(2)/norm2(pInt),	")"
		!write(600,'(a,e16.9,a,f16.12,a,f16.12,a)')	"pIon= ",norm2(pIon)," * (", &	
		!														pIon(1)/norm2(pIon)	,	", ",	pIon(2)/norm2(pIon),	")"
		!write(600,'(a,e16.9,a,f16.12,a,f16.12,a)')	"pTot= ",norm2(pTot)," * (", &	
		!														pTot(1)/norm2(pTot),	", ",	pTot(2)/norm2(pTot),	")"
		!
		!
		write(600,*)"**************PERTURBATION:"
		if( norm2(Bext) > machineP ) then
			write(600,'(a,f16.12,a,f16.12,a,f16.12,a,f16.12,a)')	"Bext= ",norm2(Bext) ," * (", &
											Bext(1)/norm2(Bext),	", ",	Bext(2)/norm2(Bext),", ", Bext(3)/norm2(Bext),	")"
		else
			write(600,'(a,f16.8,a,f16.8,a,f16.8,a)')	"Bext= (", 	Bext(1),	", ",	Bext(2),", ", Bext(3),	")"
		end if
		write(600,*)"*"
		!			
		!	
		write(600,*)"**************FIRST ORDER POL:"
		!NIU
		write(600,'(a,f16.8)') "F3 prefactor = ",prefactF3
		write(600,'(a,f16.7,a,f16.7,a,f16.7,a,a,f16.7,a,f16.7,a)')	"pNiuF2= (", 	pNiuF2(1),	", ",	pNiuF2(2),", ", pNiuF2(3),	")",&
														" moded=(",dmod(pNiuF2(1),aX/vol),", ",dmod(pNiuF2(2),aY/vol),")."
		!
		write(600,'(a,f16.7,a,f16.7,a,f16.7,a,a,f16.7,a,f16.7,a)')	"pNiuF3= (", 	pNiuF3(1),	", ",	pNiuF3(2),", ", pNiuF3(3),	")",&
														" moded=(",dmod(pNiuF3(1),aX/vol),", ",dmod(pNiuF3(2),aY/vol),")."
		!
		pNiu(:)	= pNiuF2(:) + pNiuF3(:)
		!												
		write(600,'(a,f16.7,a,f16.7,a,f16.7,a,a,f16.7,a,f16.7,a)')	"pNiu  = (", 	pNiu(1),	", ",	pNiu(2),", ", pNiu(3),	")",&
														" moded=(",dmod(pNiu(1),aX/vol),", ",dmod(pNiu(2),aY/vol),")."
		
		!PEIERLS													
		write(600,'(a,f16.7,a,f16.7,a,f16.7,a,a,f16.7,a,f16.7,a)')	"pPei  = (", 	pPei(1),	", ",	pPei(2),", ", pPei(3),	")",&
																" moded=(",dmod(pPei(1),aX/vol),", ",dmod(pPei(2),aY/vol),")."
		close(600)
		!
		return
	end subroutine







end module


