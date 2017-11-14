module postW90
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP, & 
								aUtoAngstrm, aUtoEv, & 
								myExp, myLeviCivita, eigSolver
	use sysPara
	use polarization,	only:	calcPolWannCent, calcPolViaA
	use semiclassics,	only:	calcFirstOrdP
	use output,			only:	writeInterpBands

	implicit none

	private
	public	::					effTBmodel
	

	character(len=3)				:: 	seed_name
	integer							:: 	num_kpts, num_wann, nrpts, R_null
	integer,		allocatable		:: 	wigStzDegn(:), R_vect(:,:)
	real(dp),		allocatable		:: 	krel(:,:), En_vec(:,:), R_real(:,:), wCent(:,:)
	complex(dp),	allocatable		::	H_tb(:,:,:), r_tb(:,:,:,:), &
										A_mat(:,:,:,:), Om_tens(:,:,:,:,:),	H_mat(:,:,:), Ha_mat(:,:,:,:), &
										U_mat(:,:,:), Om_mat(:,:,:,:), v_mat(:,:,:,:) 


	contains

	subroutine effTBmodel( pWann, pConn, pNiuF2, pNiuF3, pPei)
		real(dp),		intent(out)		:: pWann(3), pConn(3), pNiuF2(3), pNiuF3(3), pPei(3)

		seed_name= seedName
		if(	filesExist() ) then
			!
			!
			!SETUP THE MODEL
			call readTB()
			write(*,*)	"[effTBmodel]: done reading eff tb matrices"
			call wannInterp()
			write(*,*)	"[effTBmodel]: done interpolating to k mesh with nK=",nK
			!call gaugeTrafo()
			!write(*,*)	"[effTBmodel]: done gauging back to (H) gauge"
			!
			!
			!POLARIZATION CALC
			write(*,*)	"[effTBmodel]: start calculating electric polarization:"
			write(*,'(a,f6.3,a,f6.3,a,f6.3,a)')	"R_real(:,R_null)= (",&
													R_real(1,R_null),", ",R_real(2,R_null),", ",R_real(3,R_null),")."
			!0th order pol
			call calcPolWannCent( wCent, pWann )
			write(*,'(a,f10.5,a,f10.5,a,f10.5,a)')	"[effTBmodel]: pWann=(",pWann(1),", ",pWann(2),", ",pWann(3),")."
			call calcPolViaA(A_mat, pConn)
			write(*,'(a,f10.5,a,f10.5,a,f10.5,a)')	"[effTBmodel]: pConn=(",pConn(1),", ",pConn(2),", ",pConn(3),")."
			!1st order pol
			if( doNiu ) then
				call calcFirstOrdP(Om_mat,A_mat,v_mat,En_vec,pNiuF2, pNiuF3) !calcFirstOrdP(Fcurv, Aconn, Velo, En, p1F2, p1F3)
				write(*,'(a,f10.5,a,f10.5,a,f10.5,a)')	"[effTBmodel]: pNiuF2=(",pNiuF2(1),", ",pNiuF2(2),", ",pNiuF2(3),")."
				write(*,'(a,f10.5,a,f10.5,a,f10.5,a)')	"[effTBmodel]: pNiuF3=(",pNiuF3(1),", ",pNiuF3(2),", ",pNiuF3(3),")."
			end if
			!1st order peierls
			if(	doPei ) then
				!todo
				write(*,'(a,f10.5,a,f10.5,a,f10.5,a)')	"[effTBmodel]: pPei=(",pPei(1),", ",pPei(2),", ",pPei(3),")."
			end if
			!
			!output file
			call writePw90pol( pWann, pConn, pNiuF2, pNiuF3, pPei)
			!
		else
			write(*,*)	"[effTBmodel]: did not find all neccessary input files, skip..."
		end if


		return
	end subroutine


	logical function filesExist()
		integer						:: hrExist, rExist
		open(unit=305,iostat=hrExist, file=seed_name//'_hr.dat', status='old',action='read')
		close(305)
		open(unit=306,iostat= rExist, file=seed_name//'_r.dat' , status='old',action='read')
		close(306)
		!	
		filesExist = ( hrExist == 0 ) 	.and.	( rExist == 0 )
		!
		return
	end function







	subroutine readTB()
		integer					:: stat, dumI(3), line15(15), n, m, cnt, offset, i, R, mn(2)
		real(dp)				:: val(2), coord(6)
		!
		!READ H(R) data
		open(unit=310,iostat=stat, file=seed_name//'_hr.dat', status='old',action='read')
		if( stat /= 0)	write(*,*)	"[readTB]: warning did not file _hr.dat file"
		read(310,*)
		read(310,*) num_wann
		read(310,*)	nrpts
		if(nrpts /=  num_kpts	)	write(*,*)	"[readTB]: warning number of unit cells does not match amount of k points"
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
			do while ( cnt < nrpts )		!read 15 entries per line till nrpts values are read 
				if( nrpts - cnt >= 15	) then
					read(310,*)		line15
					cnt	= cnt + 15
					do i = 1 , 15
						if(offset+i <= size(wigStzDegn))	wigStzDegn(offset+i)	= line15(i) 
					end do
					offset=	offset + 15
				else
					read(310,*)	wigStzDegn( (offset+1):(offset+(nrpts-cnt)) )
					cnt = cnt + (nrpts-cnt)
					offset = offset + (nrpts-cnt)
				end if
			end do
		end if
		!read matrix elements
		do R = 1, nrpts
			do n = 1, num_wann
				do m = 1, num_wann
					read(310,*)	dumI(1:3), mn(1:2), val(1:2)
					!get matrix element
					H_tb(mn(1),mn(2),R)	= dcmplx(val(1)) + i_dp * dcmplx(val(2))
					!get R vector( make sure input file has R vectors as major order)
					if( n == 1 .and. m==1 ) then
						R_vect(1:3,R)	= dumI(1:3)
						if( R_vect(1,R)==0 .and. R_vect(2,R)==0 .and. R_vect(3,R)==0 ) R_null = R
					else	
						if( R_vect(1,R)	/= dumI(1) .or. R_vect(2,R)	/= dumI(2) .or. R_vect(3,R)	/= dumI(3)  ) then 
							write(*,*)	"[readTB]: warning R vect. order not given (read _hr.dat)"
						end if
					end if
				end do
			end do
		end do
		close(310)
		!
		!CONVERT UNITS TO A.U.
		H_tb	= H_tb /	aUtoEv 

	

		!READ R MATRIX 
		open(unit=320,iostat=stat, file=seed_name//'_r.dat', status='old',action='read')
		read(320,*)
		read(320,*)	dumI(3)
		if( num_wann /= dumI(3)	) write(*,*)	"[readTB]: warning _r.dat has wrong number of wannier functions"
		do R = 1 , nrpts
			do n = 1, nWfs
				do m = 1, nWfs
					read(320,*)	dumI(1:3), mn(1:2), coord(1:6)
					!get matrix element
					r_tb(1,mn(1),mn(2),R)	= dcmplx(coord(1)) + i_dp * dcmplx(coord(2))
					r_tb(2,mn(1),mn(2),R)	= dcmplx(coord(3)) + i_dp * dcmplx(coord(4))
					r_tb(3,mn(1),mn(2),R)	= dcmplx(coord(5)) + i_dp * dcmplx(coord(6))
					! todo read in all 3 components of r_tb!!!!!!

					!test if R vector ordering is given
					if( R_vect(1,R)	/= dumI(1) .or. R_vect(2,R)	/= dumI(2) .or. R_vect(3,R)	/= dumI(3)  ) then 
						write(*,*)	"[readTB]: warning R vect. order not given (read _r.dat)"
					end if

				end do
			end do
		end do

		!CONVERT UNITS TO ATOMIC UNITS 
		!r_tb	= r_tb / aUtoAngstrm

		!get centers
		do n = 1, num_wann
			wCent(1:3,n)	= r_tb(1:3,n,n,R_null)
		end do
		

		!calculate real R vector
		do R = 1 , nrpts
			R_real(1,R)	= R_vect(1,R)	* aX
			R_real(2,R)	= R_vect(2,R)	* aY
			R_real(3,R)	= R_vect(3,R)	* 0
		end do

	
		return
	end subroutine



	subroutine wannInterp()
		integer						:: ki, R, a, b, c, n, m
		complex(dp)					:: phase

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
		call writeInterpBands(En_vec)
		




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
		!VELOCITIES
		do ki = 1, nK
			do m = 1, num_wann
				do n = 1, num_wann
					v_mat(1:3,n,m,ki)	=  Ha_mat(1:3,n,m,ki)	- i_dp * dcmplx( En_vec(m,ki) - En_vec(n,ki) ) * A_mat(1:3,n,m,ki) 
				end do
			end do
		end do



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
				tmp(:,:)		= matmul(	A_mat(a,:,:,ki) 	, 	U_mat(:,:,ki)	)
				AaBar(a,:,:)	= matmul(	Ucjg				,	tmp(:,:)		)
				!HaBar
				tmp(:,:)		= matmul(	Ha_mat(a,:,:,ki) 	, 	U_mat(:,:,ki)	)
				HaBar(a,:,:)	= matmul(	Ucjg				,	tmp(:,:)		)
				!DaH
				do n = 1, num_wann
					do m = 1, num_wann
						if(n/=m)	DaH(a,n,m)	= HaBar(a,n,m) /	(	En_vec(m,ki) - En_vec(n,ki)	)
						if(n==m)	DaH(a,n,m)	= dcmplx(0.0_dp)
					end do
				end do
				!OmBar
				do b = 1, 3
					tmp(:,:)		= matmul(	Om_tens(a,b,:,:,ki)	,	U_mat(:,:,ki)	)		
					OmBar(a,b,:,:)	= matmul(	Ucjg				,	tmp(:,:)		)
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


	subroutine writePw90pol( pWann, pConn, pNiuF2, pNiuF3, pPei)
		real(dp),		intent(in)		::	pWann(3), pConn(3), pNiuF2(3), pNiuF3(3), pPei(3)
		real(dp)						:: 	pNiu(3)

			open(unit=600,file='pW90pol.txt',action='write')
		write(600,*)"**************POST W90 POLARIZATION**********************"
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

		return
	end subroutine





end module


