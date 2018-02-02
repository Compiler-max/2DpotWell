module postW90
	use omp_lib
	use mathematics,	only:	dp, PI_dp, i_dp, acc, machineP,  & 
								aUtoAngstrm, aUtoEv, & 
								myExp, myLeviCivita, isHermitian
	use sysPara,		only:	nK, writeBin, seedname, prefactF3, aX, aY, Bext, vol, w90_dir, info_dir ! aX, aY, Bext, vol only needed for output file here
	use wannInterp,		only:	DoWannInterpol
	use polarization,	only:	calcPolWannCent, calcPolViaA
	use semiclassics,	only:	calcFirstOrdP
	use output,			only:	writeVeloEffTB

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
										U_int(:,:,:), Om_mat(:,:,:,:), v_mat(:,:,:,:)


	contains

	

!public
	subroutine effTBmodel()
		real(dp)					:: 	pWann(3), pConn(3), pNiuF2(3), pNiuF3(3), pPei(3)
		integer						::	ki
		logical						::	foundFile
		!
		pWann 	= 0.0_dp
		pConn	= 0.0_dp
		pNiuF2	= 0.0_dp
		pNiuF3	= 0.0_dp 
		pPei	= 0.0_dp
		!read seed_name from eStructure input file
		seed_name= seedName
		!If TB file found do calc
		call readTBsingle( foundFile )
		if( foundFile ) then
			write(*,*)	"[effTBmodel]: done reading eff tb matrices"
			allocate(	A_mat(		3,		num_wann,	num_wann,	nK	)	)		
			allocate(	En_vec(						num_wann	,	nK	)	)
			allocate(	U_int(				num_wann, 	num_wann, 	nK	) 	)
			allocate(	v_mat(		3,		num_wann,	num_wann,	nK	)	)
			allocate(	Om_mat(		3,		num_wann,	num_wann,	nK	)	)
			!
			!
			do ki = 1, nK
				call DoWannInterpol(ki, r_tb, H_tb, R_real, En_vec(:,ki), U_int(:,:,ki), A_mat(:,:,:,ki), Om_mat(:,:,:,ki), v_mat(:,:,:,ki))
			end do
			write(*,*)	"[effTBmodel]: done interpolating to k mesh with nK=",nK
			!
			!calc all desired polarizations
			call polWrapper(pWann, pConn, pNiuF2, pNiuF3, pPei)
		else
			write(*,*)	"[effTBmodel]: did not find input file, no calculations performed"
		end if
		!
		!output file
		call writePw90pol( pWann, pConn, pNiuF2, pNiuF3, pPei)
		call writeConnTxt( A_mat )
		!
		!
		return
	end subroutine










!privat
	subroutine readTBsingle( readSuccess )
		logical,	intent(out)		::	readSuccess
		integer						:: 	stat, cnt, offset, R, n, m, i, mn(2), dumI(3), line15(15)
		real(dp)					::	real2(2), real6(6), real3(3)
		!try opening file
		open(unit=310, iostat=stat, file=w90_dir//seed_name//'_tb.dat', status='old', action='read' )
		if( stat /= 0)  then
			write(*,*) "[readTBsingle]: warning, file seedname_tb.dat not found"
			readSuccess 	= .false.
			recip_latt		= 0.0_dp
			R_real			= 0.0_dp
			H_tb			= dcmplx(0.0_dp)
			r_tb			= dcmplx(0.0_dp)
		else
			readSuccess	= .true.
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
			close(310)
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
				!!DEBUG
				!if( abs(R_real(1,R)-Rcell(1,R)) > machineP ) then
				!	write(*,*) "[readTB]: warning  Rcell and R_real dont match (x comp)" 
				!	write(*,*) "			R_real(x)=",R_real(1,R)," Rcell=",Rcell(1,R)
				!	write(*,*) "			R_real(y)=",R_real(2,R)," Rcell=",Rcell(2,R)
				!end if
				!!if( abs(R_real(2,R)-Rcell(2,R)) > machineP ) write(*,*) "[readTB]: warning Rcell and R_real dont match(y comp)"
			end do
			!
			!
		end if
		!
		!
		return
	end subroutine



	subroutine polWrapper(pWann, pConn, pNiuF2, pNiuF3, pPei)
		real(dp),		intent(out)		:: pWann(3), pConn(3), pNiuF2(3), pNiuF3(3), pPei(3)
		real(dp),		allocatable		:: dummy(:,:), niu_polF2(:,:), niu_polF3(:,:)
		!
		allocate(	dummy(3,num_wann)	)
		!
		pWann	= 0.0_dp
		pConn	= 0.0_dp
		pNiuF2	= 0.0_dp
		pNiuF3	= 0.0_dp
		pPei	= 0.0_dp
		!
		!POLARIZATION CALC
		write(*,*)	"[effTBmodel]: start calculating electric polarization:"
		write(*,'(a,f6.3,a,f6.3,a,f6.3,a)')	"R_real(:,R_null)= (",&
												R_real(1,R_null),", ",R_real(2,R_null),", ",R_real(3,R_null),")."
		!0th order pol
		call calcPolWannCent( wCent, pWann )
		write(*,'(a,f10.4,a,f10.4,a,f10.4,a)')	"[effTBmodel]: pWann=(",pWann(1),", ",pWann(2),", ",pWann(3),")."
		call calcPolViaA(A_mat, dummy)
		pConn(1:3) = sum(dummy(1:3,:))
		write(*,'(a,f10.5,a,f10.5,a,f10.5,a)')	"[effTBmodel]: pConn=(",pConn(1),", ",pConn(2),", ",pConn(3),")."
		!1st order pol
		call calcFirstOrdP(Om_mat,A_mat,v_mat,En_vec,niu_polF2, niu_polF3) !calcFirstOrdP(Fcurv, Aconn, Velo, En, p1F2, p1F3)
		write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[effTBmodel]: pNiuF2=(",pNiuF2(1),", ",pNiuF2(2),", ",pNiuF2(3),")."
		write(*,'(a,e17.10,a,e17.10,a,e17.10,a)')	"[effTBmodel]: pNiuF3=(",pNiuF3(1),", ",pNiuF3(2),", ",pNiuF3(3),")."
		!
		!
		return
	end subroutine



	subroutine writePw90pol( pWann, pConn, pNiuF2, pNiuF3, pPei)
		real(dp),		intent(in)		::	pWann(3), pConn(3), pNiuF2(3), pNiuF3(3), pPei(3)
		real(dp)						:: 	pNiu(3)
		!
		open(unit=600,file=info_dir//'pW90pol.txt',action='write')
		write(600,*)"**************POST W90 POLARIZATION**********************"
		write(600,*)"via effective tight binding model"
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"*"
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



	subroutine writeConnTxt( A_mat )
		complex(dp),	intent(in)		::	A_mat(:,:,:,:)
		integer							::	qi, n, m
		!
		open(unit=350,file=info_dir//'AconnTB.txt',action='write', status='replace')
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





end module postW90

