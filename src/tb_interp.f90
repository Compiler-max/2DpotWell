module tb_interp


	use omp_lib
	use util_math,		only:	myExp, i_dp, dp, myLeviCivita, aUtoAngstrm
	use util_sysPara,	only:	Bext, prefactF3, nSolve, nK, atPos
								
	
	use util_w90Interf,	only:	read_band_interp, read_tb_basis
	use pol_Niu,		only:	calcFirstOrdP
	
	implicit none

	private
	public	::					tb_method

	real(dp),		parameter	::	centiMet		= 1e+8_dp		!converts pol from 1/angsroem to 1/cm
	integer						::	num_wann, num_kpts, num_sc


	contains


!public:
	subroutine tb_method()

		real(dp),		allocatable		::	A_interp(:,:,:,:), curv_interp(:,:,:,:), En_interp(:,:), &
											pol0(:,:), centF2(:,:), centF3(:,:)
		complex(dp),	allocatable		::	v_interp(:,:,:,:)
		integer							::	dim, n
		real(dp)						::	polQuantum

		call tb_interpolator(A_interp, curv_interp, En_interp, v_interp)
		!
		allocate(	pol0(	3,	num_wann)	)
		allocate(	centF2(	3,	num_wann)	)
		allocate(	centF3(	3,	num_wann)	)
		!
		!get zero order
		do n = 1, num_wann
			do dim = 1, 3
				pol0(dim,n)	= sum(	A_interp(dim,n,n,:)	) / real(size(A_interp,4),dp)
			end do		
		end do
		pol0 	= pol0 		*	aUtoAngstrm
		write(*,*)	"num_wann | <r> (angs)"
		write(*,*)	"----------------------------"
		do n = 1, num_wann
			write(*,'(i3,a,f16.8,a,f16.8,a,f16.8)')		n," | ",pol0(1,n)," ",pol0(2,n)," ",pol0(3,n)
		end do


		!get first order
		!calcFirstOrdP(polQuantum, centiMet, Bext, prefactF3, Fcurv, Aconn, Velo, En, centers_F2, centers_F3)
		call calcFirstOrdP(polQuantum, centiMet, Bext, prefactF3, curv_interp, A_interp, v_interp, En_interp, centF2, centF3 )


		!convert centers to angstroem
		centF2	= centF2	*	aUtoAngstrm
		centF3	= centF3	*	aUtoAngstrm

		
		!write output file
		call writePolTB(polQuantum, centiMet, pol0, centF2, centF3)

		return
	end subroutine




!private:
	subroutine tb_interpolator(A_interp, curv_interp, En_interp, v_interp)
		real(dp),		intent(out), allocatable		::	A_interp(:,:,:,:), curv_interp(:,:,:,:), En_interp(:,:)
		complex(dp),	intent(out), allocatable		::	v_interp(:,:,:,:)
		real(dp),		allocatable						::	kpts(:,:), en_deriv(:,:,:), Rcell(:,:)
		complex(dp),	allocatable						::	tHopp(:,:,:), rHopp(:,:,:,:)
		integer											::	kpt
		!	
		
		!
		!
		v_interp			= 0.0_dp
		curv_interp 		= 0.0_dp
		write(*,*)	"[tb_interpolator]: WARNING curvature set to zero (irrelevant for 2d)"
		!
		!
		!todo: get interpolated connection
		call read_tb_basis(Rcell, tHopp, rHopp)
		num_wann	= size(tHopp,1)
		num_sc		= size(tHopp,3)
		num_kpts	= num_sc
		write(*,*)	"[tb_interpolator]: read w90 tight binding basis"
		!
		allocate(	kpts(		3,							nK	)		)
		allocate(	A_interp(	3,	num_wann, num_wann, 	nK	)		)
		allocate(	curv_interp(3,	num_wann, num_wann,		nK	)		)
		allocate(	v_interp(	3,	num_wann, num_wann,		nK	)		)

		allocate(	En_interp(			nSolve, 			nK			)		)
		allocate(	en_deriv(	3,		nSolve,				nK			)		)
		!
		!get bands & derivs
		call read_band_interp(kpts, en_interp, en_deriv) !todo: read interp mesh
		write(*,*)	"[tb_interpolator]: read w90 interpolated energies & band derivatives"

		
		!interpolate connection & curvature
		do kpt = 1, num_kpts
			call getConn(kpt, Rcell, kpts,	rHopp, A_interp(	1:3,:,:,kpt)	)
			call getCurv(kpt, Rcell, kpts,	rHopp, curv_interp(	1:3,:,:,kpt)	)
			!
			write(*,'(a,f6.2,a,f6.2,a,f6.2,a)')	"[tb_interpolatior]: done interpolating to k= (",kpts(1,kpt),", ",kpts(2,kpt),", ",kpts(3,kpt),") [a.u.]. "
		end do

		!get velocities
		call calcVeloBLOUNT(A_interp, en_interp, en_deriv, v_interp)
		!
		return 
	end subroutine



subroutine getConn(kpt, Rcell, kpts, rHopp, A_conn)
	!fourier transform of positional operator elements
	!
	!	A_k = sum_{R}	exp[i k.R] <0|ḩat{r}|R>
	!
	integer,		intent(in)		::	kpt
	real(dp),		intent(in)		:: 	Rcell(:,:), kpts(:,:)
	complex(dp),	intent(in)		::	rHopp(:,:,:,:)
	real(dp),		intent(out)		::	A_conn(:,:,:)
	complex(dp)						::	phase
	complex(dp),	allocatable		::	Atmp(:,:)
	integer							::	dim, cell, m, n
	!
	A_conn	= 0.0_dp
	allocate(	Atmp(			num_wann, num_wann			)		)
	!
	!GET VECTOR
	do dim = 1,3
		!SUM OVER CELLS
		Atmp = dcmplx(0.0_dp)	
		do cell = 1, num_sc
			phase				= myExp( 	dot_product(kpts(1:3,kpt),Rcell(1:3,cell))		) 
			Atmp(:,:)			= Atmp(:,:) +  	phase * rHopp(dim,:,:,cell)	
		end do
		!
		!real conversion
		A_conn(dim,:,:)	= dreal(Atmp)
		!DEBUG (real conversion)
		do m = 1, num_wann
			do n = 1, num_wann
				if( abs(dimag(Atmp(n,m)))	> 1e-8_dp ) write(*,*)	"[getConn]: WARNING non vanishing imag of connection at #kpt=",kpt	
			end do
		end do
	end do
	!
	!
	return
end subroutine



subroutine getCurv(kpt, Rcell, kpts, rHopp, F_curv)
	!fourier transform of positional operator elements
	!
	!	OmTens_k^{a,b} = sum_{R}	exp[i k.R]  (	 i R_a <0|ḩat{r_b}|R>   - i R_b <0|ḩat{r_a}|R>  )			EQ(40,vanderbilt2006-interpolation)
	!
	integer,		intent(in)		::	kpt
	real(dp),		intent(in)		:: 	Rcell(:,:), kpts(:,:)
	complex(dp),	intent(in)		:: 	rHopp(:,:,:,:)
	real(dp),		intent(out)		::	F_curv(:,:,:)
	complex(dp)						::	phase
	complex(dp),	allocatable		::	Om_Tens(:,:,:,:), Om_tmp(:,:)
	integer							::	a, b,c,  cell, m, n
	!
	F_curv = 0.0_dp
	allocate(	Om_Tens(3,3, 	num_wann, num_wann)		)
	allocate(	Om_tmp(			num_wann, num_wann)		)
	!
	!GET THE TENSOR
	do b = 1, 3
		do a= 1,3
			!
			!sum over cells
			Om_Tens = dcmplx(0.0_dp)
			do cell = 1, num_sc
				phase				= myExp( 	dot_product(kpts(1:3,kpt),Rcell(1:3,cell))		) 
				Om_Tens(a,b,:,:)	= Om_Tens(a,b,:,:) 	+ i_dp * Rcell(a,cell)* phase * rHopp(b,:,:,cell)
				Om_Tens(a,b,:,:)	= Om_Tens(a,b,:,:) 	- i_dp * Rcell(b,cell)* phase * rHopp(a,:,:,cell)
			end do
			!
			!
		end do
	end do

	!MAP TO VECTOR (eq.(6,vanderbild2006-interpolation))
	do c = 1,3	
		Om_tmp 	= dcmplx(0.0_dp)
		do b = 1, 3
			do a = 1, 3
				Om_tmp(:,:) = Om_tmp(:,:) 	+	myLeviCivita(a,b,c) * Om_Tens(a,b,:,:)
			end do
		end do

		!real conversion
		F_curv(c,:,:)	= Om_tmp(:,:)

		!DEBUG (real conversion)
		do m = 1, num_wann
			do n= 1, num_wann
				if( abs(dimag(Om_tmp(n,m)))	> 1e-8_dp	)		write(*,*)		"[getCurv]: WARNING non zero imag curvature at #kpt=",kpt
			end do
		end do
	end do
	!
	!
	return
end subroutine


subroutine calcVeloBLOUNT(A_conn, En_vec , en_deriv,  v_mat)
		!use Blount 1962 formu.la
		! 
		real(dp),		intent(in)		::	A_conn(:,:,:,:)
		real(dp),		intent(in)		::	En_vec(:,:), en_deriv(:,:,:)
		complex(dp),	intent(out)		::	v_mat(:,:,:,:)
		real(dp),		allocatable		::	v_Band(:,:,:), en_interp(:,:)
		integer							::	n, m, qi
		!
		allocate(	v_Band(	3,	num_wann,	num_kpts)		)
		allocate(	en_interp(	num_wann,	num_kpts)		)
		!
		v_mat	= dcmplx(0.0_dp)
		!
		!DEBUG
		if( size(A_conn,1) /= size(v_mat,1)	)	write(*,*)	"[calcVeloBLOUNT]: A_conn and v_mat have different real space dimensionality"
		if( size(A_conn,2) /= size(v_mat,2) )	stop		"[calcVeloBLOUNT]: A_conn and v_mat have different amount of states covered"
		if( size(A_conn,3) /= size(v_mat,3) )	stop		"[calcVeloBLOUNT]: A_conn and v_mat have different amount of states covered"
		if( size(A_conn,4) /= size(v_mat,4) )	stop		"[calcVeloBLOUNT]: A_conn and v_mat live on different k meshes"
		!
		!FILL MATRIX
		do qi = 1, size(A_conn,4)
			do m = 1, size(v_mat,3)
				do n = 1, size(v_mat,2)
					if(n==m)	v_mat(1:3,n,n,qi)	= dcmplx(	en_deriv(1:3,n,qi) )
					if(n/=m) 	v_mat(1:3,n,m,qi)	= dcmplx(		0.0_dp,		-1.0_dp * (En_vec(m,qi)-En_vec(n,qi)) * A_conn(1:3,n,m,qi)		 )
				end do
			end do
		end do	
		write(*,'(a,i5,a)')	"[calcVeloBLOUNT]: velocities calculated at ",size(A_conn,4)," kpts"
		!
		!
		return
	end subroutine





	subroutine	writePolTB(polQuantum, centiMet, pol0, centF2, centF3)
		real(dp),		intent(in)			::	polQuantum, centiMet, pol0(:,:), centF2(:,:), centF3(:,:)
		integer								::	n, at
		real(dp),		allocatable			::	polC(:,:)

		allocate(	polC(size(pol0,1),size(pol0,2)) )
		polC = 0.0_dp
		!substract atom centers
		do n = 1, size(pol0,2)
			at = mod(n,size(atPos,2))
			if( at== 0) at = size(atPos,2)
			polC(1:2,n)		= pol0(1:2,n) - atPos(1:2,at)*aUtoAngstrm
			!b_H_final(1:2,n)	= b_H_gauge(1:2,n) - atPos(1:2,at)*aUtoAngstrm
			!b_W_final(1:2,n)	= b_W_gauge(1:2,n) - atPos(1:2,at)*aUtoAngstrm
			!
			!ToDo: need niu cent as well ?
		end do



		open(unit=100,file='polInterp.txt',action='write')
		!
		write(100,*)	"# centers are given in angstroem. Electric polarization given in muC/cm"
		!
		!
		!RAW CENTERS
		write(100,*)	"begin centers"
		!----------------------------------------------------------------------------------------------------------------------------------------		
		write(100,*)
		write(100,*)	"nWf	|  <r> - r_atomCenter (ang)"
		write(100,*)	"begin zero_order_centers"
		do n = 1, size(pol0,2)
			write(100,'(i4,a,f16.8,a,f16.8,a,f16.8)')	n," ",pol0(1,n)," ",pol0(2,n)," ",pol0(3,n)
		end do
		write(100,*)	"end zero_order_centers"
		!----------------------------------------------------------------------------------------------------------------------------------------
		!+++
		!----------------------------------------------------------------------------------------------------------------------------------------
		write(100,*)
		write(100,*)	"begin f2_centers"
		do n = 1, size(centF2,2)
			write(100,'(i4,a,f16.8,a,f16.8,a,f16.8)')	n," ",centF2(1,n)," ",centF2(2,n)," ",centF2(3,n)
		end do
		write(100,*)	"end f2_centers"
		!----------------------------------------------------------------------------------------------------------------------------------------
		!+++
		!----------------------------------------------------------------------------------------------------------------------------------------
		write(100,*)
		write(100,*)	"begin f3_centers"
		do n = 1, size(centF3,2)
			write(100,'(i4,a,f16.8,a,f16.8,a,f16.8)')	n," ",centF3(1,n)," ",centF3(2,n)," ",centF3(3,n)
		end do
		write(100,*)	"end f3_centers"
		!----------------------------------------------------------------------------------------------------------------------------------------
		write(100,*)	"end centers"


		!MODIFIED POLARIZATION
		write(100,*)	"begin pol"
		write(100,*)	"		todo"
		write(100,*)	"end pol"

			!TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

		close(100)
		!
		!
		return
	end subroutine



end module tb_interp
