module util_output
	!module contains several routines for printing and writing data
	use mpi
	use omp_lib

	use util_sysPara 
	use util_math,		only:	dp, PI_dp, machineP, aUtoEv, aUtoAngstrm, aUtoTesla




	implicit none
	private

	public ::	writeMeshInfo, write_K_lattices,					& 
				writePolFile, 										&
				writeEnTXT, readEnTXT, 								& 
				writeVeloHtxt, writeConnTxt, writeVeloEffTB, 		&
				input_info_printer, printTiming,  printBasisInfo
			



	contains















!public:
	subroutine	input_info_printer()
		integer						:: wf, at, at_per_cell, count, cell, wf_per_cell
		!
		!
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"[main]:**************************numerics info*************************"
		write(*,'(a,i3,a,i3,a)')				"[main]: q mesh nQx (input)=",nQx," nQy(derived)=",nQy
		write(*,'(a,f8.3,a,f8.3,a)')			"[main]: dqx=",dqx, " dqy=",dqy," [1/a0]"
		write(*,'(a,i3,a,i6)')					"[main]: interpolation mesh( multiplicator=",k_mesh_multiplier,")        nK=",nK
		write(*,*)								"[main]: basis cutoff parameter  Gcut=",Gcut," [1/a0]"
		write(*,'(a,i7,a,i7,a)')				"[main]: basis function   maximum  nG=",GmaxGLOBAL," of ",nG," trial basis functions"
		write(*,*)								"[main]: only solve for the first nSolve=",nSolve," bands (solution subspace)"

		
	
		
		!
		!
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"[main]:**************************atom info*************************"
		!
		at_per_cell	= nAt / supCx
		write(*,*)								"[main]: total nAt   =", nAt
		write(*,*)								"[main]: nAt per cell=", at_per_cell
		
		write(*,*)								"	#cell |  #atom | x_rel |  y_rel	  |	 radius_x[a0] | radius_y[a0] | 	pot[Hartr]"
		write(*,*)"---------------------------------------------------------------------------"
		count 	= 0
		cell 	= 0
		do at = 1, nAt
			write(*,'(a,i3,a,i3,a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3)')		"	",cell+1," | ",at,"	 | ", relXpos(at)," | ", relYpos(at)," |   ",&
																											 atRx(at), " | ",atRy(at)," | 	",atPot(at)
																													
			!
			count = count + 1
			if( count == at_per_cell) then
				write(*,*)"---------------------------------------------------------------------------"
				cell 	= cell + 1 
				count 	= 0 
			end if
		end do

		!
		!
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"[main]:**************************Trial orbital info*************************"
		wf_per_cell	= nWfs / supCx
		!
		write(*,*)								"[main]: nBands=", nBands
		write(*,*)								"[main]: nWfs  =", nWfs
		write(*,'(a,i3,a)')						"[main]: project ",nWfs/nAt," states onto each atom"
		write(*,*)								"[main]: w90 seed_name= ", seedName
		write(*,*)								"	#cell | #wf | assoc. atom | 	nX| 	nY "
		write(*,*)								"-------------------------------------------------"
		count 	= 0
		cell	= 0
		do wf = 1, nWfs
			write(*,'(a,i3,a,i3,a,i3,a,i3,a,i3)')		"	",cell+1," | ",wf,"	 | ",proj_at(wf),"	 | ",proj_nX(wf)," | ",proj_nY(wf)

			count = count + 1
			if( count == wf_per_cell) then
				write(*,*)"---------------------------------------------------------------------------"
				cell 	= cell + 1 
				count 	= 0 
			end if
		end do
		!
		!
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"[main]:**************************Hamiltonian info*************************"
		!
		write(*,'(a,i3,a,e12.4,a)')				"[main]:	featuring ",nAt," well potentials, deepest well	",minval(atPot(:)),			" Hartree"


		if( doRashba)	write(*,'(a,e12.4,a)')	"[main]:	featuring a Rashba term with prefact 		",aRashba,							" Hartree a0"
		if (doRashba) then
			if( use_px_rashba)	write(*,*)			"[main]:	rashba will enter via aRashba 	* p_x	(x-component of Gvec)"
			if( .not. use_px_rashba)	write(*,*)	"[main]:	rashba will enter via aRashba 	* p_y	(y-component of Gvec)"
		end if
		if( doZeeman)	write(*,'(a,e12.4,a)')	"[main]:	featuring a Zeeman term with prefact 		",0.5_dp*Bext(3),					" Hartree"
		if( doMagHam)	write(*,'(a,e12.4,a)')	"[main]:	featuring a osc. mag. field  prefact 		",0.5_dp*Bext(3)*aX/(2.0_dp*PI_dp),	" Hartree a0"
		write(*,*)								"*"
		if( debugHam)	write(*,*)				"[main]:	debug on: will test Ham for hermiticity"
		!
		!
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"[main]:**************************external field*************************"
		write(*,'(a,f6.3,a,f6.3,a,f6.3,a)')		"[main]: Bext = (",Bext(1)*aUtoTesla,", ",Bext(2)*aUtoTesla, ", ",Bext(3)*aUtoTesla,") (T)"
		!
		!try to print some WARNINGs for to small Gcut
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"[main]:**************************BASIS SET DEBUG*************************"
		call printBasisInfo()
		write(*,*)								"[main]: ...wrote basis set debug info"
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"*"
		write(*,*)								"*"


		return
	end subroutine









	subroutine writeMeshInfo()
		!writes the generated meshes to readable txt file for debuggin purpose
		!
		!
		integer		:: i
		open(unit=100,file=info_dir//'meshInfo.txt',action='write')
		!
		!!R MESH
		!write(100,*)"*******************R POINT MESH******************************"
		!do i = 1, nR
		!	write(100,'(a,i6,a,f15.7,a,f15.7,a,f15.7,a)')				"r(",i,") = (", rpts(1,i) ,", ",rpts(2,i),",",rpts(3,i)," )"
		!end do		
		!
		!G VECTOR
		!write(100,*)"*******************G VECTOR  test MESH******************************"
		!write(100,*)"nG0=",nG0
		!do i = 1, nG
		!	write(100,'(a,i6,a,f12.4,a,f12.4,a)')				"G(",i,") = (", Gtest(1,i) , ", ", Gtest(2,i) , ")"
		!end do
		!!
		!K MESH
		write(100,*)"*******************K POINT MESH******************************"
		do i = 1, size(qpts,2)
			write(100,'(a,i6,a,f15.7,a,f15.7,a)')				"q(",i,") = (", qpts(1,i) , "," , qpts(2,i), " )"
		end do
		!!
		!K 	INTERPOLATION
		write(100,*)"*******************K INTERPOLATION MESH******************************"
		do i = 1, size(kpts,2)
			write(100,'(a,i6,a,f15.7,a,f15.7,a)')				"k(",i,") = (", kpts(1,i) , "," , kpts(2,i), " )"
		end do
		!!
		!ATOM POSITION AND RADII
		write(100,*)"*******************ATOMS******************************"
		do i = 1, size(atPos,2)
			write(100,'(a,i3,a,f8.3,a,f8.3,a,a,f8.3,a,f8.3,a)')	"atom(",i,") at (",atPos(1,i),", ",atPos(2,i),") ",&
																"w. rad = (", atR(1,i), ", ", atR(2,i), " )"	
		end do
		!
		!
		close(100)
		return
	end subroutine



	subroutine write_K_lattices()
		integer					:: stat
		!
		!Q MESH
		open(unit=320, iostat=stat, file=raw_dir//'qpts.dat', status='old')
		if (stat == 0) close(320, status='delete')
		open(unit=320,file='rawData/qpts.dat',form='unformatted',access='stream',action='write')
		write(320) qpts
		close(320)
		!
		!K MESH (INTERPOLATION)
		open(unit=325, iostat=stat, file=raw_dir//'kpts.dat', status='old')
		if (stat == 0) close(325, status='delete')
		open(unit=325,file='rawData/kpts.dat',form='unformatted',access='stream',action='write')
		write(325) kpts
		close(325)
		!
		return 
	end subroutine







	subroutine writeVeloHtxt(velo)
		complex(dp),	intent(in)		:: velo(:,:,:,:) !veloK(3		, 	nWfs	,	nWfs	,	nK		)
		integer		qi, n, m
		complex(dp)						:: vSum
		real(dp)						:: Ederiv

		Ederiv	= aUtoEv * aUtoAngstrm

		open(unit=800,file=info_dir//'veloBerry.txt',action='write')
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


	subroutine writeEnTXT(en)
		real(dp),	intent(in)		:: 	en(:,:)
		integer						::	qi, n
		!
		open(unit=810,file=info_dir//'enABiN.txt',action='write', form='formatted', status='replace')
		write(810,*)		"#	abinitio energies and qpts"
		write(810,'(a,f13.6,a,f13.6,a,f13.6,a)')		"#	B_ext= ",Bext(1)*aUtoTesla, " ",Bext(2)*aUtoTesla," ", Bext(3)*aUtoTesla," T"
		write(810,*)		"#	q_idx	| 			qx(1/ang)		 qy(1/ang)		 qz(1/ang)		|		Energy (eV)"
		do qi = 1, size(en,2)
			do n = 1, size(en,1)
				write(810,'(i5,a,f15.10,a,f15.10,a,f15.10,a,f15.10)')		qi," ",qpts(1,qi)/aUtoAngstrm," ",qpts(2,qi)/aUtoAngstrm," ",0.0_dp," " ,en(n,qi)*aUtoEv
			end do
		end do
		!
		close(810)
		!
		!
		return
	end subroutine


	subroutine readEnTXT(en)
		real(dp),	intent(out)		:: 	en(:,:)
		integer						::	qi, n, f_qi
		real(dp)					::	f_qpt(3)
		!
		open(unit=815,file=info_dir//'enABiN.txt',action='read', form='formatted')
		read(815,*)
		read(815,*)
		read(815,*)
		!
		do qi = 1, size(en,2)
			do n = 1, size(en,1)
				read(815,*) 	 f_qi, f_qpt(1:3), en(n,qi)
			end do
		end do
		!
		close(815)
		!
		!CONVERT FROM EV TO AU
		en = en / aUtoEv
		!
		return
	end subroutine

	subroutine writeVeloEffTB(velo)
		complex(dp),	intent(in)		:: 	velo(:,:,:,:) !veloK(3		, 	nWfs	,	nWfs	,	nK		)
		integer							::	qi, n, m
		complex(dp)						:: vSum

		open(unit=800,file=info_dir//'veloTB.txt',action='write')
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







	subroutine writeConnTxt( A_mat )
		real(dp),		intent(in)		::	A_mat(:,:,:,:)
		integer							::	qi, n, m
		!
		open(unit=350,file=info_dir//'AconnBerry.txt',action='write', status='replace')
		write(350,*)	"connection calculated via berryMethod"
		write(350,*)	"n m real(A_x) imag(A_x) real(A_y) imag(A_y) real(A_z) imag(A_z)"
		do qi = 1, size(A_mat,4)
			write(350,*)	"qi=",	qi
			do m = 1, size(A_mat,3)
				do n = 1, size( A_mat,2)
					write(350,'(i3,a,i3,a,i5,a,e15.6,a,e15.6,a,e15.6,a)')	n," ",m," ",qi," (",A_mat(1,n,m,qi),", ",A_mat(2,n,m,qi),", ",A_mat(3,n,m,qi),")"
				end do
			end do
		end do
		close(350)
		!
		return
	end subroutine



	subroutine writePolFile(polQuantum, centiMet, w_centers, b_H_gauge, b_W_gauge, niu_centF2, niu_centF3, essin_centF3)	!writePolFile(pWann, pBerry, pNiu, pPei )
		!
		real(dp),		intent(in)		::	polQuantum, centiMet, w_centers(:,:),  b_H_gauge(:,:), b_W_gauge(:,:), niu_centF2(:,:), niu_centF3(:,:), essin_centF3(:,:) !CENTERS IN ANGSTROEM
		real(dp)						:: 	pWann(3), pBerryH(3),pBerryW(3), &
											pNiuF2(3), pNiuF3(3), pNiu(3), pFirst(3), Btesla(3)
		real(dp),		allocatable		::	w_final(:,:), b_H_final(:,:), b_W_final(:,:)
		integer							::	n, x
		


		!GET CORRECTED VALUES
		allocate(		w_final( 	size(w_centers,1),	size(w_centers,2) )		)
		allocate(		b_H_final(	size(b_H_gauge,1),	size(b_H_gauge,2) )		)
		allocate(		b_W_final(	size(b_W_gauge,1),	size(b_W_gauge,2) )		)		

		!zero init
		w_final 	= 0.0_dp
		b_H_final	= 0.0_dp
		b_W_final	= 0.0_dp


		!substract atom centers
		do n = 1, size(w_centers,2)
			w_final(1:2,n)		= w_centers(1:2,n) - atPos(1:2,proj_at(n))*aUtoAngstrm
			b_H_final(1:2,n)	= b_H_gauge(1:2,n) - atPos(1:2,proj_at(n))*aUtoAngstrm
			b_W_final(1:2,n)	= b_W_gauge(1:2,n) - atPos(1:2,proj_at(n))*aUtoAngstrm
			!
			!ToDo: need niu cent as well ?
		end do


		!SUM OVER STATES
		do x = 1, 3
			pWann(	x)	= sum(	w_final(	x, :)		)	 * polQuantum * centiMet ! muC/cm 
			pBerryH(x)	= sum(	b_H_final(	x, :)		)	 * polQuantum * centiMet ! muC/cm
			pBerryW(x)	= sum(	b_W_final(	x, :)		)	 * polQuantum * centiMet ! muC/cm
			pNiuF2(	x)	= sum(	niu_centF2(	x, :)		)	 * polQuantum * centiMet ! muC/cm
			pNiuF3(	x)	= sum(	niu_centF3(	x, :)		)	 * polQuantum * centiMet ! muC/cm
		end do
		!
		pNiu 	=	( pNiuF2 	+ pNiuF3 	)
		pFirst 	=	( pBerryW 	+ pNiu 		) 


	
		!	
		open(unit=600,file=info_dir//'polOutput.txt',action='write')
		write(600,*)"**************POLARIZATION OUTPUT FILE**********************"
		write(600,*)" via wavefunction method"
		if(doMagHam )	write(600,*)"MAGNETIC HAMILTONIAN USED!!!"
		write(600,*)"maximum amount of basis function basis =",GmaxGlobal
		write(600,*)"nQx=",nQx,"; nQy=",nQy
		write(600,*)"pol quantum = e/vol =",polQuantum
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"begin gCut"
		write(600,*)	Gcut
		write(600,*)"end gCut"

		write(600,*)"begin mp_grid"
		write(600,'(i3,a,i3,a,i3)')	nQx," ",nQy," ",1
		write(600,*)"end mp_grid"
		!
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"**************UNIT CELL (Angstroem): "
		write(600,*)"begin unit_cell"
		write(600,*)" ",aX*aUtoAngstrm," ",aY*aUtoAngstrm
		write(600,*)"end unit_cell"
		!
		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"**************ATOMS:"
		write(600,*)		"	at | centers [Å] | V [eV]"
		do n = 1, size(atPos,2)
				write(600,'(i3,a,f14.6,a,f14.6,a ,f14.6)')	n," | ",atPos(1,n)*aUtoAngstrm,", ",atPos(2,n)*aUtoAngstrm," 	| ",atPot(n)*aUtoEv
		end do
		write(600,*)"begin atPot"
		write(600,"(100f15.5)") ( atPot(n)*aUtoEv, n=1,size(atPot) )
		write(600,*)"end atPot"

		write(600,*)"*"
		write(600,*)"*"
		write(600,*)"**************PERTURBATION:"
		!
		!
		Btesla(1:3)	= Bext(1:3) * aUtoTesla
		write(600,*)"begin magnetic_field"
		write(600,'(a,e13.6,a,e13.6,a,e13.6)')	"  ", 	Btesla(1),	" ",	Btesla(2)," ", Btesla(3)
		write(600,*)"end magnetic_field"
		write(600,*)"*"
		write(600,*)"begin alpha_rashba"
		write(600,*)	aRashba*aUtoEv*aUtoAngstrm
		write(600,*)"end alpha_rashba"
		write(600,*)"*"
		write(600,*) "states considered for perturbation:" 
		write(600,*) "begin nSolve"
		write(600,*)	nSolve
		write(600,*) "end nSolve" 
		if( norm2(Bext) > machineP ) then
			write(600,*)	"no magnetic field applied ( norm of field is zero)"
		else
			write(600,'(a,f6.2,a,f6.2,a,f6.2,a)')	"Bext= (", 	Bext(1),	", ",	Bext(2),", ", Bext(3),	") [a.u.]"			
			write(600,'(a,f6.2,a,f6.2,a,f6.2,a)')	"Bext= (", 	Btesla(1),	", ",	Btesla(2),", ", Btesla(3),	") [T]"
		end if
		write(600,*)"*"

		!
		write(600,*)"**************POL:"
		write(600,*)		" method | 	P[muC/cm]	"
		write(600,'(a,e13.4,a,e13.4,a,e13.4,a)')		" w90    | 	(",	pWann(1)	,", ",	pWann(2)	,", ",	pWann(3)	,")"
		write(600,'(a,e13.4,a,e13.4,a,e13.4,a)')		" B_H    | 	(",	pBerryH(1)	,", ",	pBerryH(2)	,", ",	pBerryH(3)	,")"
		write(600,'(a,e13.4,a,e13.4,a,e13.4,a)')		" B_W    | 	(",	pBerryW(1)	,", ",	pBerryW(2)	,", ",	pBerryW(3)	,")"
		write(600,'(a,e13.4,a,e13.4,a,e13.4,a)')		" nF2    | 	(",	pNiuF2(1)	,", ",	pNiuF2(2)	,", ",	pNiuF2(3)	,")"
		write(600,'(a,e13.4,a,e13.4,a,e13.4,a)')		" nF3    | 	(",	pNiuF3(1)	,", ",	pNiuF3(2)	,", ",	pNiuF3(3)	,")"
		write(600,'(a,e13.4,a,e13.4,a,e13.4,a)')		" niu    | 	(",	pNiu(1)		,", ",	pNiu(2)		,", ",	pNiu(3)		,")"
		write(600,'(a,e13.4,a,e13.4,a,e13.4,a)')		" p0+p1  | 	(",	pFirst(1)	,", ",	pFirst(2)	,", ",	pFirst(3)	,")"



		!STATE RESOLVED
		write(600,*)	"band resolved results:"
		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"********wannier90 centers***********************"
		write(600,*)		" #wf | 	<r>-atPos [Å]	| p[\{mu}C/cm]"
		do n = 1, size(w_centers,2)
				write(600,'(i3,a,f14.6,a,f14.6,a,f14.6,a,a,e13.4,a,e13.4,a)') n,"  | ",w_final(1,n),", ", w_final(2,n),",",w_final(3,n), "  | ",&
																"(",w_final(1,n)*polQuantum*centiMet,", ", w_final(2,n)*polQuantum*centiMet, ")."
		end do
		write(600,'(a,f14.6,a,f14.6,a,f14.6,a)')	"	sum	",sum(w_final(1,:))," ",&
																sum(w_final(2,:))," ",&
																sum(w_final(3,:))," | "
		write(600,*)	"begin zero_order"
		write(600,'(a,e16.9,a,e16.9,a,e16.9,a)')	"		 ",sum(w_final(1,:))*polQuantum*centiMet," ",&
															sum(w_final(2,:))*polQuantum*centiMet," ",&
															sum(w_final(3,:))*polQuantum*centiMet," #muC/cm"
		write(600,*)	"end zero_order"
		
		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"********Berry centers (k-space integral) (H-gauge)***********************"
		write(600,*)		" #state | 	<r>[Å]			| 	p[	\{mu}C/cm	]"
		do n = 1, size(b_H_gauge,2)
			write(600,'(i3,a,f6.2,a,f6.2,a,f6.2,a,a,e13.4,a,e13.4,a)') n,"  | ",b_H_final(1,n),", ", b_H_final(2,n),",",b_H_final(3,n), "  | ",&
																"(",b_H_final(1,n)*polQuantum*centiMet,", ", b_H_final(2,n)*polQuantum*centiMet, ")."
		end do
		write(600,'(a,e16.9,a,e16.9,a)')	"sum | 				(",sum(b_H_final(1,:))*polQuantum*centiMet,", ",sum(b_H_final(2,:))*polQuantum*centiMet,")."
		!		

		
		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"********Berry centers (k-space integral) (W-gauge)***********************"
		write(600,*)		" #state | 	<r>[Å]			| 	p[	\{mu}C/cm	]"
		do n = 1, size(b_W_gauge,2)
			write(600,'(i3,a,f6.2,a,f6.2,a,f6.2,a,a,e13.4,a,e13.4,a)') n,"  | ",b_W_final(1,n),", ", b_W_final(2,n),",",b_W_final(3,n), "  | ",&
																"(",b_W_final(1,n)*polQuantum*centiMet,", ", b_W_final(2,n)*polQuantum*centiMet, ")."
		end do
		
		write(600,'(a,e16.9,a,e16.9,a)')	"sum | 				",sum(b_W_final(1,:))*polQuantum*centiMet," ",sum(b_W_final(2,:))*polQuantum*centiMet," #muC/cm"






		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"********Niu F2 (first order contribution)***********************"
		write(600,*)		" #state | 	<r>[Å]			| 	p[	\{mu}C/cm	]"
		do n = 1, size(niu_centF2,2)
			write(600,'(i3,a,f14.6,a,f14.6,a,f14.6,a,a,e13.4,a,e13.4,a)') n,"  | ",niu_centF2(1,n),", ", niu_centF2(2,n),",",niu_centF2(3,n), "  | ",&
																"(",niu_centF2(1,n)*polQuantum*centiMet,", ", niu_centF2(2,n)*polQuantum*centiMet, ")."
		end do
		write(600,'(a,f14.6,a,f14.6,a,f14.6,a)')	"	sum	",sum(niu_centF2(1,:))," ",&
																sum(niu_centF2(2,:))," ",&
																sum(niu_centF2(3,:))," | "
		write(600,*)	"begin niu_f2"
		write(600,'(a,e16.9,a,e16.9,a,e16.9,a)')	"				",sum(niu_centF2(1,:))*polQuantum*centiMet," ",&
																sum(niu_centF2(2,:))*polQuantum*centiMet," ",&
																sum(niu_centF2(3,:))*polQuantum*centiMet," #muC/cm"
		write(600,*)	"end niu_f2"

		!
		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"********Niu F3 (first order contribution)***********************"
		write(600,*)		" #state | 	<r>[Å]			| 	p[	\{mu}C/cm	]"
		do n = 1, size(niu_centF3,2)
			write(600,'(i3,a,f14.6,a,f14.6,a,f14.6,a,a,e13.4,a,e13.4,a)') n,"  | ",niu_centF3(1,n),", ", niu_centF3(2,n),",",niu_centF3(3,n), "  | ",&
																"(",niu_centF3(1,n)*polQuantum*centiMet,", ", niu_centF3(2,n)*polQuantum*centiMet, ")."
		end do
		write(600,'(a,f14.6,a,f14.6,a,f14.6,a)')	"	sum	",sum(niu_centF3(1,:))," ",&
																sum(niu_centF3(2,:))," ",&
																sum(niu_centF3(3,:))," | "
		write(600,*)	"begin niu_f3"
		write(600,'(a,e16.9,a,e16.9,a,e16.9,a)')	"				",sum(niu_centF3(1,:))*polQuantum*centiMet," ",&
																sum(niu_centF3(2,:))*polQuantum*centiMet," ",&
																sum(niu_centF3(3,:))*polQuantum*centiMet," #muC/cm"
		write(600,*)	"end niu_f3"

			!
		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"*"
		write(600,*)	"********essin F3 (first order contribution)***********************"
		write(600,*)		" #state | 	<r>[Å]			| 	p[	\{mu}C/cm	]"
		do n = 1, size(essin_centF3,2)
			write(600,'(i3,a,f14.6,a,f14.6,a,f14.6,a,a,e13.4,a,e13.4,a)') n,"  | ",essin_centF3(1,n),", ", essin_centF3(2,n),",",essin_centF3(3,n), "  | ",&
																"(",essin_centF3(1,n)*polQuantum*centiMet,", ", essin_centF3(2,n)*polQuantum*centiMet, ")."
		end do
		write(600,'(a,f14.6,a,f14.6,a,f14.6,a)')	"	sum	",sum(essin_centF3(1,:))," ",&
																sum(essin_centF3(2,:))," ",&
																sum(essin_centF3(3,:))," | "
		write(600,*)	"begin essin_f3"
		write(600,'(a,e16.9,a,e16.9,a,e16.9,a)')	"				",sum(essin_centF3(1,:))*polQuantum*centiMet," ",&
																sum(essin_centF3(2,:))*polQuantum*centiMet," ",&
																sum(essin_centF3(3,:))*polQuantum*centiMet," #muC/cm"
		write(600,*)	"end essin_f3"

		close(600)
		!
		!
		return
	end subroutine



	


	subroutine printBasisInfo()
		!nG to Gcut RATIO
		if(			nG		< 		vol * Gcut * dsqrt(Gcut) /	(2.0_dp*PI_dp**2)		) then
			stop	"[BasisInfo]:	volume: increase Gcut"
		else
			write(*,*)	"[BasisInfo]: 	volume: accepted Gcut"
		end if
		!
		return 
	end subroutine



	subroutine printTiming(alloT,hamT,berryT,pwT,mastT) 	!call printTiming(alloT,hamT,berryT, postWT,mastT)
		real,	intent(in)	:: alloT, hamT, pwT, berryT, mastT
		!
		print '    ("r&alloc  time spend     = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									alloT 				, 100.0_dp*alloT		   	/mastT
		print '    ("ham-solver time spend   = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									hamT 				, 100.0_dp*hamT		   	/mastT
		print '    ("berry method            = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									berryT				, 100.0_dp*berryT		   	/mastT			
		print '    ("post w90 - eff TB       = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									pwT				, 100.0_dp*pwT		   	/mastT
		print '    ("other    time spend     = ",f16.4," seconds = ",f16.4,"% of overall time")',& 
									(mastT-alloT-hamT-berryT) 	, 100.0_dp*(mastT-alloT-hamT-berryT)/mastT
		
		print '    ("overall  time spend     = ",f16.4," seconds.")', mastT
		!
		return
	end subroutine




end module util_output


