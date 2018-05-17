module util_sysPara
	!this modules reads the input file and generates the meshes accordingly
	use mpi
	use util_math, 		only:	dp, PI_dp, setAcc, acc, machineP, aUtoTesla,  aUtoAngstrm, aUtoEv
	use m_config
	implicit none
	!#include "mpif.h"

	private
				!public routines:
	public :: 	readInp, 																& 
				!public para:
				!unit cell:
				dim, aX, aY, vol, recpLatt, supCx,										& 
				!atoms:
				nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot,atPos, atR,  			&
				!cutoff:
				nG, nGq, Gmax, GmaxGLOBAL, GminGLOBAL, Gcut, nSolve, glob_min_gap,		&
				!grids:
				qpts, kpts, 															&
				nQ, nQx, nQy, 															&
				k_mesh_multiplier, nKx, nKy, nK,  										&
				nSC, nSCx, nSCy, 														&
				dqx, dqy, dkx, dky,														&
				nR, nRx, nRy, nRz,   													&	
				!w90:
				nShells, nw90it, shells, 												&
				nBands, nWfs, proj_at, proj_nX, proj_nY,  								&
				seedName, w90_dir, info_dir, mkdir, raw_dir,							&
				!perturbations:	
				aRashba, use_px_rashba, Bext, prefactF3, 								&
				!switches:
				debugHam, doSolveHam, doMagHam, doRashba, doZeeman, doVdesc, 			&
				useBloch, doPw90, pw90GaugeB, 	do_w90plot,	 							&
				doBerry, doNiu, doGaugBack, fastConnConv,	 							&
				!mpi:
				myID, nProcs, root, ierr, qChunk


	!
	integer  										:: 	dim=2, nAt=0, nAt_inp=0, supCx=1,												& 
														nG, nGx,nGy, Gmax, GmaxGLOBAL, GminGLOBAL,nSolve=20,							&  
														nQx=1, nQy=1,nQ , 																&
														nSCx=1, nSCy=1, nSC, 															& 
														nKx=1, nKy=1, nK,k_mesh_multiplier=1,   										&
														nShells, nw90it, 																&
														nRx=10, nRy=10, nRz=2, nR, nBands=1,nBands_inp=1,nWfs=1,nWfs_inp=1, 			&
														myID, nProcs, ierr, qChunk
	integer,	parameter							::	root=0
	real(dp) 										::	aX=0.0_dp,aX_inp=0.0_dp, aY=0.0_dp, vol=0.0_dp, Gcut=2*PI_dp, thres,			& 
														glob_min_gap,																	&
														dqx, dqy, dkx, dky,																& 
														aRashba=0.0_dp, B0, Bext(3)	, prefactF3, recpLatt(2,2)
	character(len=3)								::	seedName										
	character(len=9)								::	w90_dir	="w90files/"
	character(len=7)								::	info_dir="output/"
	character(len=8)								::	raw_dir	="rawData/", mkdir="mkdir ./"	!to use with system(mkdir//$dir_path) 
	integer,	allocatable,	dimension(:)		::	nGq, shells, &
														proj_at, proj_nX, proj_nY, 														&
														proj_at_inp, proj_nX_inp, proj_nY_inp		
	!
	real(dp),	allocatable,	dimension(:)		::	relXpos_inp, relYpos_inp, atRx_inp, atRy_inp, atPot_inp, dVpot_inp, 			&
														relXpos, relYpos, atRx, atRy, atPot, dVpot
	!
	real(dp),	allocatable,	dimension(:,:)		::	atPos_inp, atR_inp, 															&
														atPos, atR, 																	&
														qpts, kpts 

	logical											::	debugHam,  																		&
														doSolveHam, doVdesc,  doRashba, doZeeman, doMagHam, 							&
														use_px_rashba, 																	&
														useBloch, doPw90, pw90GaugeB, 	do_w90plot,										& 
														doBerry, doNiu, doGaugBack, fastConnConv


	contains
!public:
	subroutine readInp()
		!root reades input parameters from input.txt file 
		!calls mesh generation subroutines
		!and bcasts everything around
		!the first array index is always the x,y value of the vector
		logical				:: dir_exists
		!
		!read inputs and allocates arrays
		if( myID == root ) then
			write(*,*)	"hello from root, will read the input file now"
			call rootRead()
		end if
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		call bcastPARAM() 
		!
		!qpts per mpi process
		qChunk = nQ/nProcs
		!
		!fill arrays
		call setAcc(thres)
		call qmeshGen()
		!call rmeshGen()
		call popGvec()
		call popAtPos()
		call popAtR()
		call kWmeshGen()
		!
		!MAKE DIRECTORIES
		if( myID == root ) then
			inquire(directory=w90_dir, exist=dir_exists)
			if( .not. dir_exists )	call system(mkdir//w90_dir)
			!
			inquire(directory=info_dir, exist=dir_exists)
			if( .not. dir_exists )	call system(mkdir//info_dir)
			!
			inquire(directory=raw_dir, exist=dir_exists) 
			if( .not. dir_exists )	call system(mkdir//raw_dir)
		end if

		!
		return
	end subroutine



!privat	
!READ & DISTRIBUTION ROUTINES
	subroutine rootRead()
		type(CFG_t) :: my_cfg
		
		!OPEN FILE
		call CFG_read_file(my_cfg,"input.txt")
		
		!READ SCALARS
		![unitCell]
		call CFG_add_get(my_cfg,	"unitCell%aX"      	,	aX_inp  	  ,	"length of unit cell in agnstroem"		)
		call CFG_add_get(my_cfg,	"unitCell%aY"      	,	aY  	   	,	"length of unit cell in agnstroem"		)
		call CFG_add_get(my_cfg,	"unitCell%supCx"    ,	supCx  	   	,	"length of unit cell in agnstroem"		)
		![atoms]
		call CFG_add_get(my_cfg,	"atoms%nAt"			,	nAt_inp		,	"number of atoms per unit cell"			)
		![wann]
		call CFG_add_get(my_cfg,	"wann%nBands"		,	nBands_inp	 ,	"# of bands to project onto trial orbs"	)
		call CFG_add_get(my_cfg,	"wann%nWfs"			,	nWfs_inp	 ,	"# wannier functions to generate"		)
		![field]	
		call CFG_add_get(my_cfg,	"field%B0"			,	B0			,	"scaling fact. of ext. magnetic field"	)
		call CFG_add_get(my_cfg,	"field%Bext"		,	Bext		,	"vector of ext. magnetic field"			)
		![rashba]
		call CFG_add_get(my_cfg,	"rashba%use_px_rashba",	use_px_rashba,	"add px term to Ham. else add py "		)
		call CFG_add_get(my_cfg,	"rashba%aRashba"	,	aRashba		,	"rashba coeff in eV Ang"				)
		![numerics]
		call CFG_add_get(my_cfg,	"numerics%Gcut"		,	Gcut	    ,	"k space cut of parameter"				)
		call CFG_add_get(my_cfg,	"numerics%nSolve"	,	nSolve	    ,	"number of eigenstates to find"			)
		call CFG_add_get(my_cfg,	"numerics%nQx"     	,	nQx      	,	"amount of k points used"				)
		!call CFG_add_get(my_cfg,	"numerics%nQy"     	,	nQy      	,	"amount of k points used"				)
		call CFG_add_get(my_cfg,	"numerics%thres"    ,	thres      	,	"threshold for overlap WARNINGs"		)
		![ham]
		call CFG_add_get(my_cfg,	"ham%debugHam"		, 	debugHam	,	"switch for debuging tests in solveHam"	)	
		call CFG_add_get(my_cfg,	"ham%doVdesc"		,	doVdesc		,	"switch on lin. desc. pot inside wells"	)
		call CFG_add_get(my_cfg,	"ham%doZeeman"		,	doZeeman	,	"include B-field via Zeeman in wells"	)
		call CFG_add_get(my_cfg,	"ham%doMagHam"		,	doMagHam	,	"include osc. B-field via peierls sub"	)
		call CFG_add_get(my_cfg,	"ham%doRashba"		,	doRashba	,	"switch for rashba term"				)
		![methods]
		call CFG_add_get(my_cfg,	"methods%doSolveHam",	doSolveHam	,	"solve electronic structure or read in"	)
		call CFG_add_get(my_cfg,	"methods%doPw90"	,	doPw90		,	"read in the matrices in wann base	"	)	
		call CFG_add_get(my_cfg,	"methods%doBerry"	,	doBerry		,	"switch on/off 	berry( unk) method "	)
		![w90]
		call CFG_add_get(my_cfg,	"w90%seedName"		, 	seedName	,	"seedName for wannier files(char len=3)")
		call CFG_add_get(my_cfg,	"w90%useBloch"		,	useBloch	,	"use bloch phase for projections	"	)
		call CFG_add_get(my_cfg,	"w90%nShells"		,	nShells		,	"number of shells to use for FD k-space")
		call CFG_add_get(my_cfg,	"w90%nw90it"		, 	 nw90it		,	"number of iterations for wannnierisat,")
		call CFG_add_get(my_cfg,	"w90%pw90GaugeB"	,	pw90GaugeB	,	"logical for switching gauge trafo	   ")
		![w90plot]
		!call CFG_add_get(my_cfg,	"w90plot%nSCx"	     ,	nSCx   		,	"#			supercells "				)
		!call CFG_add_get(my_cfg,	"w90plot%nSCy"	     ,	nSCy   		,	"#			supercells "				)
		call CFG_add_get(my_cfg,	"w90plot%do_w90plot",	do_w90plot	,	"switch for plotting wannier & unk"		)
		call CFG_add_get(my_cfg,	"w90plot%nRx"    	,	nRx      	,	"spacing of real space plot grid"		)
		call CFG_add_get(my_cfg,	"w90plot%nRy"    	,	nRy      	,	"spacing of real space plot grid"		)
		call CFG_add_get(my_cfg,	"w90plot%nRz"		,	nRz			,	"spacing of real space plot grid"		)
		![w90interp]
		call CFG_add_get(my_cfg,	"w90interp%k_mesh_multiplier"		,	k_mesh_multiplier		,	"# k x points of interpolation mesh"	)		
		!call CFG_add_get(my_cfg,	"w90interp%nKx"		,	nKx			,	"# k x points of interpolation mesh"	)
		!call CFG_add_get(my_cfg,	"w90interp%nKy"		,	nKy			,	"# k x points of interpolation mesh"	)
		![berry]
		call CFG_add_get(my_cfg,	"berry%fastConnConv",fastConnConv	,	"try faster converging fd formula"		)
		call CFG_add_get(my_cfg,	"berry%doNiu"		,	doNiu		,	"switch for nius first order pol"		)
		call CFG_add_get(my_cfg,	"berry%doGaugBack"	,	doGaugBack	,	"switch for trafo: Wann to Ham gauge"	)
		![semiclassics]
		call CFG_add_get(my_cfg,	"semiclassics%prefactF3"	,	prefactF3,	"real prefactor for F3 "			)
		



		!reformulate super cell
		aX				= real(supCx,dp) * aX_inp
		nAt				= supCx * nAt_inp
		nBands			= supCx * nBands_inp
		nWfs			= supCx * nWfs_inp


		!
		if( mod(aX/aY,1.0_dp) > 1e-7_dp)	stop	"wARNING unit length (aX/aY) should be integer multiples from each other "

		write(*,*)	"[rootRead]: scale nQy with int(aX/aY)=",int(aX/aY)
		nQy				= nQx * int(aX/aY) 
		nSCx			= nQx + 1
		nSCy			= nQy + 1
		nKx				= nQx * k_mesh_multiplier + 1
		nKy				= nQy * k_mesh_multiplier + 1
		write(*,*)		"[rootRead]:	WARNING nQy, nSCx, nSCy, nKx, nKy overwritten (handles uniform q/k mesh automatically)"


		!get derived quantities
		call getTestGridSize(nGx, nGy)	
		nG				= 	nGx*nGy
		vol				=	aX 		* 	aY
		nR 				= 	nRx 	*	nRy   * nRz
		nQ 				= 	nQx 	*	nQy
		nSC 			=	nSCx	*	nSCy
		nK 				=	nKx		*	nKy
		!
		!

		!scale the field
		Bext=	B0 		* 	Bext	/ aUtoTesla
		!
		!convert rashba coeff to a.u.
		aRashba	= aRashba /	( aUtoEv * aUtoAngstrm	)		!from (eV Ang) to (Eh a0)
		!
		!setup reciprocal lattice
		recpLatt		= 0.0_dp
		recpLatt(1,1)	= 2.0_dp * PI_dp * aY / vol
		recpLatt(2,2)	= 2.0_dp * PI_dp * aX / vol
		!
		!
		!atoms (input)
		allocate(	relXpos_inp(nAt_inp)		)
		allocate(	relYpos_inp(nAt_inp)		)
		allocate(	atRx_inp(nAt_inp)			)
		allocate(	atRy_inp(nAt_inp)			)
		allocate(	atPot_inp(nAt_inp)			)
		allocate(	dVpot_inp(nAt_inp)			)
		!
		!wann (input)
		allocate(	proj_at_inp(nWfs_inp)		)
		allocate(	proj_nX_inp(nWfs_inp)		)
		allocate(	proj_nY_inp(nWfs_inp)		)
		!
		call allocateArrays()
		!
		!read the arrays:
		![atoms]
		call CFG_add_get(my_cfg,	"atoms%relXpos"		,	relXpos_inp		,	"relative positions in unit cell"		)
		call CFG_add_get(my_cfg,	"atoms%relYpos"		,	relYpos_inp		,	"relative positions in unit cell"		)
		call CFG_add_get(my_cfg,	"atoms%atRx"		,	atRx_inp		,	"radius of each atom in angstroem"		)
		call CFG_add_get(my_cfg,	"atoms%atRy"		,	atRy_inp		,	"radius of each atom in angstroem"		)
		call CFG_add_get(my_cfg,	"atoms%atPot"		,	atPot_inp		,	"potential depth in hartree"			)
		call CFG_add_get(my_cfg,	"atoms%dVpot"		,	dVpot_inp		,	"potential gradient"					)
		![wann]
		call CFG_add_get(my_cfg,	"wann%projAt"		,	proj_at_inp		,	"list of atoms to project to (init U)"	)
		call CFG_add_get(my_cfg,	"wann%projnX"		,	proj_nX_inp		,	"on which state to project x-dir"		)
		call CFG_add_get(my_cfg,	"wann%projnY"		,	proj_nY_inp		,	"on which state to project y_dir"		)
		![w90]
		call CFG_add_get(my_cfg,	"w90%shells"		, 	shells		,	"list of shells used for FD connection"	)
	
		!project to super cell
		call superCell_projector()
		!
		!
		return
	end subroutine



	subroutine superCell_projector()
		!projects onto a superCell in x direction
		integer							::	cell, at, wf, index


		!!atoms (input)
		!allocate(	relXpos_inp(nAt_inp)		)
		!allocate(	relYpos_inp(nAt_inp)		)
		!allocate(	atRx_inp(nAt_inp)			)
		!allocate(	atRy_inp(nAt_inp)			)
		!allocate(	atPot_inp(nAt_inp)			)
		!allocate(	dVpot_inp(nAt_inp)			)
		
		!!wann (super Cell)
		!allocate(	proj_at_inp(nWfs_inp)		)
		!allocate(	proj_nX_inp(nWfs_inp)		)
		!allocate(	proj_nY_inp(nWfs_inp)		)
		!!allocate(	atPos_inp(dim,nAt_inp)		

		do cell = 0, supCx-1
			!project atom arrays
			do at = 1, nAt_inp
				index				= cell*nAt_inp + at
				!
				relXpos(index)		= ( real(cell,dp) + relXpos_inp( at ) )	*	aX_inp 	/ aX 	! = 	(cell offset + local position) / superCellSize 
				!
				relYpos(index)		= relYpos_inp(at)
				atRx(index)			= atRx_inp(at)
				atRy(index)			= atRy_inp(at)
				atPot(index)		= atPot_inp(at)
				dVpot(index)		= dVpot_inp(at)
			end do
			!
			!
			!project Wannier functions
			do wf = 1, nWfs_inp
				proj_at(	cell*nWfs_inp + wf)	=	cell*nAt_inp + proj_at_inp(wf)
				!
				proj_nX(	cell*nWfs_inp + wf)	=	proj_nX_inp(wf)
				proj_nY(	cell*nWfs_inp + wf)	=	proj_nY_inp(wf)
			end do
		end do
		write(*,*)	"[superCell_projector]: finished projection"
		!
		return
	end subroutine






	subroutine bcastPARAM()
		integer							:: ierr
		!
		!
		![unitCell]
		call MPI_Bcast( aX_inp		,		1	,	MPI_DOUBLE_PRECISION	, 	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( aX			,		1	,	MPI_DOUBLE_PRECISION	, 	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( aY			,		1	,	MPI_DOUBLE_PRECISION	, 	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( supCx		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		![atoms]
		call MPI_Bcast( doVdesc		, 		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)	
		call MPI_Bcast( nAt_inp		, 		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)	
		call MPI_Bcast( nAt			, 		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)		
		![wann]
		call MPI_Bcast( nBands_inp	,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)		
		call MPI_Bcast( nBands		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)		
		call MPI_Bcast( nWfs_inp	,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( nWfs		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)		
		![field]
		call MPI_Bcast(	B0 			,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)		
		call MPI_Bcast(	Bext 		,		3	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		![rashba]
		call MPI_Bcast( use_px_rashba,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast( aRashba		,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		![numerics]
		call MPI_Bcast(	nGx			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nGy			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)		
		call MPI_Bcast(	nG			,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	Gcut		,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nSolve		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nQx			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nQy			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nQ			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nSCx		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nSCy		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nSC			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nKx			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nKy			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( nK			,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast( nR			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast(	thres		,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		![ham]
		call MPI_Bcast( debugHam	,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)		
		call MPI_Bcast(	doVdesc		,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	doZeeman	,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	doMagHam	,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( doRashba	,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)
		![methods]
		call MPI_Bcast( doSolveHam	,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( doPw90		,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( doBerry		,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)
		![w90]
		call MPI_Bcast( seedName	,		3	, 		MPI_CHARACTER		,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( useBloch	,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( nShells		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( nw90it		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)		
		call MPI_Bcast( pw90GaugeB	,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)
		![w90plot]
		call MPI_Bcast( do_w90plot ,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast(	nRx			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nRy			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( nRz			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD,	ierr)
		![berry]
		call MPI_Bcast(fastConnConv	,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast( doNiu		,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast( doGaugBack	,		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD,	ierr)
		![semiclassics]
		call MPI_Bcast( prefactF3	,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD,	ierr)
		!
		!derived scalars
		call MPI_Bcast( vol			,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast(	recpLatt	, 	dim**2	, 	MPI_DOUBLE_PRECISION	,	root, 	MPI_COMM_WORLD, ierr)
		!
		!
		if( myID /= root ) call allocateArrays()
		!
		!
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	shells		,	nShells	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)

		!call MPI_Bcast(	atPot_inp	,	nAt_inp	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	atPot		,	nAt		,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)


		!call MPI_Bcast(	proj_at_inp	,	nWfs_inp	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD,	ierr)
		!call MPI_Bcast( proj_nX_inp	,	nWfs_inp	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD,	ierr)
		!call MPI_Bcast( proj_nY_inp	,	nWfs_inp	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD,	ierr)

		call MPI_Bcast(	proj_at		,	nWfs	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast( proj_nX		,	nWfs	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast( proj_nY		,	nWfs	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD,	ierr)


		call MPI_BARRIER(MPI_COMM_WORLD, ierr)	!not necessary, helps against paranoia though
		!
		return
	end subroutine




	subroutine allocateArrays()
		!allocates all needed arrays
		!
		!atoms (super Cell)
		allocate(	relXpos(nAt)		)
		allocate(	relYpos(nAt)		)
		allocate(	atRx(nAt)			)
		allocate(	atRy(nAt)			)
		allocate(	atPot(nAt)			)
		allocate(	dVpot(nAt)			)
		allocate(	atPos(dim,nAt)		)
		allocate(	atR(dim,nAt) 		)
		!wann (super Cell)
		allocate(	proj_at(nWfs)		)
		allocate(	proj_nX(nWfs)		)
		allocate(	proj_nY(nWfs)		)
		!
		!
		!meshes
		allocate(	qpts(dim,nQ)		)
		allocate(	kpts(dim,nK)		)
		!w90
		allocate(	shells(nShells)		)
		!
		!
		return
	end subroutine





	logical function insideAt(at, r)
		!return true if real point is inside the potential of atom at
		integer,	intent(in)		:: at
		real(dp),	intent(in)		:: r(2)
		real(dp)					:: posX, posY, radX, radY, xRel, yRel
		logical						:: insideX, insideY
		!
		!PROJECT TO FIRST UNIT CELL
		xRel	=	dmod(	r(1),	aX)
		yRel	=	dmod(	r(2),	aY)
		!CHECK IF INSIDE
		insideX	= .false.
		insideY = .false.
		posX	= atPos(1,at)
		posY	= atPos(2,at)
		radX	= atR(1,at)
		radY	= atR(2,at)
		!
		if( posX-radX <= xRel .and. xRel <= posX+radX ) then
			insideX = .true.
		end if
		if( posY-radY <= yRel .and. yRel <= posY+radY) then
			insideY = .true.
		end if

		insideAt = insideX .and. insideY
		!
		return
	end function










!privat:
	subroutine qmeshGen()
		!generates the (coarse) k point mesh for solving electronic structure
		integer		:: qIx, qIy, qI
		real(dp)	:: qxMin, qyMin
		!
		if( myID == root ) then
			qxMin	= -1.0_dp * PI_dp * aY 	/ 		vol
			dqx		=  2.0_dp * PI_dp * aY  /	(vol * nQx)
			qyMin	= -1.0_dp * PI_dp * aX	/		vol
			dqy		=  2.0_dp * PI_dp * aX	/	(vol * nQy)
			!
			if( abs(dqx - dqy) > 1e-7_dp) then
				write(*,*)	'[qmeshGen]: non uniform mesh; nQx=',nQx,', nQy=',nQy
				stop 	'q mesh setup failed'
			end if


			do qIy = 1, nQy
				do qIx = 1, nQx
					qI	=	(qIy-1) * nQx + qIx
					qpts(1,qI)	=	qxMin + (qIx-1) * dqx  		!x component
					qpts(2,qI)	=	qyMin + (qIy-1) * dqy  		!y component
				end do
			end do
		end if
		call MPI_Bcast(qpts, dim*nQ, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
		!
		return
	end subroutine


	subroutine kWmeshGen()
		!generates the (fine) k point mesh for wannier interpolation
		integer		:: kIx, kIy, kI
		real(dp)	:: kxMin, kyMin
		!
		if( myID == root ) then
			kxMin	= -1.0_dp * PI_dp * aY / vol
			dkx		=  2.0_dp * PI_dp * aY / (vol * (nKx-1))
			kyMin	= -1.0_dp * PI_dp * aX / vol
			dky		=  2.0_dp * PI_dp * aX / (vol* (nKy-1))
			!
			do kIy = 1, nKy
				do kIx = 1, nKx
					kI	=	(kIy-1) * nKx + kIx
					kpts(1,kI)	=	kxMin + (kIx-1) * dkx  		!x component
					kpts(2,kI)	=	kyMin + (kIy-1) * dky  		!y component
				end do
			end do
		end if
		!
		call MPI_Bcast(kpts, dim*nK, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
		!
		return
	end subroutine


	!subroutine rmeshGen()
	!	!generates the real space mash
	!	!meshes from (0,0) to (aX*nSC,0) (0,aY*nSC) (aX*nSC,aY*nSC)
	!	integer		:: rIx, rIy, rIz, rI
	!	real(dp)	:: rxMin, ryMin, rzMin, dx, dy, dz
	!	!
	!	if( myID == root ) then
	!		rxMin	= 0.0_dp
	!		dx		= aX / real(nRx-1,dp)
	!		ryMin	= 0.0_dp
	!		dy		= aY / real(nRy-1,dp)
	!		rzMin	= 0.0_dp
	!		dz		= 0.0_dp
	!		if( nRz > 1 )	dz		= min(aX,aY) * 0.99_dp / real(nRz-1,dp)
	!		!
	!		do rIz = 1, nRz
	!			do rIy = 1, nRy
	!				do rIx = 1, nRx
	!					rI	=	 (	(rIz-1) * nRy +(rIy-1) ) * nRx + rIx
	!					rpts(1,rI)	= rxMin + (rIx-1) * dx		!x component
	!					rpts(2,rI)	= ryMin + (rIy-1) * dy		!y component
	!					rpts(3,rI)	= rzMin + (riZ-1) * dz
	!				end do
	!			end do
	!		end do
	!	end if
	!	call MPI_Bcast(rpts, size(rpts,1)*nR, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
	!	!
	!	return
	!end subroutine

!G SET
	subroutine getTestGridSize(nGx, nGy)
		integer,	intent(out)			:: nGx, nGy
		!
		!old:
		!nGx	= ceiling(	aX		 *	( Gcut + 2.0_dp	)	/ PI_dp		)
		!nGy = ceiling(	aY		 *	( Gcut + 2.0_dp	)	/ PI_dp		)
		!	
		!
		!new:
		!		nGi * 
		!
		nGx 	= ceiling(  (aX*Gcut/PI_dp + 1.0_dp)		) + 1
		nGy 	= ceiling(  (aY*Gcut/PI_dp + 1.0_dp)		) + 1

		!
		!make sure Grid is symmetric (needs to be odd number)
		if(mod(nGx,2)==0) nGx = nGx + 1 
		if(mod(nGy,2)==0) nGy = nGy + 1
		!
		write(*,*)	"[getTestGridSize]: test grid should contain",nGx*nGy," basis functions"
		!
		return
	end subroutine





!
	subroutine popGvec()
		integer						:: qi, ix, iy, inside,tot, qLoc, Gmin
		real(dp)					:: kg(2), Gtest(2)
		real(dp),	allocatable		:: Gtemp(:,:,:)
		character(len=20)			::	filename
		character(len=1024)			::	format='(a,i7.7)'
		!
		allocate(	nGq(					qChunk		)		)
		allocate(	Gtemp(	dim,	nG ,	qChunk		)		)
		!
		qLoc = 1
		do qi = myID*qChunk+1, myID*qChunk+qChunk
			nGq(qLoc)	= 0
			inside 		= 0
			tot 		= 0
			!LOOP TEST GRID
			do ix = 1, nGx
				do iy = 1, nGy
					!GET CURRENT GRID POINT
					Gtest(:)	= (ix-1-nGx/2) * recpLatt(:,1) + (iy-1-nGy/2) * recpLatt(:,2)
					kg(:)	= qpts(:,qi) + Gtest(:)
					tot		= tot + 1
					!
					!IF IN SPHERE APPEND TO GVEC
					if( norm2(kg) < Gcut ) then
						nGq(qLoc) = nGq(qLoc) + 1
						Gtemp(:,nGq(qLoc),qLoc) = kg(:)
						inside = inside + 1
						if( ix==1 .or. ix==nGx ) 	stop	"[;popGvec]: ERROR, hit Gx boundary! Gcut-sphere hits boundary of test Grid. Developer should increase test grid size"
						if( iy==1 .or. iy==nGy )	stop	"[;popGvec]: ERROR, hit Gy boundary! Gcut-sphere hits boundary of test Grid. Developer should increase test grid size"
					end if
				end do
			end do
			!DEBUG INFO
			if(nGq(qLoc) > nG) write(*,'(a,i3,a,i4,a,i6)')	"[#",myID,&
											";popGvec]: WARNING, somehow counted more basis functions at qi=",qi," limit nG=",nG	
			qLoc = qLoc + 1	
		end do
		!
		!DEBUG OUTPUT
		Gmax = maxval(nGq(:))
		Gmin = minval(nGq(:))
		write(*,'(a,i3,a,i3,a,i7,a,i7)')	"[#",myID,";popGvec]: my qChunk=",qChunk, " Gmax=",Gmax,";	  Gmin=",Gmin
		call MPI_ALLREDUCE(Gmax, GmaxGLOBAL, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr) 
		call MPI_ALLREDUCE(Gmin, GminGLOBAL, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
		if( myID == root ) write(*,'(a,i3,a,i7,a,i7)') "[#",myID,";popGvec]: global Gmax=",GmaxGLOBAL, ";	global Gmin=",GminGLOBAL
		!
		!
		!WRITE TO FILE
		do qLoc = 1, qChunk
			!get global q index & according filename
			qi = myID*qChunk + qLoc
			write(filename, format) raw_dir//'gVec.',qi
			!
			!write to file
			open(unit=210, file=filename		, form='unformatted', access='stream', action='write',status='replace') 
			write(210)	Gtemp(1:dim,1:GmaxGLOBAL,qLoc)
			close(210)
		end do
		!
		!
		return
	end subroutine




	subroutine testG(a1, a2, b1, b2)
		real(dp),		intent(in)		:: a1(2), a2(2), b1(2), b2(2)

		if( abs( dot_product(a1,b1)  - 2.0_dp*PI_dp) > acc   ) then
			write(*,*)"[testG]: problem with reciprocal lattice setup  (cond1 not fullfilled)"
		end if

		if( abs( dot_product(a2,b2)  - 2.0_dp*PI_dp) > acc   ) then
			write(*,*)"[testG]: problem with reciprocal lattice setup  (cond2 not fullfilled)"
		end if

		if( abs( dot_product(a1,b2)  ) > acc   ) then
			write(*,*)"[testG]: problem with reciprocal lattice setup  (cond1 not fullfilled)"
		end if

		if( abs( dot_product(a2,b1)  ) > acc   ) then
			write(*,*)"[testG]: problem with reciprocal lattice setup  (cond1 not fullfilled)"
		end if


		return
	end subroutine




	subroutine popAtPos()
		!populates the atom position vector
		integer		:: at
		!
		if(myID == root ) then
			

			do at = 1, nAt
				!test if inside unit cell
				if( relXpos(at) < 0.0_dp .or. relXpos(at) > 1.0_dp ) stop '[popAtPos]: relXpos out of bounds'
				if( relYpos(at) < 0.0_dp .or. relYpos(at) > 1.0_dp ) stop '[popAtPos]: relYpos out of bounds'
				atPos(1,at)	= relXpos(at) * aX	!x component
				atPos(2,at)	= relYpos(at) * aY	!y component
			end do




		end if
		call MPI_Bcast(atPos, dim*nAt, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
		!
		return
	end subroutine


	subroutine popAtR()
		!populates the atom radius vector
		integer		:: at
		!
		if(myID == root ) then
			do at = 1, nAt
				atR(1,at)	= atRx(at)
				atR(2,at)	= atRy(at)

				!Test if well exceeds unit cell
				if( (atPos(1,at)-atR(1,at)) < 0.0_dp .or. (atPos(1,at)+atR(1,at)) > aX ) stop'[popAtR]: wells exceed unit cell (x-dir)'
				if( (atPos(2,at)-atR(2,at)) < 0.0_dp .or. (atPos(2,at)+atR(2,at)) > aY ) stop'[popAtR]: wells exceed unit cell (x-dir)'
			end do
			!ToDo: test if wells overlap

		end if
		call MPI_Bcast(atR, dim*nAt, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
		!
		return
	end subroutine








end module util_sysPara



















!	subroutine calcRcell()
!		integer		:: nI, nJ, n
!		!
!		R0 = 1
!		if(myID == root ) then
!			do nJ = 1, nSCy
!				do nI = 1, nSCx
!					n = (nJ-1) * nSCx + nI  !rI	=	(rIy-1) * nRy + rIx
!					Rcell(1,n)	= real((nI-1),dp) * aX -real(nSCx-1,dp) * aX / 2.0_dp
!					Rcell(2,n)	= real((nJ-1),dp) * aY -real(nSCy-1,dp) * aY / 2.0_dp
!					if(abs(Rcell(1,n))< machineP .and. abs(Rcell(2,n))< machineP ) R0 = n
!				end do
!			end do
!		end if
!		call MPI_Bcast(Rcell, dim*nSC, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
!		!
!		return
!	end subroutine
