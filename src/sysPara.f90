module sysPara
	!this modules reads the input file and generates the meshes accordingly
	use mpi
	use mathematics, only: dp, PI_dp, setAcc, acc, machineP
	use m_config
	implicit none

	private
	public :: 	readInp, readGvec, insideAt, getRindex, getRleftX, getRrightX, getRleftY, getRrightY,& 
				getKindex, getGammaPoint, getPot, &
				dim, aX, aY, vol, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot, &
				nG, nGq, nG0, Gmax, GmaxGLOBAL, Gcut, Gvec, Gtest, R0, nSolve, &
				nQ, nQx, nQy, nKx, nKy, nK, nSC, nSCx, nSCy, dqx, dqy, dkx, dky, &
				nR, nRx, nRy,  dx, dy,nw90it, shell, &
				nBands, nWfs,   &
				atPos, atR, qpts, rpts, Rcell, kpts, Zion, recpLatt, &
				Bext, prefactF3, &
				seedName, w90_dir, info_dir, mkdir, raw_dir,&
				debugProj, debugHam, debugWann, doSolveHam, doMagHam, useBloch, doPw90, pw90GaugeB, doVdesc,  &
				doBerry, useRot, doWanni, doVeloNUM, doNiu, doPei, doGaugBack, writeBin, &
				myID, nProcs, root, ierr, qChunk


	!
	integer  										:: 	dim=2, nAt=0, nG, nGdim, Gmax, GmaxGLOBAL,nSolve=20, nG0,&  
														nQx=1, nQy=1,nQ , nSCx=1, nSCy=1,& 
														nKx=1, nKy=1, nK, R0,  &
														nw90it, shell, &
														nRx=10, nRy=10, nR, nBands=1,nWfs=1, nSC, &
														myID, nProcs, root, ierr, qChunk
	real(dp) 										::	aX=0.0_dp, aY=0.0_dp,vol=0.0_dp, Gcut=2*PI_dp, thres,& 
														dx, dy, dqx, dqy, dkx, dky, B0, Bext(3)	, prefactF3, recpLatt(2,2)
	character(len=3)								::	seedName										
	character(len=9)								::	w90_dir	="w90files/"
	character(len=7)								::	info_dir="output/"
	character(len=8)								::	raw_dir	="rawData/", mkdir="mkdir ./"	!to use with system(mkdir//$dir_path) 
	integer,	allocatable,	dimension(:)		::	nGq
	real(dp),	allocatable,	dimension(:)		::	relXpos, relYpos, atRx, atRy, atPot, dVpot, Zion
	real(dp),	allocatable,	dimension(:,:)		::	Gtest , atPos, atR, qpts, rpts, Rcell, kpts 

	real(dp),	allocatable,	dimension(:,:,:)	::	Gvec
	logical											::	debugHam, debugWann, debugProj, &
														doSolveHam, doMagHam, doPw90, pw90GaugeB, useBloch, doVdesc , &
														doBerry, useRot, doWanni, doVeloNUM, doNiu, doPei, doGaugBack, &
														writeBin 


	contains
!public:
	subroutine readInp()
		!root reades input parameters from input.txt file 
		!calls mesh generation subroutines
		!and bcasts everything around
		!the first array index is always the x,y value of the vector
		logical				:: dir_exists
		!
		dim = 	2
		!
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
		call rmeshGen()
		call popGtest()
		call popGvec()
		call popAtPos()
		call popAtR()
		call calcRcell()
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


	subroutine readGvec()
		!
		open(unit=800, file=raw_dir//"nGq.dat",form='unformatted',access='stream',action='read')
		read(800) nGq
		close(800)
		write(*,'(a,i3,a)')	"[#",myID,";readGvec]: read nGq"
		!
		!
		open(unit=805, file=raw_dir//"Gvec.dat",form='unformatted',access='stream',action='read')
		read(805) Gvec
		close(805)
		write(*,'(a,i3,a)')	"[#",myID,";readGvec]: read Gvec"
		!
		!
		return
	end subroutine



	
!READ & DISTRIBUTION ROUTINES
	subroutine rootRead()
		type(CFG_t) :: my_cfg
		
		!OPEN FILE
		call CFG_read_file(my_cfg,"input.txt")
		
		!READ SCALARS
		![unitCell]
		call CFG_add_get(my_cfg,	"unitCell%aX"      	,	aX  	   	,	"length of unit cell in agnstroem"		)
		call CFG_add_get(my_cfg,	"unitCell%aY"      	,	aY  	   	,	"length of unit cell in agnstroem"		)
		![atoms]
		call CFG_add_get(my_cfg,	"atoms%doVdesc"		,	doVdesc		,	"switch on/off linear descending pot"	)
		call CFG_add_get(my_cfg,	"atoms%nAt"			,	nAt			,	"number of atoms per unit cell"			)
		![wann]
		call CFG_add_get(my_cfg,	"wann%nBands"		,	nBands	 	,	"# of bands to project onto trial orbs"	)
		call CFG_add_get(my_cfg,	"wann%nWfs"			,	nWfs	 	,	"# wannier functions to generate"		)
		![perturbation]
		call CFG_add_get(my_cfg,	"perturbation%B0"	,	B0			,	"scaling fact. of ext. magnetic field"	)
		call CFG_add_get(my_cfg,	"perturbation%Bext"	,	Bext		,	"vector of ext. magnetic field"			)
		![numerics]
		call CFG_add_get(my_cfg,	"numerics%Gcut"		,	Gcut	    ,	"k space cut of parameter"				)
		call CFG_add_get(my_cfg,	"numerics%nSolve"	,	nSolve	    ,	"number of eigenstates to find"			)
		call CFG_add_get(my_cfg,	"numerics%nQx"     	,	nQx      	,	"amount of k points used"				)
		call CFG_add_get(my_cfg,	"numerics%nQy"     	,	nQy      	,	"amount of k points used"				)
		call CFG_add_get(my_cfg,	"numerics%nSCx"     ,	nSCx   		,	"#			supercells "				)
		call CFG_add_get(my_cfg,	"numerics%nSCy"     ,	nSCy   		,	"#			supercells "				)
		call CFG_add_get(my_cfg,	"numerics%nKx"		,	nKx			,	"# k x points of interpolation mesh"	)
		call CFG_add_get(my_cfg,	"numerics%nKy"		,	nKy			,	"# k x points of interpolation mesh"	)
		call CFG_add_get(my_cfg,	"numerics%nRx"     	,	nRx      	,	"amount of r points used"				)
		call CFG_add_get(my_cfg,	"numerics%nRy"     	,	nRy      	,	"amount of r points used"				)
		call CFG_add_get(my_cfg,	"numerics%thres"    ,	thres      	,	"threshold for overlap warnings"		)
		![methods]
		call CFG_add_get(my_cfg,	"methods%doSolveHam",	doSolveHam	,	"solve electronic structure or read in"	)
		call CFG_add_get(my_cfg,	"methods%doMagHam"	,	doMagHam	,	"include B-field via peierls in ham."	)
		call CFG_add_get(my_cfg,	"methods%useBloch"	,	useBloch	,	"use bloch phase for projections	"	)
		call CFG_add_get(my_cfg,	"methods%doPw90"	,	doPw90		,	"read in the matrices in wann base	"	)	
		call CFG_add_get(my_cfg,	"methods%doBerry"	,	doBerry		,	"switch on/off 	berry( unk) method "	)
		call CFG_add_get(my_cfg,	"methods%writeBin"	,	writeBin	,	"switch for writing binary files"		)
		![pw90]
		call CFG_add_get(my_cfg,	"pw90%seedName"		, 	 seedName	,	"seedName for wannier files(char len=3)")
		call CFG_add_get(my_cfg,	"pw90%shell"		, 	 shell		,	"manually set the shell to use for FD  ")
		call CFG_add_get(my_cfg,	"pw90%nw90it"		, 	 nw90it		,	"number of iterations for wannnierisat,")
		call CFG_add_get(my_cfg,	"pw90%pw90GaugeB"	,	pw90GaugeB	,	"logical for switching gauge trafo	   ")
		![berry]
		call CFG_add_get(my_cfg,	"berry%useRot"		,	useRot		,	"logical for switching on/off rotation	")
		call CFG_add_get(my_cfg,	"berry%doVeloNUM"	,	doVeloNUM	,	"if true tb velocities, else analyitcal")
		call CFG_add_get(my_cfg,	"berry%doNiu"		,	doNiu		,	"switch for nius first order pol"		)
		call CFG_add_get(my_cfg,	"berry%doPei"		,	doPei		,	"switch for  peierls first order pol"	)
		call CFG_add_get(my_cfg,	"berry%doWanni"		,	doWanni		,	"switch on/off 	wannier( wnf ) method"	)
		call CFG_add_get(my_cfg,	"berry%doGaugBack"	,	doGaugBack	,	"switch for trafo: Wann to Ham gauge"	)
		
		![semiclassics]
		call CFG_add_get(my_cfg,	"semiclassics%prefactF3"	,	prefactF3,	"real prefactor for F3 "			)
		![debug]
		call CFG_add_get(my_cfg,	"debug%debugProj"	, 	debugProj	,	"switch for debuging tests in solveHam"	)
		call CFG_add_get(my_cfg,	"debug%debugHam"	, 	debugHam	,	"switch for debuging tests in solveHam"	)
		call CFG_add_get(my_cfg,	"debug%debugWann"	, 	debugWann	,	"switch for debuging in wannier"		)
	
		!SET
		nGdim = getTestGridSize()	
		nG	= nGdim**2
		vol	=	aX 		* 	aY
		nR 	= 	nRx 	*	nRy
		nQ 	= 	nQx 	*	nQy
		nSC =	nSCx	*	nSCy
		nK =	nKx		*	nKy
		Bext=	B0 		* 	Bext
		!
		recpLatt		= 0.0_dp
		recpLatt(1,1)	= 2.0_dp * PI_dp * aY / vol
		recpLatt(2,2)	= 2.0_dp * PI_dp * aX / vol
		!
		call allocateArrays()
		!
		![atoms]
		call CFG_add_get(my_cfg,	"atoms%relXpos"		,	relXpos		,	"relative positions in unit cell"		)
		call CFG_add_get(my_cfg,	"atoms%relYpos"		,	relYpos		,	"relative positions in unit cell"		)
		call CFG_add_get(my_cfg,	"atoms%atRx"		,	atRx		,	"radius of each atom in angstroem"		)
		call CFG_add_get(my_cfg,	"atoms%atRy"		,	atRy		,	"radius of each atom in angstroem"		)
		call CFG_add_get(my_cfg,	"atoms%atPot"		,	atPot		,	"potential depth in hartree"			)
		call CFG_add_get(my_cfg,	"atoms%dVpot"		,	dVpot		,	"potential gradient"					)
		call CFG_add_get(my_cfg,	"atoms%Zion"		,	Zion		,	"effective charge of the ions"			)
		!
		!
		return
	end subroutine


	subroutine bcastPARAM()
		integer							:: ierr
		!
		!
		![unitCell]
		call MPI_Bcast( aX			,		1	,	MPI_DOUBLE_PRECISION	, 	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( aY			,		1	,	MPI_DOUBLE_PRECISION	, 	root,	MPI_COMM_WORLD, ierr)
		![atoms]
		call MPI_Bcast( doVdesc		, 		1	,		MPI_LOGICAL			,	root,	MPI_COMM_WORLD, ierr)	
		call MPI_Bcast( nAt			, 		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)		
		![wann]
		call MPI_Bcast( nBands		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)		
		call MPI_Bcast( nWfs		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		![perturbation]
		call MPI_Bcast(	B0 			,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)		
		call MPI_Bcast(	Bext 		,		3	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		![numerics]
		call MPI_Bcast(	Gcut		,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nSolve		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nQx			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nQy			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nSCx		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nSCy		,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nKx			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nKy			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nRx			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	nRy			,		1	,		MPI_INTEGER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	thres		,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		![methods]
		call MPI_Bcast( doSolveHam	,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	doMagHam	,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( useBloch	,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( doPw90		,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( doBerry		,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( writeBin	,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD, ierr)
		![pw90]
		call MPI_Bcast( seedName	,		3	, 	MPI_CHARACTER			,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( shell		,		1	,	MPI_INTEGER				,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( nw90it		,		1	,	MPI_INTEGER				,	root,	MPI_COMM_WORLD, ierr)		
		call MPI_Bcast( pw90GaugeB	,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD, ierr)
		![berry]
		call MPI_Bcast( useRot		,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast( doVeloNUM	,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast( doNiu		,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast( doPei		,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast( doWanni		,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD,	ierr)
		call MPI_Bcast( doGaugBack	,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD,	ierr)
		![semiclassics]
		call MPI_Bcast( prefactF3	,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD,	ierr)
		![debug]		
		call MPI_Bcast( debugProj	,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD, ierr)		
		call MPI_Bcast( debugHam	,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( debugWann	,		1	,	MPI_LOGICAL				,	root,	MPI_COMM_WORLD, ierr)		
		!
		!
		!derived scalars
		call MPI_Bcast( nGdim		,		1	,	MPI_INTEGER				,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( vol			,		1	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD,	ierr)
		!
		!
		if( myID /= root ) then
			nG 	= 	nGdim**2
			nR 	= 	nRx 	*	nRy
			nQ 	= 	nQx 	*	nQy
			nSC =	nSCx	*	nSCy
			nK 	=	nKx		*	nKy
			Bext=	B0 		* 	Bext
		end if
		!
		!
		if( myID /= root ) call allocateArrays()
		!
		return
	end subroutine




	subroutine allocateArrays()
		!allocates all needed arrays
		!
		!basis
		allocate(	Gtest(	dim,	nG				)		)
		!atoms
		allocate(	relXpos(nAt)		)
		allocate(	relYpos(nAt)		)
		allocate(	atRx(nAt)			)
		allocate(	atRy(nAt)			)
		allocate(	atPot(nAt)			)
		allocate(	dVpot(nAt)			)
		allocate(	atPos(dim,nAt)		)
		allocate(	atR(dim,nAt) 		)
		allocate(	Zion(nAt)			)
		!meshes
		allocate(	qpts(dim,nQ)		)
		allocate(	rpts(dim,nR)		)
		allocate(	Rcell(dim,nSC)		)
		allocate(	kpts(dim,nK)		)
		!
		!
		return
	end subroutine



	subroutine bcastArrays()
		!basis
		call MPI_Bcast(	Gtest		, dim*nG	,	MPI_DOUBLE_PRECISION	, 	root,	MPI_COMM_WORLD, ierr)
		!atoms
		call MPI_Bcast(	relXpos		, nAt 		,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	relYpos		, nAt		,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	atRx		, nAt		,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	atRy		, nAt		,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	atPot		, nAt		,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	dVpot		, nAt		,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	atPos		, dim*nAt 	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	atR			, dim*nAt	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	Zion		, nAt		,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		!meshes
		call MPI_Bcast( qpts		, dim*nQ	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)		
		call MPI_Bcast( rpts		, dim*nR	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast( Rcell		, dim*nSc	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)
		call MPI_Bcast(	kpts		, dim*nK	,	MPI_DOUBLE_PRECISION	,	root,	MPI_COMM_WORLD, ierr)		
		!
		if(myID == root ) write(*,*) "finished broadcasting input parameters"
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


	integer function getRindex(xi,yi)
		integer,	intent(in)		:: xi, yi
		!
		getRindex = (yi-1) * nRx + xi
		return
	end function



	integer function getRleftX(xi,yi)
		integer,		intent(in)		:: xi, yi
		!
		if( xi == 1) then
			getRleftX = getRindex(nRx,yi)
		else
			getRleftX = getRindex(xi-1,yi)
		end if	 
		!
		return
	end function


	integer function getRrightX(xi,yi)
		integer,		intent(in)		:: xi, yi
		!
		if( xi == nRx) then
			getRrightX = getRindex(1,yi)
		else
			getRrightX = getRindex(xi+1,yi)
		end if	 
		!
		return
	end function



	integer function getRleftY(xi,yi)
		integer,		intent(in)		:: xi, yi
		!
		if( yi == 1) then
			getRleftY = getRindex(xi,nRy)
		else
			getRleftY= getRindex(xi,yi-1)
		end if	 
		!
		return
	end function


	integer function getRrightY(xi,yi)
		integer,		intent(in)		:: xi, yi
		!
		if( yi == nRy) then
			getRrightY = getRindex(xi,1)
		else
			getRrightY = getRindex(xi,yi+1)
		end if	 
		!
		return
	end function

	integer function getKindex(qx,qy)
		integer,	intent(in)		:: qx, qy
		!
		getKindex = (qy-1) * nQx + qx
		return
	end function


	integer function getGammaPoint()
		getGammaPoint	= getKindex(1+nQx/2,1+nQy/2)
		return
	end function


	real(dp) function getPot(ri)
		integer,	intent(in)	:: ri
		integer					:: at
		logical					:: found
		found		= .false.
		getPot 	= 0.0_dp
		at			= 1
		do while(at <= nAt .and. .not. found )
			if( insideAt(at,rpts(:,ri)) ) then
				getPot	= atPot(at)
				found		= .true.
			end if
			at = at +1 
		end do
		!
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
			qxMin	= -1.0_dp * PI_dp * aX 	/ 		vol
			dqx		=  2.0_dp * PI_dp * aX  /	(vol * nQx)
			qyMin	= -1.0_dp * PI_dp * aY	/		vol
			dqy		=  2.0_dp * PI_dp * aY	/	(vol * nQy)
			!
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
			kxMin	= -1.0_dp * PI_dp * aX / vol
			dkx		=  2.0_dp * PI_dp * aX / (vol * (nKx-1))
			kyMin	= -1.0_dp * PI_dp * aY / vol
			dky		=  2.0_dp * PI_dp * aY / (vol* (nKy-1))
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


	subroutine rmeshGen()
		!generates the real space mash
		!meshes from (0,0) to (aX*nSC,0) (0,aY*nSC) (aX*nSC,aY*nSC)
		integer		:: rIx, rIy, rI
		real(dp)	:: rxMin, ryMin
		!
		if( myID == root ) then
			rxMin	=  -real(nSCx-1,dp) * aX / 2.0_dp
			dx		= real(nSCx,dp) * aX / real(nRx-1,dp)
			ryMin	= -real(nSCy-1,dp) * aX / 2.0_dp
			dy		= real(nSCy,dp) * aY / real(nRy-1,dp)
			!
			do rIy = 1, nRy
				do rIx = 1, nRx
					rI	=	(rIy-1) * nRx + rIx
					rpts(1,rI)	= rxMin + (rIx-1) * dx		!x component
					rpts(2,rI)	= ryMin + (rIy-1) * dy		!y component
				end do
			end do
		end if
		call MPI_Bcast(rpts, dim*nR, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
		!
		return
	end subroutine

!G SET
	integer function getTestGridSize()
		integer				:: nGrid
		!
		!nGrid = ceiling( 	dmax1(aX,aY)*Gcut/PI_dp 	+	dsqrt(2.0_dp)		)

		!the addition of 1 makes shure that k point shifts are handled
		nGrid = ceiling(	dmax1(aX,aY) *	( Gcut + 1.0_dp	)	 / PI_dp		)
		!
		!make sure Grid is symmetric (needs to be odd number)
		if(mod(nGrid,2)==0) nGrid = nGrid + 1 
		!
		getTestGridSize = nGrid
		write(*,*)	"[getTestGridSize]: test grid should contain",getTestGridSize**2," basis functions"
		!
		return
	end function



	subroutine popGtest()
	!populates the G vector (basis vector)
		integer		:: i, ix, iy
		real(dp)	:: thres, b1(2), b2(2), a1(2), a2(2)
		!
		thres	= 1e-15_dp
		
		
		a1(1)	= aX
		a1(2)	= 0.0_dp
		a2(1)	= 0.0_dp
		a2(2)	= aY


		b1(1)	= 2.0_dp * PI_dp * aY / vol
		b1(2)	= 0.0_dp
		b2(1)	= 0.0_dp
		b2(2)	= 2.0_dp * PI_dp * aX / vol

		if( myID == root ) then
			call testG( a1, a2, b1, b2)
			do ix = 1, nGdim
				do iy = 1, nGdim
					i	= (iy-1) * nGdim + ix
					Gtest(:,i)	= (ix-1-nGdim/2) * recpLatt(:,1) + (iy-1-nGdim/2) * recpLatt(:,2)
				end do
			end do
		end if
		call MPI_Bcast(Gtest, size(Gtest,1)*size(Gtest,2), MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
		!
		return
	end subroutine


!
	subroutine popGvec()
		integer						:: qi, gi, inside,tot, qLoc
		real(dp)					:: kg(2)
		!^
		!fill array
		write(*,'(a,i3,a,i6)')	"[#",myID,";popGvec]: my qChunk=",qChunk
		allocate(	nGq(					qChunk		)		)
		allocate(	Gvec(	dim,	nG ,	qChunk		)		)
		!
		qLoc = 1
		do qi = myID*qChunk+1, myID*qChunk+qChunk
			nGq(qLoc)	= 0
			inside 		= 0
			tot 		= 0
			do gi = 1, nG
				kg(:)	= qpts(:,qi) + Gtest(:,gi)
				tot		= tot + 1
				if( norm2(kg) < Gcut ) then
					nGq(qLoc) = nGq(qLoc) + 1
					Gvec(:,nGq(qLoc),qLoc) = kg(:)
					inside = inside + 1
					if( gi == 1)	write(*,*)	"[popGvec]: warning hit boundary of Gtest grid"
				end if
			end do
			!DEBUG INFO
			if(nGq(qLoc) > nG) write(*,'(a,i3,a,i4,a,i6)')	"[#",myID,&
											";popGvec]: warning, somehow counted more basis functions at qi=",qi," limit nG=",nG	
			write(*,'(a,i3,a,i6,a,i4)')	"[#",myID,";popGvec]: using ",nGq(qLoc), "basis functions at qi=",qi
			qLoc = qLoc + 1	
		end do
		!
		Gmax = maxval(nGq)
		call MPI_ALLREDUCE(Gmax, GmaxGLOBAL, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
		write(*,'(a,i3,a,i7)')	"[#",myID,";popGvec]: maximum amount of basis functions is",Gmax
		call MPI_BARRIER( MPI_COMM_WORLD, ierr)
		if( myID == root ) write(*,'(a,i3,a,i7)') "[#",myID,";popGvec]: global Gmax=",GmaxGLOBAL
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
			end do
		end if
		call MPI_Bcast(atR, dim*nAt, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
		!
		return
	end subroutine


	subroutine calcRcell()
		integer		:: nI, nJ, n
		!
		R0 = 1
		if(myID == root ) then
			do nJ = 1, nSCy
				do nI = 1, nSCx
					n = (nJ-1) * nSCx + nI  !rI	=	(rIy-1) * nRy + rIx
					Rcell(1,n)	= real((nI-1),dp) * aX -real(nSCx-1,dp) * aX / 2.0_dp
					Rcell(2,n)	= real((nJ-1),dp) * aY -real(nSCy-1,dp) * aY / 2.0_dp
					if(abs(Rcell(1,n))< machineP .and. abs(Rcell(2,n))< machineP ) R0 = n
				end do
			end do
		end if
		call MPI_Bcast(Rcell, dim*nSC, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
		!
		return
	end subroutine








end module sysPara