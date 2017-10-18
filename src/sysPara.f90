module sysPara
	!this modules reads the input file and generates the meshes accordingly
	use mathematics, only: dp, PI_dp, setAcc, acc
	use m_config
	implicit none
	private
	public :: 	readInp, insideAt, getRindex, getRleftX, getRrightX, getRleftY, getRrightY,& 
				getKindex, getGammaPoint, getPot, &
				dim, aX, aY, vol, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot, &
				nG, nGdim, nG0, Gcut, nQ, nQx, nQy, nKx, nKy, nK, nSC, nSCx, nSCy, nR, nRx, nRy, R0,  dx, dy, dqx, dqy, dkx, dky, &
				gaugeSwitch, nBands, nWfs, connSwitch,  &
				Gvec, atPos, atR, qpts, rpts, Rcell, kpts, trialOrbVAL, trialOrbSw, Zion, &
				Bext, &
				debugProj, debugHam, debugWann, doSolveHam, doVdesc, doProj, doBerry, doWanni, doNiu, doPei, doGaugBack, writeBin


	!
	integer  										:: 	dim=2, nAt=0, nG, nGdim=16, nG0,  nQx=1, nQy=1,nQ , nSCx=1, nSCy=1,& 
														nKx=1, nKy=1, nK, connSwitch=0, &
														nRx=10, nRy=10, nR, R0=1, nBands=1,nWfs=1, nSC, gaugeSwitch, trialOrbSw
	real(dp) 										::	aX=0.0_dp, aY=0.0_dp,vol=0.0_dp, Gcut=2*PI_dp, thres,& 
														dx, dy, dqx, dqy, dkx, dky, B0, Bext(3)											
	real(dp),	allocatable,	dimension(:)		::	relXpos, relYpos, atRx, atRy, atPot, dVpot, trialOrbVAL, Zion
	real(dp),	allocatable,	dimension(:,:)		::	Gvec, atPos, atR, qpts, rpts, Rcell, kpts 
	logical											::	debugHam, debugWann, debugProj, &
														doSolveHam, doVdesc , doProj , &
														doBerry, doWanni, doNiu, doPei, doGaugBack, &
														writeBin















	contains
!public:
	subroutine readInp()
		!reades input parameters from input.txt file
		!calls mesh generation subroutines
		!the first array index is always the x,y value of the vector
		!
		type(CFG_t) :: my_cfg
		
		!read in config
		call CFG_read_file(my_cfg,"input.txt")
		
		![unitCell]
		call CFG_add_get(my_cfg,	"unitCell%aX"      	,	aX  	   	,	"length of unit cell in agnstroem"		)
		call CFG_add_get(my_cfg,	"unitCell%aY"      	,	aY  	   	,	"length of unit cell in agnstroem"		)
		![atoms]
		call CFG_add_get(my_cfg,	"atoms%nAt"			,	nAt			,	"number of atoms per unit cell"			)
		![numerics]
		call CFG_add_get(my_cfg,	"numerics%nGdim"    ,	nGdim	    ,	"amount of G_n used"					)
		call CFG_add_get(my_cfg,	"numerics%Gcut"		,	Gcut	    ,	"k space cut of parameter"				)
		call CFG_add_get(my_cfg,	"numerics%nQx"     	,	nQx      	,	"amount of k points used"				)
		call CFG_add_get(my_cfg,	"numerics%nQy"     	,	nQy      	,	"amount of k points used"				)
		call CFG_add_get(my_cfg,	"numerics%nSCx"     ,	nSCx   		,	"#			supercells "				)
		call CFG_add_get(my_cfg,	"numerics%nSCy"     ,	nSCy   		,	"#			supercells "				)
		call CFG_add_get(my_cfg,	"numerics%nRx"     	,	nRx      	,	"amount of r points used"				)
		call CFG_add_get(my_cfg,	"numerics%nRy"     	,	nRy      	,	"amount of r points used"				)
		call CFG_add_get(my_cfg,	"numerics%thres"    ,	thres      	,	"threshold for overlap warnings"		)
		![methods]
		
		call CFG_add_get(my_cfg,	"methods%doSolveHam",	doSolveHam	,	"solve electronic structure or read in"	)
		call CFG_add_get(my_cfg,	"methods%doVdesc"	,	doVdesc		,	"switch on/off linear descending pot"	)
		call CFG_add_get(my_cfg,	"methods%doProj"	,	doProj		,	"switch on/off 	projections onto trial"	)
		call CFG_add_get(my_cfg,	"methods%doBerry"	,	doBerry		,	"switch on/off 	berry( unk) method "	)
		call CFG_add_get(my_cfg,	"methods%doWanni"	,	doWanni		,	"switch on/off 	wannier( wnf ) method"	)
		call CFG_add_get(my_cfg,	"methods%doNiu"		,	doNiu		,	"switch for nius first order pol"		)
		call CFG_add_get(my_cfg,	"methods%doPei"		,	doPei		,	"switch for  peierls first order pol"	)
		call CFG_add_get(my_cfg,	"methods%doGaugBack",	doGaugBack	,	"switch for trafo: Wann to Ham gauge"	)
		![output]
		call CFG_add_get(my_cfg,	"output%writeBin"	,	writeBin	,	"switch for writing binary files"		)
		![debug]
		call CFG_add_get(my_cfg,	"debug%debugProj"	, 	debugProj	,	"switch for debuging tests in solveHam"	)
		call CFG_add_get(my_cfg,	"debug%debugHam"	, 	debugHam	,	"switch for debuging tests in solveHam"	)
		call CFG_add_get(my_cfg,	"debug%debugWann"	, 	debugWann	,	"switch for debuging in wannier"		)
		![wannier]
		call CFG_add_get(my_cfg,	"wann%gaugeSwitch"	,	gaugeSwitch ,	"switch gauge the basis coeff directly"	)
		call CFG_add_get(my_cfg,	"wann%nBands"		,	nBands	 	,	"# of bands to project onto trial orbs"	)
		call CFG_add_get(my_cfg,	"wann%nWfs"			,	nWfs	 	,	"# wannier functions to generate"		)
		call CFG_add_get(my_cfg,	"wann%trialOrbSw"	,	trialOrbSw	,	"switch different trial orbitals"		)
		call CFG_add_get(my_cfg,	"wann%R0"			,	R0			,	"home unit cell used for wann cent calc")
		call CFG_add_get(my_cfg,	"wann%nKx"			,	nKx			,	"# k x points of interpolation mesh"	)
		call CFG_add_get(my_cfg,	"wann%nKy"			,	nKy			,	"# k x points of interpolation mesh"	)
		call CFG_add_get(my_cfg,	"wann%connSwitch"	,	connSwitch	,	"connection via K or via R space"		)
		
		![perturbation]
		call CFG_add_get(my_cfg,	"perturbation%B0"	,	B0			,	"scaling fact. of ext. magnetic field"	)
		call CFG_add_get(my_cfg,	"perturbation%Bext"	,	Bext		,	"vector of ext. magnetic field"			)

		
		dim = 	2
		nG	=	nGdim**2 !makes sure that bove dimensions get equal basis vectors
		vol	=	aX 		* 	aY
		nR 	= 	nRx 	*	nRy
		nQ 	= 	nQx 	*	nQy
		nSC =	nSCx	*	nSCy
		nK =	nKx		*	nKy
		Bext=	B0 		* 	Bext

		!basis
		allocate(	Gvec(dim,nG)		)
		!atoms
		allocate(	relXpos(nAt)		)
		allocate(	relYpos(nAt)		)
		allocate(	atRx(nAt)			)
		allocate(	atRy(nAt)			)
		allocate(	atPot(nAt)			)
		allocate(	dVpot(nAt)			)
		allocate(	atPos(dim,nAt)		)
		allocate(	atR(dim,nAt) 		)
		allocate(	trialOrbVAL(nAt)	)
		allocate(	Zion(nAt)			)
		!meshes
		allocate(	qpts(dim,nQ)		)
		allocate(	rpts(dim,nR)		)
		allocate(	Rcell(dim,nSC)		)
		allocate(	kpts(dim,nK)		)
		
		![atoms]
		call CFG_add_get(my_cfg,	"atoms%relXpos"		,	relXpos		,	"relative positions in unit cell"		)
		call CFG_add_get(my_cfg,	"atoms%relYpos"		,	relYpos		,	"relative positions in unit cell"		)
		call CFG_add_get(my_cfg,	"atoms%atRx"		,	atRx		,	"radius of each atom in angstroem"		)
		call CFG_add_get(my_cfg,	"atoms%atRy"		,	atRy		,	"radius of each atom in angstroem"		)
		call CFG_add_get(my_cfg,	"atoms%atPot"		,	atPot		,	"potential depth in hartree"			)
		call CFG_add_get(my_cfg,	"atoms%dVpot"		,	dVpot		,	"potential gradient"					)
		call CFG_add_get(my_cfg,	"atoms%Zion"		,	Zion		,	"effective charge of the ions"			)
		![wann]
		call CFG_add_get(my_cfg,	"wann%trialOrbVAL"	,	trialOrbVAL	,	"weight of trial orbital for each atom"	)

		!calculate desired quantities from input: aLatt,kpts,rpts(aLatt),gVec
		call setAcc(thres)
		call kmeshGen()
		call rmeshGen()
		call popGvec()
		call popAtPos()
		call popAtR()
		call calcRcell()
		call kWmeshGen()

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
	subroutine kmeshGen()
		!generates the (coarse) k point mesh for solving electronic structure
		integer		:: qIx, qIy, qI
		real(dp)	:: qxMin, qyMin
		!
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
		!
		return
	end subroutine


	subroutine kWmeshGen()
		!generates the (fine) k point mesh for wannier interpolation
		integer		:: kIx, kIy, kI
		real(dp)	:: kxMin, kyMin
		!
		kxMin	= -1.0_dp * PI_dp * aX / vol
		dkx		=  2.0_dp * PI_dp * aX / (vol * nKx)
		kyMin	= -1.0_dp * PI_dp * aY / vol
		dky		=  2.0_dp * PI_dp * aY / (vol* nKy)
		!
		do kIy = 1, nKy
			do kIx = 1, nKx
				kI	=	(kIy-1) * nKx + kIx
				kpts(1,kI)	=	kxMin + (kIx-1) * dkx  		!x component
				kpts(2,kI)	=	kyMin + (kIy-1) * dky  		!y component
			end do
		end do
		!
		return
	end subroutine


	subroutine rmeshGen()
		!generates the real space mash
		!meshes from (0,0) to (aX*nSC,0) (0,aY*nSC) (aX*nSC,aY*nSC)
		integer		:: rIx, rIy, rI
		real(dp)	:: rxMin, ryMin
		!
		rxMin	= 0.0_dp
		dx		= real(nSCx,dp) * aX / real(nRx,dp)
		ryMin	= 0.0_dp
		dy		= real(nSCy,dp) * aY / real(nRy,dp)
		!
		do rIy = 1, nRy
			do rIx = 1, nRx
				rI	=	(rIy-1) * nRx + rIx
				rpts(1,rI)	= rxMin + (rIx-1) * dx		!x component
				rpts(2,rI)	= ryMin + (rIy-1) * dy		!y component
			end do
		end do
		!
		return
	end subroutine


	subroutine popGvec()
		!populates the G vector (basis vector)
		integer		:: i, ix, iy, nxMin, nyMin
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

		call testG( a1, a2, b1, b2)

		do ix = 1, nGdim
			do iy = 1, nGdim
				i	= (iy-1) * nGdim + ix
				Gvec(:,i)	= (ix-1-nGdim/2) * b1(:) + (iy-1-nGdim/2) * b2(:)
			end do
		end do


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
		do at = 1, nAt
			atPos(1,at)	= relXpos(at) * aX	!x component
			atPos(2,at)	= relYpos(at) * aY	!y component
		end do
		!
		return
	end subroutine


	subroutine popAtR()
		!populates the atom radius vector
		integer		:: at
		!
		do at = 1, nAt
			atR(1,at)	= atRx(at)
			atR(2,at)	= atRy(at)
		end do
		!
		return
	end subroutine


	subroutine calcRcell()
		integer		:: nI, nJ, n
		!
		do nJ = 1, nSCy
			do nI = 1, nSCx
				n = (nJ-1) * nSCx + nI  !rI	=	(rIy-1) * nRy + rIx
				Rcell(1,n)	= real((nI-1),dp) * aX
				Rcell(2,n)	= real((nJ-1),dp) * aY
			end do
		end do

		return
	end subroutine




end module sysPara