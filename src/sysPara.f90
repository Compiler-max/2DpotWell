module sysPara
	!this modules reads the input file and generates the meshes accordingly
	use mathematics, only: dp, PI_dp
	use m_config
	implicit none
	private
	public :: 	readInp, insideAt, getRindex, getKindex, &
				dim, aX, aY, vol, nAt, relXpos, relYpos, atRx, atRy, atPot,&
				nG, nG0, Gcut, nK, nKx, nKy, nSC, nSCx, nSCy, nR, nRx, nRy, R0,  dx, dy, dkx, dky, &
				gaugeSwitch, nWfs, &
				Gvec, atPos, atR, kpts, rpts, Rcell, trialOrbVAL, trialOrbSw, Zion


	!
	integer  										:: 	dim=2, nAt=0, nG=11, nG0,  nKx=1, nKy=1,nK , nSCx=1, nSCy=1,& 
														nRx=10, nRy=10, nR, R0=1, nWfs=1, nSC, gaugeSwitch, trialOrbSw
	real(dp) 										::	aX=0.0_dp, aY=0.0_dp,vol=0.0_dp, Gcut=2*PI_dp,& 
														dx, dy, dkx, dky											
	real(dp),	allocatable,	dimension(:)		::	relXpos, relYpos, atRx, atRy, atPot, trialOrbVAL, Zion
	real(dp),	allocatable,	dimension(:,:)		::	Gvec, atPos, atR, kpts, rpts, Rcell 















	contains
!public:
	subroutine readInp()
		!reades input parameters from input.txt file
		!calls mesh generation subroutines
		!the first array index is always the x,y value of the vector
		!
		type(CFG_t) :: my_cfg
		integer 	:: i
		
		!read in config
		call CFG_read_file(my_cfg,"input.txt")
		
		![unitCell]
		call CFG_add_get(my_cfg,	"unitCell%aX"      	,	aX  	   	,	"length of unit cell in agnstroem"		)
		call CFG_add_get(my_cfg,	"unitCell%aY"      	,	aY  	   	,	"length of unit cell in agnstroem"		)
		![atoms]
		call CFG_add_get(my_cfg,	"atoms%nAt"			,	nAt			,	"number of atoms per unit cell"			)
		![numerics]
		call CFG_add_get(my_cfg,	"numerics%nG"       ,	nG	     	,	"amount of G_n used"					)
		call CFG_add_get(my_cfg,	"numerics%Gcut"		,	Gcut	    ,	"k space cut of parameter"				)
		call CFG_add_get(my_cfg,	"numerics%nKx"     	,	nKx      	,	"amount of k points used"				)
		call CFG_add_get(my_cfg,	"numerics%nKy"     	,	nKy      	,	"amount of k points used"				)
		call CFG_add_get(my_cfg,	"numerics%nSCx"     ,	nSCx   		,	"#			supercells "				)
		call CFG_add_get(my_cfg,	"numerics%nSCy"     ,	nSCy   		,	"#			supercells "				)
		call CFG_add_get(my_cfg,	"numerics%nRx"     	,	nRx      	,	"amount of r points used"				)
		call CFG_add_get(my_cfg,	"numerics%nRy"     	,	nRy      	,	"amount of r points used"				)
		![wannier]
		call CFG_add_get(my_cfg,	"wann%gaugeSwitch"	,	gaugeSwitch ,	"switch gauge the basis coeff directly"	)
		call CFG_add_get(my_cfg,	"wann%nWfs"			,	nWfs	 	,	"# wannier functions to generate"		)
		call CFG_add_get(my_cfg,	"wann%trialOrbSw"	,	trialOrbSw	,	"switch different trial orbitals"		)
		call CFG_add_get(my_cfg,	"wann%R0"			,	R0			,	"home unit cell used for wann cent calc")

		dim = 	2
		vol	=	aX * aY
		nG	=	floor(	sqrt(real(nG,dp))	)**2 !makes sure that bove dimensions get equal basis vectors
		nK 	= 	nKx 	*	nKy
		nR 	= 	nRx 	*	nRy
		nSC =	nSCx	*	nSCy
		!basis
		allocate(	Gvec(dim,nG)		)
		!atoms
		allocate(	relXpos(nAt)		)
		allocate(	relYpos(nAt)		)
		allocate(	atRx(nAt)			)
		allocate(	atRy(nAt)			)
		allocate(	atPot(nAt)			)
		allocate(	atPos(dim,nAt)		)
		allocate(	atR(dim,nAt) 		)
		allocate(	trialOrbVAL(nAt)	)
		allocate(	Zion(nAt)			)
		!meshes
		allocate(	kpts(dim,nK)		)
		allocate(	rpts(dim,nR)		)
		allocate(	Rcell(dim,nSC)		)
		
		![atoms]
		call CFG_add_get(my_cfg,	"atoms%relXpos"		,	relXpos		,	"relative positions in unit cell"		)
		call CFG_add_get(my_cfg,	"atoms%relYpos"		,	relYpos		,	"relative positions in unit cell"		)
		call CFG_add_get(my_cfg,	"atoms%atRx"		,	atRx		,	"radius of each atom in angstroem"		)
		call CFG_add_get(my_cfg,	"atoms%atRy"		,	atRy		,	"radius of each atom in angstroem"		)
		call CFG_add_get(my_cfg,	"atoms%atPot"		,	atPot		,	"potential depth in hartree"			)
		call CFG_add_get(my_cfg,	"atoms%Zion"		,	Zion		,	"effective charge of the ions"			)
		![wann]
		call CFG_add_get(my_cfg,	"wann%trialOrbVAL"	,	trialOrbVAL	,	"weight of trial orbital for each atom"	)

		!calculate desired quantities from input: aLatt,kpts,rpts(aLatt),gVec
		
		call kmeshGen()
		call rmeshGen()
		call popGvec()
		call popAtPos()
		call popAtR()
		call calcRcell()

		return
	end


	logical function insideAt(at, x, y)
		!return true if real point is inside the potential of atom at
		integer,	intent(in)		:: at
		real(dp),	intent(in)		:: x, y
		integer						:: i
		real(dp)					:: posX, posY, radX, radY
		logical						:: insideX, insideY
		!
		insideX	= .false.
		insideY = .false.
		posX	= atPos(1,at)
		posY	= atPos(2,at)
		radX	= atR(1,at)
		radY	= atR(2,at)
		!
		if( posX-radX <= x .and. x <= posX+radX ) then
			insideX = .true.
		end if
		if( posY-radY <= y .and. y <= posY+radY) then
			insideY = .true.
		end if

		insideAt = insideX .and. insideY
		!
		return
	end


	integer function getRindex(xi,yi)
		integer,	intent(in)		:: xi, yi
		!
		getRindex = (yi-1) * nRx + xi
		return
	end


	integer function getKindex(kx,ky)
		integer,	intent(in)		:: kx, ky
		!
		getKindex = (ky-1) * nKx + kx
		return
	end













!privat:
	subroutine kmeshGen()
		!generates the k point mesh
		integer		:: kIx, kIy, kI
		real(dp)	:: kxMin, kyMin
		!
		kxMin	= -1.0_dp * PI_dp /  aX
		dkx		=  2.0_dp * PI_dp / (aX * nKx)
		kyMin	= -1.0_dp * PI_dp /  aY
		dky		=  2.0_dp * PI_dp / (aY * nKy)
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
	end


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
	end 


	subroutine popGvec()
		!populates the G vector (basis vector)
		integer		:: i, ix, iy, nxMin, nyMin, nGdim
		real(dp)	:: thres
		!
		thres	= 1e-15_dp
		nGdim	= int (		sqrt(  real(nG,dp)	)		)
		nxMin	= -(nGdim-1)/2
		nyMin	= nxMin	
		!
		do iy = 1, nGdim
			do ix = 1, nGdim
				i = (iy-1) * nGdim + ix
				Gvec(1,i)	= (nxMin+ix-1) * 2*PI_dp / aX		!x component 
				Gvec(2,i)	= (nyMin+iy-1) * 2*PI_dp / aY		!y component
				if( 	abs( Gvec(1,i) ) < thres 	.and.	abs( Gvec(2,i) ) < thres ) then
					nG0 = i
				end if
			end do
		end do
		!
		return
	end


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
	end

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
	end


	subroutine calcRcell()
		integer		:: nI, nJ, n

		do nJ = 1, nSCy
			do nI = 1, nSCx
				n = (nJ-1) * nSCx + nI  !rI	=	(rIy-1) * nRy + rIx
				Rcell(1,n)	= (nI-1) * aX
				Rcell(2,n)	= (nJ-1) * aY
			end do
		end do

		return
	end




end module sysPara