module potWellModel
	!this modules sets up the Hamiltonian matrix and solves it for each k point
	!	in the process the wannier functions are generated aswell with routines from wannGen module
	use mathematics,	only:	dp, PI_dp,i_dp, myExp, eigSolver, nIntegrate, isHermitian
	use sysPara,		only: 	readInp, getKindex, &
									dim, aX, aY,vol, nAt, atR, atPos, atPot,&
									nG, nG0, Gcut, nK, nKx, nKy, nWfs, nSC, nR, nRx, nRy, dx, dy, dkx, dky, &
									Gvec, atPos, atR, kpts, rpts, gaugeSwitch
	use wannGen,		only:	projectBwf, genWannF
	use output,			only:	printMat
	implicit none	
	
	private
	public ::					solveHam, calcVeloMat



		

	


	contains
!public:
	subroutine solveHam(U, wnF, unkW, En, veloBwf)
		!solves Hamiltonian at each k point
		!also generates the Wannier functions on the fly (sum over all k)
		complex(dp),	intent(out)		::	U(:,:,:), wnF(:,:,:), unkW(:,:,:), veloBwf(:,:,:)		!Uh(  nWfs	, nWfs,		nK		)	
																				!wnF( nR	, nSupC,	nWfs	)	
																				!unkW(nR	, nKpts,	nWfs	)
																				!veloBwf(nR,2*nG,nK)
		real(dp),		intent(out)		::	En(:,:)																	
		complex(dp),	allocatable		::	Hmat(:,:), bWf(:,:), lobWf(:,:), gnr(:,:)
		real(dp),		allocatable		::	EnT(:), bwfR(:,:), bwfI(:,:)	 
		integer							:: 	ki, xi , n
		real(dp)						::	kVal(dim)
		!
		allocate(	Hmat(	nG,	nG	)			)
		allocate(	EnT(		nG		)			)
		allocate(	bWf(	nR, nG	)			)
		allocate(	bWfR(	nR, nG	)			)
		allocate(	bWfI(	nR, nG	)			)
		allocate(	lobWf(	nR, nWfs)			)
		allocate(	gnr(	nR,	nWfs)			)
		
		wnF		=	dcmplx(0.0_dp)
		bWf		=	dcmplx(0.0_dp)
		unkW	=	dcmplx(0.0_dp)
		veloBwf	=	dcmplx(0.0_dp)
		!
		
		!call genTrialOrb(gnr) !ToDo
		!
		open(unit=200, file='rawData/bandStruct.dat', form='unformatted', access='stream', action='write')
		open(unit=210, file='rawData/bwfR.dat'		, form='unformatted', access='stream', action='write')
		open(unit=211, file='rawData/bwfI.dat'		, form='unformatted', access='stream', action='write')
		do ki = 1, nK
			write(*,*)"[solveHam]: ki=",ki
			kVal	=	kpts(:,ki)
			!
			call populateH(kVal, Hmat)	
			call eigSolver(Hmat, EnT)
			write(200)	EnT
			En(:,ki) = EnT(1:nWfs) 
			!
			call gaugeCoeff(kVal, Hmat)
			call genBlochWf(ki, Hmat, bWf)		
			call calcVeloBwf(ki,Hmat, veloBwf)
			bWfR 	= dreal(bWf)                 !Todo fix that
			bWfI	= dimag(bWf)
			write(210)bWfR
			write(211)bwfI
			!
			call projectBwf(ki, bWf, loBwf, U(:,:,ki))
			
			call genWannF(ki, lobWf, wnF)

			call genUnk(ki, bWf, unkW(:, ki, :))
			!
		end do
		close(200)
		close(210)
		close(211)
		!
		return
	end




	subroutine calcVeloMat(unk, veloBwf, Velo)
		complex(dp),		intent(in)		:: unk(:,:,:), veloBwf(:,:,:)	!	unk(nR, nK, nWfs) , veloBwf(nR,2*nWfs,nK)
		complex(dp),		intent(out)		:: Velo(:,:,:,:)   !Velo(		3,			nWfs	, nwFs,	nK)		
		integer								:: ki, n,m, ri
		complex(dp),		allocatable		:: fx(:), fy(:)

		allocate(	fx(nR)	)
		allocate(	fy(nR)	)

		do ki = 1, nK
			do m = 1, nWfs
				do n = 1, nWfs
					!FILL INTEGRATION ARRAY
					do ri = 1, nWfs
						fx(ri)	= myExP( 	dot_product( kpts(:,ki), rpts(:,ri) )		)	*unk(ri,ki,n)	* veloBwf(ri,	m		,ki)
						fy(ri)	= myExP( 	dot_product( kpts(:,ki), rpts(:,ri) )		)	*unk(ri,ki,n)	* veloBwf(ri,	nWfs+m	,ki)
					end do
					!INTEGRATE
					Velo(1,n,m,ki)	= nIntegrate(nR, nRx, nRy, dx, dy, fx)
					Velo(2,n,m,ki)	= nIntegrate(nR, nRx, nRy, dx, dy, fy)
					Velo(3,n,m,ki) 	= dcmplx(0.0_dp)
				end do
			end do
		end do


		return
	end












!private:
	!POPULATION OF H MATRIX
	subroutine populateH(k, Hmat)
		!populates the Hamiltonian matrix by adding 
		!	1. kinetic energy terms
		!	2. potential terms
		!and checks if the resulting matrix is still hermitian( for debugging)
		real(dp)   , intent(in)    :: k(dim)
		complex(dp), intent(inout) :: Hmat(:,:)
		!init to zero
		Hmat = dcmplx(0.0_dp) 
		!
		!FEATURES:
		call Hkin(k,Hmat) 
		call Hpot( Hmat)
		!
		!DEBUGGING:
		if ( .not.	isHermitian(Hmat)	) then
			write(*,*)"[populateH]: Hamiltonian matrix is not Hermitian :"
			call printMat(nG, Hmat)
		else
			write(*,*)"[populateH]: Hmat is hermitian"
		end if
		return
	end


	subroutine Hkin(k, Hmat)
		!add kinetic energy to Hamiltonian
		real(dp)   , intent(in)	   :: k(2)
		complex(dp), intent(inout) :: Hmat(:,:)
		integer 				   :: i
		real(dp)				   :: fact, kg(2)
		fact = 0.5_dp	!hbar*hbar/(2*me) in a.u.
		do i = 1, nG
			kg(:) = k(:) + Gvec(:,i)
			Hmat(i,i) = Hmat(i,i) +  dcmplx(	fact * dot_product(kg,kg) , 0.0_dp	) 
		end do
		return
	end


	subroutine Hpot( Hmat)
		!Add potential to Hamiltonian
		complex(dp) , intent(inout) :: Hmat(:,:)
		integer 					:: i,j
		complex(dp)					:: tmp
		do i = 1,nG
			do j = 1,nG
				tmp			= V(i,j)
				!write(*,'(a,i2,a,i2,a,f15.10,a,f15.10)')	"[Hpot]: V(",i,",",j,") =",dreal(tmp),"+i*",dimag(tmp)
				Hmat(i,j)	= Hmat(i,j) + tmp
			end do
		end do
		!
		return
	end


	complex(dp) function V(i,j)
		!calc potential matrix elements
		!the integrals were solved analytical and are hard coded in this function
		integer,	intent(in)	::	i, j
		integer					::	at
		complex(dp)				::	Vx, Vy, Vpot
		real(dp)				::  xL, yL, xR, yR, dGx, dGy, thres
		!
		V 		= dcmplx(0.0_dp)
		dGx		= Gvec(1,j) - Gvec(1,i)
		dGy		= Gvec(2,j) - Gvec(2,i)
		thres	= 1e-15_dp
		!
		do at = 1, nAt
			Vpot	=	dcmplx(atPot(at))
			xL	=	atPos(1,at) - atR(1,at)
			xR	=	atPos(1,at) + atR(1,at) 
			yL	=	atPos(2,at) - atR(2,at)
			yR	=	atPos(2,at) + atR(2,at) 
			!
			if( abs(dGx) < thres ) then
				!write(*,*)"zero x difference i=",i," j=",j
				Vx	= ( xR - xL )
			else
				Vx	= -i_dp	*	(	myExp(dGx*xR) - myExp(dGx*xL)	)	/ dcmplx( dGx )
			end if
			

			if( abs(dGy) < thres ) then
				!write(*,*)"zero y difference i=",i," j=",j
				Vy	= ( yR - yL )
			else
				Vy	= -i_dp *	(	myExp(dGy*yR) - myExp(dGy*yL)	)	/ dcmplx( dGy )
			end if
			!
			!
			V	= V +	Vpot*Vx*Vy
			!write(*,'(a,i2,a,i2,a,f15.10)')	"i=",i," j=",j,": atPot=",Vpot
			!write(*,'(a,i2,a,i2,a,f15.10,a,f15.10)')	"i=",i," j=",j,": V=",dreal(V),"+i*",dimag(V)
		end do
		!
		return
	end



	!GAUGING BASIS COEFFICIENTS
	subroutine gaugeCoeff(k, basCoeff)
		!method gauges the basis coefficients directly
		!	controled via the gaugeSwitch from the input file
		!	this should only be used for testing & understanding purposes,
		!	since this gauge trafo is not saved in the U matrix (U matrix rotates between ham & wann gauge)
		real(dp),		intent(in)		:: k(dim)
		complex(dp),	intent(inout) 	:: basCoeff(:,:)
		!
		select case (gaugeSwitch)
			case (1)
				write(*,*) "[gaugeCoeff]: gauging with G=0 phase"
				call gaugeCoeff1(basCoeff)
			!case (2)
			!	write(*,*) "[gaugeCoeff]: gauging with G(nG0-1) for k<=0, and G(nG0+1) for k>0. (G(nG0)=0)"
			!	call gaugeCoeff2(k,basCoeff)
			case (0)
				!write(*,*) "[gaugeCoeff]: use gauge convention from solver"
				! do nothing
			case default
				write(*,*) "[gaugeCoeff]: unknown gaugeSwitch, eigVecs wont be gauged"
		end select
		!
		return
	end


	subroutine gaugeCoeff1(basCoeff)
		!Gauges basCoeff with the G=0 component, effectively removes k-dependence of basCoeff gauge
		complex(dp)	, intent(inout) :: basCoeff(:,:)
		integer						:: i1,i2
		complex(dp)				  	:: phase0
		real(dp)				   	:: delta
		!
		delta=1E-14_dp
		!
		do i2=1,nG !loop eigVecs
			if(	abs( basCoeff(nG0,i2) )	> delta	) then

				phase0 = basCoeff(nG0,i2) / abs(basCoeff(nG0,i2)) 
				do i1=1,nG !gauge each vector
						basCoeff(i1,i2) = basCoeff(i1,i2)/ phase0
				end do
			endif
		end do
		!
		return
	end


	!subroutine gaugeCoeff2(k, basCoeff)
	!	!what does k < 0 mean here ?!, makes only sense in 1D
	!	!for negative k gauge with Gmin, for positive with Gmax
	!	real(dp),		intent(in)		:: k(dim)
	!	complex(dp), 	intent(inout) 	:: basCoeff(:,:)
	!	integer 				   		:: i1,i2
	!	complex(dp)				   		:: phase0
	!	!
	!	do i2=1,nG
	!		if( k <= 0.0_dp) then
	!			phase0 = basCoeff(nG0-1 ,i2) / abs(basCoeff(nG0-1 ,i2))
	!		else
	!			phase0 = basCoeff(nG0+1,i2) / abs(basCoeff(nG0+1,i2))
	!		end if
	!		do i1=1,nG
	!			basCoeff(i1,i2) = basCoeff(i1,i2)/ phase0
	!		end do
	!	end do
	!	!
	!	return
	!end



	!GENERATING BLOCH WAFEFUNCTIONS
	subroutine genBlochWf(ki,basCoeff, bWf)
		!generates the bloch wavefunctions, with  the basCoeff from eigSolver
		integer		, intent(in)	:: ki
		complex(dp)	, intent(in)	:: basCoeff(:,:)
		complex(dp)	, intent(out)	:: bWf(:,:)	!bWf(nRpts,nG)			
		complex(dp)	, allocatable	:: basVec(:)
		integer 				 	:: xi,n
		allocate(	basVec(nG)	)
		!
		do xi = 1, nR
				call calcBasVec(ki,xi, basVec)
				bWf(xi,:) = matmul(	 basVec , basCoeff	)  /  dsqrt(vol)
		end do
		!
		return 
	end


	subroutine calcBasVec(ki, ri, basVec)
		!calculates the basis vectors e^i(k+G).r
		!	if |k+G| is larger then the cutoff the basis vector is set to zero
		!	the cutoff enforces a symmetric base at each k point
		integer,	 intent(in)  :: ki, ri
		complex(dp), intent(out) :: basVec(:)
		real(dp)				 :: tmp(dim)
		integer 				 ::	i 
		!

		do i =1, nG
			tmp(:) = kpts(:,ki) + Gvec(:,i)
			!
			if( norm2(tmp) < Gcut ) then
				basVec(i) = myExp( 		dot_product( tmp, rpts(:,ri) )			)
			else
				basVec(i) = dcmplx( 0.0_dp )
			end if
		end do
		!
		return
	end



	subroutine calcVeloBwf(ki, basCoeff, veloBwf)
		!calculates the bwf with changed basis vectors.
		!	basis vetors have the \nabla operator applied to the plane waves included
		!
		!	|veloBwf> = \hat{v} |bwf>
		!
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: basCoeff(:,:)
		complex(dp),	intent(out)		:: veloBwf(:,:,:)	!veloBwf(nR,2*nG,nK)
		complex(dp),	allocatable		:: basVec(:), tmp(:)
		integer							:: xi, min, max

		allocate(	basVec(2*nG)	)
		allocate(	tmp(nG)		)


		do xi = 1, nR
			call calcVeloBasVec(ki,xi,basVec)
			!X COMPONENT
			min 					= 1
			max 					= nG
			tmp						= matmul(	basVec(min:max),	basCoeff) / dsqrt(vol)
			veloBwf(xi,1:nWfs,ki)	= tmp(1:nWfs)
			!Y COMPONENT
			min						= nG+1
			max						= 2*nG
			tmp						= matmul(	basVec(min:max),	basCoeff) / dsqrt(vol)
			min						= nWfs +1
			max						= 2*nWfs
			veloBwf(xi,min:max,ki)	= tmp(1:nWfs)
		end do
		!
		return
	end


	subroutine calcVeloBasVec(ki,ri,basVec)
		integer,		intent(in)		:: ki, ri
		complex(dp),	intent(out)		:: basVec(:)
		real(dp)				 :: tmp(2)
		integer 				 ::	i 

		do i =1, nG
			tmp(:) = kpts(:,ki) + Gvec(:,i)
			!
			if( norm2(tmp) < Gcut ) then
				!X COMPONENT
				basVec(i) 		= i_dp * (	kpts(1,ki) + Gvec(1,i)	) * myExp( 		dot_product( tmp, rpts(:,ri) )			)
				!Y COMPONENT
				basVec(i+nG)	= i_dp * (	kpts(2,ki) + Gvec(2,i)	) * myExp(		dot_product( tmp, rpts(:,ri) )			)
			else
				basVec(i) = dcmplx( 0.0_dp )
			end if
		end do


		return
	end


























	subroutine genUnk(ki, bWf, unk)
		! generates the lattice periodic part from given bloch wave functions
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: bWf(:,:) !lobWf(	nR, nWfs)
		complex(dp),	intent(out)		:: unk(:,:)   !unk(	nR, nWfs)
		integer							:: xi, n
		complex(dp)						:: phase
		!
		do n = 1, nWfs
			do xi = 1, nR
				phase = myExp( -1.0_dp 	*	 dot_product( kpts(:,ki) , rpts(:,xi)	) 			)
				unk(xi,n) = phase * Bwf(xi,n)
			end do
		end do
		!
		return
	end


	logical function isLattSym(bWf)
		!checks if bwf(k) = bwf(k+G)
		complex(dp),	intent(in)		:: bWf(:,:,:) !nR, nK , nG or nWfs
		integer							:: k00, k10, k01, k11, n ! edge point indices

		isLattSym = .true.
		k00 = getKindex(	1	, 1		)
		k10	= getKindex(	nKx	, 1		)
		k01 = getKindex(	1	, nKy	)
		k11 = getKindex(	nKx , nKy	)
		write(*,'(a,i3,a,i3,a,i3,a,i3)')"[isLattSym]: k00 =",k00,", k10=",k10,", k01=",k01,", k11=",k11 


		do n = 1, size(bwf,3) ! loop states

		end do

		return
	end








end module potWellModel