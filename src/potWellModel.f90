module potWellModel
	!this modules sets up the Hamiltonian matrix and solves it for each k point
	!	in the process the wannier functions are generated aswell with routines from wannGen module
	use mathematics,	only:	dp, PI_dp,i_dp, myExp, myLeviCivita, eigSolver, nIntegrate, isHermitian
	use sysPara,		only: 	readInp, getKindex, getRindex, &
									dim, aX, aY,vol, nAt, atR, atPos, atPot,&
									nG, nG0, Gcut, nK, nKx, nKy, nWfs, nSC, nSCx, nSCy, nR, nRx, nRy, dx, dy, dkx, dky, &
									Gvec, atPos, atR, kpts, rpts, gaugeSwitch
	use wannGen,		only:	projectBwf, genWannF
	use output,			only:	printMat
	implicit none	
	
	private
	public ::					solveHam, calcVeloMat, calcConn, calcCurv



		

	


	contains
!public:
	subroutine solveHam(U, wnF, unkW, En, veloBwf)
		!solves Hamiltonian at each k point
		!also generates the Wannier functions on the fly (sum over all k)
		complex(dp),	intent(out)		::	U(:,:,:), wnF(:,:,:), unkW(:,:,:), veloBwf(:,:,:)		!Uh(  nWfs	, nWfs,		nK		)	
																				!wnF( nR	, nSupC,	nWfs	)	
																				!unkW(nR	, nKpts,	nWfs	)
																				!veloBwf(nR,nK,2*nG)
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
			En(ki,:) = EnT(1:nWfs) 
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
		!calculates matrix elements of the velocity operator
		!	velocity operator is analytically applied to the plane wave basis 
		!	and then weighted by the basCoeff obtained from the solver and stored in veloBwf
		!	
		complex(dp),		intent(in)		:: unk(:,:,:), veloBwf(:,:,:)	!	unk(nR, nK, nWfs) , veloBwf(nR,nK ,2*nWfs)
		complex(dp),		intent(out)		:: Velo(:,:,:,:)   !Velo(3,	nK,	nWfs	, nwFs)		
		integer								:: ki, n,m, ri
		complex(dp),		allocatable		:: fx(:), fy(:)

		allocate(	fx(nR)	)
		allocate(	fy(nR)	)

		
		do m = 1, nWfs
			do n = 1, nWfs
				do ki = 1, nK
					!FILL INTEGRATION ARRAY
					do ri = 1, nWfs
						fx(ri)	= myExP( 	dot_product( kpts(:,ki), rpts(:,ri) )		)	*unk(ri,ki,n)	* veloBwf(ri, ki,	m		)
						fy(ri)	= myExP( 	dot_product( kpts(:,ki), rpts(:,ri) )		)	*unk(ri,ki,n)	* veloBwf(ri, ki,	nWfs+m	)
					end do
					!INTEGRATE
					Velo(1,ki,n,m)	= nIntegrate(nR, nRx, nRy, dx, dy, fx)
					Velo(2,ki,n,m)	= nIntegrate(nR, nRx, nRy, dx, dy, fy)
					Velo(3,ki,n,m) 	= dcmplx(0.0_dp)
				end do
			end do
		end do
		!	
		!
		return
	end





	subroutine calcConn(unk,nxk, nyk, A)
		!finite difference on lattice periodic unk to calculate the Berry connection A
		!	A_n(k) 	= <u_n(k)|i \nabla_k|u_n(k)>
		!		 	= i  <u_n(k)| \sum_b{ w_b * b * [u_n(k+b)-u_n(k)]}
		!			= i \sum_b{		w_b * b * [  <u_n(k)|u_n(k+b)> -  <u_n(k)|u_n(k)>]		}
		!
		! see Mazari, Vanderbilt PRB.56.12847 (1997), Appendix B
		!
		complex(dp),	intent(in)		:: unk(:,:,:)		!unk(	nR, nK/nKw, nWfs/nG	)
		integer,		intent(in)		:: nxk, nyk
		complex(dp),	intent(out)		:: A(:,:,:)			!Aconn(	3,nK, nWfs)		)	
		complex(dp)						:: Mxl, Mxr, Myl, Myr, M, one
		integer							:: n, Z, ki, kx, ky, kxl, kxr, kyl, kyr
		real(dp)						:: thres, wbx,wby, bxl(2), bxr(2), byl(2), byr(2) !for nearest neighbours, assuming cubic mesh
		!
		thres	= 1e-3_dp
		A 		= dcmplx(0.0_dp)
		Z 		= 4	!amount of nearest neighbours( 2 for 2D cubic unit cell)
		wbx 	= 3.0_dp / 		( real(Z,dp) * dkx**2 )
		wby 	= 3.0_dp /		( real(Z,dp) * dky**2 )
		!b vector two nearest X neighbours:
		bxl(1) 	= -dkx				
		bxl(2)	= 0.0_dp
		bxr(1) 	= +dkx
		bxr(2)	= 0.0_dp
		!b vector two nearest Y neighbours:
		byl(1) 	= 0.0_dp
		byl(2)	= -dky
		byr(1) 	= 0.0_dp
		byr(2)	= +dky


		do n = 1, nWfs
			do kx = 1, nxk
				kxl	= getLeft(kx,nxk)
				kxr	= getRight(kx,nxk)
				!
				do ky = 1, nyk
					kyl	= getLeft(ky,nyk)
					kyr = getRight(ky,nyk)
					ki	= getKindex(kx,ky)
					!
					!OVERLAP TO NEAREST NEIGHBOURS
					one = UNKoverlap(	n, 		ki		, 		ki					, unk	)
					Mxl	= UNKoverlap(	n, 		ki		, getKindex( kxl, ky ) 		, unk	) 
					Mxr	= UNKoverlap(	n, 		ki		, getKindex( kxr, ky )		, unk	)
					Myl	= UNKoverlap(	n, 		ki		, getKindex( kx ,kyl )		, unk	)
					Myr	= UNKoverlap(	n, 		ki		, getKindex( kx ,kyr )		, unk	)
					!
					!write(*,'(a,f15.12,a,f15.12)')"[calcConn]: Mxl=",dreal(Mxl),"+i*",dimag(Mxl)
					!FD SUM OVER NEAREST NEIGHBOURS
					A(1:2,ki,n) = A(1:2,ki,n) + wbx * bxl(:) * i_dp * ( Mxl - one )
					A(1:2,ki,n) = A(1:2,ki,n) + wbx * bxr(:) * i_dp * ( Mxr - one )
					A(1:2,ki,n) = A(1:2,ki,n) + wby * byl(:) * i_dp * ( Myl - one )
					A(1:2,ki,n) = A(1:2,ki,n) + wby * byr(:) * i_dp * ( Myr - one )
					!FD SUM OVER NEAREST NEIGHBOURS
					!A(:,ki,n) = A(:,ki,n) + wbx * bxl(:) * dimag(	log( Mxl ) )
					!A(:,ki,n) = A(:,ki,n) + wbx * bxr(:) * dimag(	log( Mxr ) )
					!A(:,ki,n) = A(:,ki,n) + wby * byl(:) * dimag(	log( Myl ) )
					!A(:,ki,n) = A(:,ki,n) + wby * byr(:) * dimag(	log( Myr ) )
					!
					!
					if(abs( abs(one) - 1.0_dp ) > thres ) then
						write(*,'(a,i2,a,i7,a,f16.8,a,f16.8)')	"[calcConn]: n=",n," unk normalization problem at ki=",ki,&
													" one=",dreal(one),"+i*",dimag(one)
					end if
				end do
			end do
		end do
		!
		!
		return
	end




	subroutine calcCurv(En, Velo, Fcurv)
	!Calculates the connection via Kubo formula on matrix elements of velocity operator
	!see Wang/Vanderbilt PRB 74, 195118 (2006) eq.(5)
	!
	!	F_n,c(k) = \sum{a,b} leviCivi(a,b,c) F_n,{a,b}(k)
	!
	!
		real(dp),		intent(in)		:: En(:,:)		!En(nK,nWfs)
		complex(dp),	intent(in)		:: Velo(:,:,:,:) !Velo(3,nK,nWfs,nWfs)
		complex(dp),	intent(out)		:: Fcurv(:,:,:)  !Fcurv(3,nK,nWfs)
		integer							:: n, ki, a,b,c
		!
		Fcurv = dcmplx(0.0_dp)
		do n = 1, nWfs
			do ki = 1, nK
				do c = 1,3
					do b= 1,3
						do a=1,3
							if( myLeviCivita(a,b,c) /= 0) then
								Fcurv(c,ki,n) = Fcurv(c,ki,n) + myLeviCivita(a,b,c) * omega(n,a,b,ki,Velo, En)
							end if
						end do 
					end do
				end do
			end do
		end do
		!
		!
		return
	end



	real(dp) function omega(n,a,b,ki,Velo,En)
		!returns curvature tensor element a,b
		!see Wang/Vanderbilt PRB 74, 195118 (2006) eq.(9)
		!
		integer,		intent(in)		:: n, a, b, ki
		complex(dp),	intent(in)		:: Velo(:,:,:,:)!Velo(3,nK,nWfs,nWfs)
		real(dp),		intent(in)		:: En(:,:)!En(nK,nWfs)
		integer							:: m
		!
		omega	= 0.0_dp
		do m = 1, nWfs
			if(m /= n) then
				omega = omega +		dimag(	Velo(a,ki,n,m) * Velo(b,ki,m,n)	) 	/ 	( En(ki,m) - En(ki,n) )**2
			end if 
		end do		
		omega	= -2.0_dp * omega
		!
		!
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




	!GENERATING BLOCH WAFEFUNCTIONS & LATT PERIOD FUNCTIONS
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




	!VELOCTIY MATRIX
	subroutine calcVeloBwf(ki, basCoeff, veloBwf)
		!calculates the bwf with changed basis vectors.
		!	basis vetors have the \nabla operator applied to the plane waves included
		!
		!	|veloBwf> = \hat{v} |bwf>
		!
		integer,		intent(in)		:: ki
		complex(dp),	intent(in)		:: basCoeff(:,:)
		complex(dp),	intent(out)		:: veloBwf(:,:,:)	!veloBwf(nR,nK,2*nG)
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
			veloBwf(xi,ki,1:nWfs)	= tmp(1:nWfs)
			!Y COMPONENT
			min						= nG+1
			max						= 2*nG
			tmp						= matmul(	basVec(min:max),	basCoeff) / dsqrt(vol)
			min						= nWfs +1
			max						= 2*nWfs
			veloBwf(xi,ki,min:max)	= tmp(1:nWfs)
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




	!CONNECTION HELPERS
	complex(dp) function UNKoverlap(n, ki, knb, unk)
		!HELPER for calcConn
		!calculates the overlap between unk at ki and at a neigbhouring k point knb
		!	integration only over the first unit cell
		!
		integer,		intent(in)		:: n, ki, knb
		complex(dp),	intent(in)		:: unk(:,:,:)  !unk(	nR, nK, nWfs/nG	)
		complex(dp),	allocatable		:: f(:)
		integer							:: xi,yi,ri,rloc, nRx1, nRy1, nR1
		!
		!Set integration range to first unit cell
		nRx1 	= int(		real(nRx,dp) / real(nSCx,dp)		)
		nRy1 	= int(		real(nRy,dp) / real(nSCy,dp)		)
		nR1		= nRx1 * nRy1 
		allocate(	f(nR1)	)
		!
		!fill integration array
		f 		= dcmplx(0.0_dp)
		do yi = 1, nRy1
			do xi = 1, nRx1
				ri		= getRindex(xi,yi)			!overall index, to get correct position from unk
				rloc 	= (yi-1) * nRx1 + xi		!for mapping to f array
				f(rloc)	= dconjg( unk(ri,ki,n) ) * unk(ri,knb,n)
				!write(*,'(a,f10.6,a,f10.6)')	"[overlap] f=",dreal(f(rloc)),"+i*",dimag(f(rloc))
			end do
		end do
		!
		!integrate
		UNKoverlap = nIntegrate(nR1, nRx1, nRy1, dx, dy, f	)
		
		!write(*,'(a,f10.6,a,f10.6)')"[overlap]=",dreal(overlap),"+i*",dimag(overlap)
		!
		!
		return
	end


	integer function getLeft(i,N)
		!HELPER for calcConn
		!gets left (lower) neighbour, using the periodicity at boundary
		!
		integer,	intent(in)	:: i,N
		if(i==1) then
			getLeft = N
		else
			getLeft = i-1
		end if
		!
		return
	end


	integer function getRight(i,N)
		!HELPER for calcConn
		!gets right (upper) neighbour, using the periodicity at boundary
		!
		integer,	intent(in)	:: i,N
		if(i==N) then
			getRight = 1
		else
			getRight = i+1
		end if
		!
		return
	end











end module potWellModel