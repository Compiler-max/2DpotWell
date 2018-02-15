module mathematics
	!module contains some basic constants, especially dp for double precision convention
	!additionally it contains the matrix operation wrappers to call intel mkl
	implicit none
	private


	public :: 	dp, PI_dp, i_dp, acc, aUtoAngstrm, aUtoEv, aUtoTesla, setAcc, machineP, & 
				myExp, Cangle, myLeviCivita,  crossP,& 
				isUnit, isIdentity, isHermitian, mySVD,&
				eigSolverFULL, eigSolverPART, myMatInvSqrt, rotMat, myCommutat


	interface nIntegrate
		module procedure nIntegrateREAL
		module procedure nIntegrateCPLX
	end interface nIntegrate

	interface crossP
		module procedure crossPreal
		module procedure crossPcplx
	end interface crossP


	!for clean double precision convention through the code
	integer, 		parameter 	:: 	dp 				= kind(0.d0)
	real(dp), 		parameter	:: 	machineP 		= 1e-15_dp
	real(dp)					:: 	acc				= 1e-14_dp
	
	!mathematical constants
	real(dp), 		parameter 	::	PI_dp 			= 4 * datan (1.0_dp)
	complex(dp),	parameter 	::	i_dp 			= dcmplx(0.0_dp, 1.0_dp)



	!physical constants
	real(dp),		parameter 	::	aUtoAngstrm 	= 0.52917721092_dp
	real(dp),		parameter 	::	aUtoEv	 		= 27.211385_dp
	real(dp),		parameter	::	aUtoTesla		= 235051.76_dp

	contains

















!public
	subroutine myFT_R_to_kdense()

		return
	end subroutine

	subroutine myFT_R_to_kcoarse()

		return
	end subroutine

	subroutine myFT_kcoarse_to_R()

		return
	end subroutine

	subroutine myFT_kdense_to_R()

		return
	end subroutine


	subroutine setAcc(thres)
		real(dp),		intent(in)	:: thres
		acc	= thres
		return
	end subroutine


	complex(dp) function myExp(x)
		!supposed to boost performance
		real(dp), intent(in) :: x
		!
		myExp = dcmplx(  dcos(x) , dsin(x) ) 
		!
		return
	end function


	real(dp) function Cangle(z)
		!returns principale value between -pi and pi
		complex(dp),	intent(in)	:: z
		!
		Cangle = datan2(dimag(z),dreal(z)) 
		if(Cangle< 0.0_dp) then
			Cangle = 2*PI_dp + Cangle 
		end if
		Cangle = Cangle / PI_dp
		!
		return
	end function


	integer function myLeviCivita(i,j,k)
		!Hard coded Levi Civita tensor
		integer,		intent(in)		:: i,j,k
		logical							:: even, odd
		!
		!
		even	= (i==1 .and. j==2 .and. k==3) .or. (i==2 .and. j==3 .and. k==1) .or. (i==3 .and. j==1 .and. k==2)
		odd		= (i==3 .and. j==2 .and. k==1) .or. (i==1 .and. j==3 .and. k==2) .or. (i==2 .and. j==1 .and. k==3)
		!
		if(even) then
			myLeviCivita	=  1
		else if (odd) then
			myLeviCivita	= -1
		else
			myLeviCivita	=  0
		end if
		!
		!DEBUGGING
		if(even .and. odd) then
			write(*,*)"[myLeviCivita]: myLeviCivita detected even and odd, contact the programer he fucked up"
		end if
		!
		return
	end function


	subroutine mySVD(Mat, U,s,Vt)
		!single value decomposition performed with zgesvd:
		!https://software.intel.com/en-us/node/469236
		complex(dp),	intent(in)		:: Mat(:,:)
		complex(dp),	intent(out) 	:: U(:,:), Vt(:,:)
		real(dp),		intent(out) 	:: s(:)
		character*1						:: jobu, jobvt
		integer							:: m, n, lda, sda, ldu, ldvt, lwork, rworkS, info
		complex(dp),	allocatable 	:: A(:,:), work(:), rwork(:)
		!
		m 		= size(	Mat	,	1						)
		n 		= size(	Mat	,	2						)
		lda 	= max(	1	,	m						)
		sda		= max(	1	,	n						)
		lwork	= max(	1	,	2*min(m, n)+max(m, n) 	)
		rworkS	= max(	1	,	5*min(m, n)				)
		ldu		= m
		ldvt	= n
		allocate(	A(  lda ,sda  )				)
		allocate(	work(  lwork ) 				)
		allocate(	rwork(  rworkS  )			)
		!
		jobu 	= 'S' !calculate all rows/columns in U matrix
		jobvt	= 'S' !calculate all rows/columns in V matrix
		A = Mat
		!
		U = dcmplx(0.0_dp)
		Vt= dcmplx(0.0_dp)

		!all zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
		call zgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, Vt, ldvt, work, lwork, rwork, info)

		!error handeling
		if (info .ne. 0) then
			if(info < 0) then
				write(*,*)"error during single value decomposition: the ",-info,"th component had and illegal value"
			else if(info > 0) then
				write(*,*)"error during single value decomposition: the ",info,"th superdiagonal did not converge"
			else
				write(*,*)"unknown error during single value decomposition, contact somebody for help. Ideally not the dev...."
			end if
		end if
		!
		return
	end subroutine


	subroutine eigSolverFULL(A, w)
		!solves standard (non general) eigenvalue problem with following mkl routine 
		!return eigVectors stored in A, eigValues in w
		!https://software.intel.com/en-us/node/469182        
		complex(dp) , intent(inout)   :: A(:,:)
		real(dp)    , intent(out)     :: w(:)
            

            !set up workspace arrays (assuming jobz='V' and using zheevd)
            complex(dp), allocatable, dimension(:)   :: work
            real(dp)   , allocatable, dimension(:)   :: rwork
            integer    , allocatable, dimension(:)   :: iwork
            character*1                              :: jobz,uplo
            integer                                  :: n, info,lwork,lrwork,liwork
            n		= size(A,1)
            if(n /= size(A,2)) then
            	write(*,*)"[eigSolver]: WARNING the matrix to solve is not a square matrix"
            end if
            lwork  	=   n*n + 2*n
            lrwork 	= 2*n*n + 5*n + 1
            liwork 	=         5*n + 3
            allocate(  work( lwork) )
            allocate( rwork(lrwork) )
            allocate( iwork(liwork) )

            jobz='V'
            uplo='U' !test this behaviour

            !solve eigenvalue problem
            !call zheevd(jobz, uplo, n, a,lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
             call zheevd(jobz, uplo, n, A, n , w, work, lwork, rwork, lrwork, iwork, liwork, info)
            
            !check if system was solved correctly
            call errCheck(n,info,jobz)
            !
            if(.true.) then
                  call normalCheck(n,A)
            end if
            !
            return
	end subroutine


	subroutine eigSolverFULL_NEW(A, w, z, m)
		!https://software.intel.com/en-us/mkl-developer-reference-fortran-heevr#6ADF761A-127A-4C9B-9A2A-1A8AA4602CE1
		!uses same routine as eigSolver2
		complex(dp),	intent(inout)			:: 	A(:,:)
		real(dp),		intent(out)				::	w(:)
		complex(dp),	intent(out)				::	z(:,:)
    	integer,		intent(out)				::	m
    	character*1	 							::	jobz, range, uplo
    	integer									::	n, lda, il, iu, ldz, lwork, lrwork, liwork,  info 
    	integer,		allocatable				::	isuppz(:) , iwork(:)
    	real(dp)								::	vl, vu, abstol
    	real(dp),		allocatable				::	rwork(:)
    	complex(dp),	allocatable				::	work(:)		

    	jobz	= 'V'
    	range	= 'A'
    	uplo	= 'U'
    	n 		= size(a,1)
    	lda		= n
    	vl		= 0.0_dp
    	vu  	= 0.0_dp
    	il		= 0
    	iu		= 0
    	abstol	= 1e-15_dp
    	ldz		= n
    	lwork	= 2*n
    	lrwork	= 24*n
    	liwork	= 10*n
    	if( size(z,1)/= ldz ) write(*,*)"[eigSolver2]: z array has wrong size"
    	if( size(w)/= n) write(*,*)"[eigSolver2]; w array has wrong size"
    	!
    	allocate( isuppz(2*iu)	)
    	allocate(  work(lwork)	)
    	allocate( rwork(lrwork)	)
    	allocate( iwork(liwork)	)	
    	
    	!
    	call zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, & 
    					m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
    	!
    	!
    	if( info < 0 ) then
    		write(*,*)	"[eigSolver2]: the ",-info,"th parameter had an illegal value"
    	else if( info > 0 ) then
    		write(*,*)	"[eigSolver2]: internal error has occured"
    	end if
    	!
    	!
		return
	end subroutine




	subroutine eigSolverPART(a, w ,z, m)
    	!https://software.intel.com/en-us/mkl-developer-reference-fortran-heevr#6ADF761A-127A-4C9B-9A2A-1A8AA4602CE1
    	complex(dp),	intent(inout)			:: 	a(:,:)
    	real(dp),		intent(out)				:: 	w(:)
    	complex(dp),	intent(out)				:: 	z(:,:)
    	integer,		intent(out)				:: 	m
    	character*1	 							:: 	jobz, range, uplo
    	integer									:: 	n, lda, il, iu, ldz, lwork, lrwork, liwork,  info 
    	integer,		allocatable				:: 	isuppz(:) , iwork(:)
    	real(dp)								:: 	vl, vu, abstol
    	real(dp),		allocatable				:: 	rwork(:)
    	complex(dp),	allocatable				:: 	work(:)
    	!
    	jobz	= 'V'
    	range	= 'I'
    	uplo	= 'U'
    	n 		= size(a,1)
    	lda		= n
    	vl		= 0.0_dp
    	vu  	= 0.0_dp
    	il		= 1
    	iu		= size(z,2)
    	abstol	= 1e-15_dp
    	ldz		= n
    	lwork	= 2*n
    	lrwork	= 24*n
    	liwork	= 10*n
    	if( size(z,1)/= ldz ) write(*,*)"[eigSolver2]: z array has wrong size"
    	if( size(w)/= n) write(*,*)"[eigSolver2]; w array has wrong size"
    	!
    	allocate( isuppz(2*iu)	)
    	allocate(  work(lwork)	)
    	allocate( rwork(lrwork)	)
    	allocate( iwork(liwork)	)	
    	
    	!
    	call zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &
    					m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
    	!
    	!
    	if( info < 0 ) then
    		write(*,*)	"[eigSolver2]: the ",-info,"th parameter had an illegal value"
    	else if( info > 0 ) then
    		write(*,*)	"[eigSolver2]: internal error has occured"
    	end if
    	!
    	!
    	return
    end subroutine


	subroutine myMatInvSqrt(Mat, minS, maxS)
		!routine calculates the inverse sqrt (Mat)^{-1/2} of the complex matrix Mat
		!
		complex(dp),	intent(inout)	:: Mat(:,:)
		real(dp),		intent(out)		:: minS, maxS
		complex(dp),	allocatable		:: U(:,:), Vt(:,:), TMP(:,:)
		real(dp),		allocatable		:: s(:)
		integer							:: i,j, m, n, k, lda, ldb, ldc
		complex(dp)						:: alpha, beta
		character*1						:: transa, transb
		
		m 		= size(	Mat	,	1 )
		n 		= size(	Mat	,	2 )	
		if(m /= n) then
			write(*,*)"[myMatInvSqrt]: WARNING not square matrix"
		end if
		allocate(	 U(  m , max(1,m)   )		)
		allocate(	Vt(  n , max(1,n)   )		)
		allocate(	s(  max(1,min(m,n))  )		)
		allocate(	TMP(  size(s), size(Vt,2)  ))
		!

		call mySVD(Mat, U,s,Vt)

		!diagnostics both values should be positive and larger then zero
		maxS	= maxval(s)
		minS	= minval(s)
		if( minS < 1e-15_dp) then
			write(*,'(a,e16.9,a,e16.9)')	"[myMatInvSqrt]: WARNING, minimum eigenvalue close to zero, minS=",&
												minS," maxS=",maxS
		else
			write(*,'(a,e16.9,a,e16.9)')	"[myMatInvSqrt]: minS=",minS," maxS=",maxS
		end if
		!
		!scaling, i.e. perform inversion and sqrt of eignevalues s(:)
		do i = 1, size(s)
				s(i) = 1.0_dp /	( dsqrt( s(i) ) + machineP )
		end do
		!
		!after scaling do the rotations
		do j = 1, size(Vt,2)
			do i = 1, size(s)
				TMP(i,j) = dcmplx( s(i) ) * Vt(i,j)
			end do
		end do
		

		!Mat = matmul(U, TMP)
		transa	= 'n'
		transb	= 'n'
		m		= size(U,1)
		n 		= size(TMP,2)
		k		= size(U,2)
		lda		= m
		ldb		= k
		ldc		= m
		alpha	= dcmplx(1.0_dp)
		beta	= dcmplx(0.0_dp)
		call zgemm(transa, transb, m, n, k, alpha, U, lda, TMP, ldb, beta, Mat, ldc)
		!
		!
		return
	end subroutine


	logical function isUnit(U)
		!return true if matrix u is unitary matrix
		!	i.e. I = U^H U
		complex(dp),	intent(in)		:: U(:,:)
		complex(dp),	allocatable		:: I(:,:)
		complex(dp)						:: alpha, beta
		integer							:: m, n, k, lda, ldb, ldc
		character*1						:: transa, transb
		!
		!
		isUnit = .false.
		n	= size(U,1)
		m	= size(U,2)
		!
		if(n /= m )	then
			write(*,*)"[isUnit]: matrix is not even a square matrix ^^"
		else
			allocate(	I( n,n)	)
			!
			transa	= 'c'
			transb	= 'n'
			k		= size(I,1)
			lda		= size(U,1)
			ldb		= size(U,2)
			ldc		= size(I,1)
			alpha	= dcmplx(1.0_dp)
			beta	= dcmplx(0.0_dp)
			call zgemm(transa, transb, m, n, k, alpha, U, lda, U, ldb, beta, I, ldc)
			!
			isUnit = isIdentity(I)
		end if
		!
		return
	end function


	logical function isIdentity(I)
		!test if I is Identity matrix
		complex(dp),	intent(in)		:: I(:,:)
		integer							:: n,m, Isize
		!
		isIdentity 	= .true.
		!
		if(size(I,1) /= size(I,2) ) then
			write(*,*)"[isIdentity]: ERROR - only implemented for square matrices (will return false)"
			isIdentity = .false.
		else
			Isize	= size(I,1)
			m = 1
			do while( m<= Isize .and. isIdentity )
				n = 1
				do while( n<= Isize .and. isIdentity )
					!
					!
					if(n == m) then 
						if(  abs(abs(dreal(I(n,n)))-1.0_dp) > acc  .or. abs(dimag(I(n,n))) > acc	) 	isIdentity = .false.
					else 
						if(	abs(I(n,m)) > acc ))														isIdentity = .false.
					end if
					n = n + 1
					!
					!	
				end do
				m = m + 1
			end do 
		end if
		!
		return
	end function


	logical function isHermitian(H)
		!tests if Matrix is hermitian
		complex(dp),	intent(in)		:: H(:,:)
		integer							:: n,i,j
		!
		isHermitian = .true.
		if(size(H,1) /= size(H,2)) then
			write(*,*)"[isHermitian]: matrix is not a square matrix ^^"
			isHermitian = .false.
		else
			n	= size(H,1)
			i	= 1
			do while ( i <= n .and. isHermitian) 
				j = 1
				do while ( j <= n .and. isHermitian) 
					if(		abs( H(i,j) - dconjg(H(j,i)) ) 	> acc		) 	isHermitian = .false.
					j = j + 1
				end do
				i = i + 1 
			end do
		end if
		return
	end function


	subroutine rotMat(U, Mat, res)
		! res = U^dagger * Mat * U
		!	1. res = Mat * U
		!	2. res = U^da * res
		complex(dp),	intent(in)		:: U(:,:), Mat(:,:)
		complex(dp),	intent(out)		:: res(:,:)
		integer							:: m, n, k, lda, ldb, ldc
		complex(dp)						:: alpha, beta
		character*1						:: transa, transb
		!
		transa	= 'n'
		transb	= 'n'
		m		= size(Mat,1)
		n		= size(U,2)
		k		= size(Mat,2)
		lda		= m
		ldb		= k
		ldc		= m
		alpha	= dcmplx(1.0_dp)
		beta	= dcmplx(0.0_dp)
		call zgemm(transa, transb, m, n, k, alpha, Mat, lda, U, ldb, beta, res, ldc)
		!
		transa	= 'c'
		transb	= 'n'
		call zgemm(transa, transb, m, n, k, alpha, U, lda, res, ldb, beta, res, ldc)
		!
		!
		return
	end subroutine


	subroutine myCommutat(A, B, res)
		!	computes the commutator of matrix M and N
		!	RES	= [A,B] = AB - BA
		!	1. RES = BA
		!	2. RES = AB - RES
		complex(dp),	intent(in)		:: A(:,:), B(:,:)
		complex(dp),	intent(out)		:: RES(:,:)
		integer							:: m, n, k, lda, ldb, ldc
		complex(dp)						:: alpha, beta
		character*1						:: transa, transb
		!
		if( size(A,2) /= size(B,2)  .or. size(B,1) /= size(A,1)	) then
			write(*,*)"[myCommutat]: Matrix ranks dont match. can not compute commutator"
		else
			transa	= 'n'
			transb	= 'n'
			m 		= size(B,1)
			n 		= size(A,2)
			k		= size(A,1)
			lda		= m
			ldb		= k
			ldc		= m
			alpha	= dcmplx(1.0_dp)
			beta	= dcmplx(0.0_dp)
			call zgemm(transa, transb, m, n, k, alpha, B, lda, A, ldb, beta, RES, ldc)
			!
			beta	= dcmplx(-1.0_dp)
			call zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, RES, ldc)
		end if
		!
		!
		return
	end subroutine





























!private:  



	complex(dp)	function nIntegrateCPLX(f)
		complex(dp),	intent(in)		:: f(:)				!f(nR)
		real(dp)						:: rel, img
		!
		rel = nIntegrateREAL(dreal(f))
		img = nIntegrateREAL(dimag(f))
		!
		nIntegrateCPLX = dcmplx(	rel , img	)
		!
		return
	end function


	real(dp) function nIntegrateREAL(f)
		!2D integration routine
			!integrates for each y point integrate over whole x direction
		real(dp),		intent(in)		:: f(:)				!f(nR)
		

		nIntegrateREAL = sum(f) / real(size(f),dp)
		!
		return
	end function

	!real(dp) function nIntegrateREAL(nR, nRx,nRy, dx,dy, f)
	!	!2D integration routine
	!		!integrates for each y point integrate over whole x direction
	!	integer,		intent(in)		:: nR, nRx, nRy
	!	real(dp),		intent(in)		:: dx, dy
	!	real(dp),		intent(in)		:: f(:)				!f(nR)
	!	!
	!	nIntegrateREAL = sum( f(:) ) / real(size(f),dp)
	!	!
	!	return
	!end function


	function crossPreal(a,b)
		!cross product of two real 3dim vectors a,b
		real(dp),		intent(in)	:: a(3), b(3)
		real(dp)					:: crossPreal(3)
		!
		crossPreal(1)	=	a(2) * b(3)		-	a(3) * b(2)	 
		crossPreal(2)	=	a(3) * b(1)		-	a(1) * b(3) 
		crossPreal(3)	=	a(1) * b(2)		-	a(2) * b(1) 
		!
		return
	end function


	function crossPcplx(a,b)
		!cross product of two complex 3dim vectors a,b
		complex(dp),	intent(in)	:: a(3), b(3)
		complex(dp)					:: crossPcplx(3)
		!
		crossPcplx(1)	=	a(2) * b(3)		-	a(3) * b(2)	 
		crossPcplx(2)	=	a(3) * b(1)		-	a(1) * b(3) 
		crossPcplx(3)	=	a(1) * b(2)		-	a(2) * b(1) 
		!
		return
	end function





	subroutine hiordqWrapper(dx,f,res)
		!uses intlib routine for equally spaced data
		real(dp), 	intent(in)		:: dx, f(:)
		real(dp),	intent(out)		:: res
		integer						:: ntab
		real(dp), 	allocatable		:: work(:)
		ntab = size(f)
		allocate(	work(2*(ntab-1))	)
		!
		call  hiordq( ntab, dx, f, work, res )
		return
	end subroutine


	!EIGSOLVER HELPERS
    subroutine errCheck(n,info, jobz)
		!small subroutine used in subroutine eigSolver
		implicit none
		integer    , intent(in) :: n,info
		character*1, intent(in) :: jobz
		!
		if(info > 0) then
		      write(*,*) '[solver/errCheck]: WARNING, Problem solving the eigenvalue problem: '
		      if(jobz .EQ. 'N') write(*,*) '[solver/errCheck]: the algorithm failed to converge; ', info ,&
		      				' off-diagonal elements of an intermediate tridiagonal form did not converge to zero;'
		      if(jobz .EQ. 'V') write(*,*) '[solver/errCheck]: the algorithm failed to compute an',&
		    ' eigenvalue while working on the submatrix lying in rows and columns ', info/(n+1), ' through ', mod(info, n+1)
		elseif(info < 0) then
		     write(*,*) '[solver/errCheck]: the ',info,'-th parameter had an illegal value.' 
		endif
		!
		return
    end subroutine


	subroutine normalCheck(n,eigVec)
    	!check if eigVec are orthonormal this should always be the case
		integer     , intent(in)      :: n
		complex(dp) , intent(in)      :: eigVec(:,:)
		integer                       :: i1,i2
		real(dp)                      :: tmp
		do i2=1,n
		      tmp = 0.0_dp
		      do i1=1,n
		            tmp = tmp + dconjg(eigVec(i1,i2))*eigVec(i1,i2)
		      end do
		      if(abs(tmp)-1.0_dp>1e-12_dp) then
		            write(*,*)"[solver/normalCheck]: WARNING eigVec not normal: ",dsqrt(tmp)
		      end if
		end do
		!
		return
	end subroutine





end module mathematics