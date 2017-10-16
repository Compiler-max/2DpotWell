module tightBind
	use omp_lib
	use mathematics,	only:	dp, PI_dp,i_dp, machineP, myExp, myLeviCivita, eigSolver, nIntegrate, isUnit, isHermitian
	use sysPara
	use blochWf,		only:	genBwfVelo, testNormUNK							
	use output,			only:	printMat, writeEnAndUNK
	implicit none	
	
	private
	public ::					solveELECTRstruct(unk, En)



		

	


	contains
	subroutine solveELECTRstruct(unk, En)

	end subroutine

end module tightBind