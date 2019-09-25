!
!  Find statistics of mean and covariance matrix of underlying latent vector.
!
!   beta is the estimated parameter vector.

!   thpr is a predictor vector.
!   tolsing is a tolerance for modified Cholesky decomposition.

!   covmatrix is the estimated covariance matrix of the underlying latent vector corresponding to beta.
!   meanvec is the estimated mean of the underlying latent vector.
!   Note that this routine is for the case of an underlying normal latent vector.
subroutine meancovstat(beta,thpr,tolsing,covmatrix,meanvec)
implicit none
interface
!	chol is used to compute a modified Cholesky decomposition of the n by n matrix sym.
!	The singularity tolerance is tolsing.
	
	function chol(sym,tolsing)
		implicit none
		real(kind=8),intent(in)::sym(:,:),tolsing
	 	real(kind=8)::chol(size(sym,1),size(sym,1))
	end function chol




!	Extract the linear component lintheta and the quadratic component quadtheta from the parameter vector beta.
	subroutine getlinquad(beta,lintheta,quadtheta)
		implicit none
		real(kind=8),intent(in)::beta(:)
		real(kind=8),intent(out)::lintheta(:,:),quadtheta(:,:,:)
	end subroutine getlinquad

!	Extract the linear component linth and the quadratic component quadth from lintheta and quadtheta for a given predictor thpr.
	subroutine getlinquadth(lintheta,quadtheta,thpr,linth,quadth)
		implicit none
		real(kind=8),intent(in)::lintheta(:,:),quadtheta(:,:,:),thpr(:)
		real(kind=8),intent(out)::linth(:),quadth(:,:)
	end subroutine getlinquadth

	
!	Compute inverse for n by n matrix with modified Cholesky
!	decomposition tri.
	function invert(tri)
		implicit none
		real(kind=8),intent(in)::tri(:,:)
		real(kind=8)::invert(size(tri,1),size(tri,1))
	end function invert
!		This program uses the modified Cholesky decomposition
!		in lu to solve the equation lusolve = b.


	function solve(lu,b)
		implicit none
		
		real(kind=8),intent(in)::lu(:,:)
		real(kind=8),intent(in)::b(:)
		real(kind=8)::solve(size(b))
	end function solve

end interface
real(kind=8),intent(in)::beta(:),thpr(:),tolsing
real(kind=8),intent(out)::covmatrix(:,:),meanvec(:)

!	covinverse is the inverse of covmatrix.
!   linth is the linear commponent of beta for thpr.
!   lintheta is a matrix version of the linear component of beta.
!   quadth is the quadratic component of beta for thpr.
!   quadtheta is an array version of the quadratic component of beta.

real(kind=8)::covinverse(size(covmatrix,1),size(covmatrix,1)),linth(size(covmatrix,1)), lintheta(size(covmatrix,1),size(thpr)),&
	quadth(size(covmatrix,1),size(covmatrix,1)),quadtheta(size(covmatrix,1),size(covmatrix,1),size(thpr))


!   Linear component of beta.
call getlinquad(beta,lintheta,quadtheta)


!   Linear component of beta for thpr.
call getlinquadth(lintheta,quadtheta,thpr,linth,quadth)

covinverse=-2.0_8*quadth

covinverse=chol(covinverse,tolsing)
!   Mean of theta.
meanvec=solve(covinverse,linth)
!	invert
covmatrix=invert(covinverse)

end subroutine meancovstat
