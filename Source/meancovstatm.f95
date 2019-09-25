!
!  Find statistics of mean and covariance matrix of underlying latent vector.
!
!   beta is the estimated parameter vector.
!	quadpoint is the array of quadrature points.
!	quadweight is the array of quadrature weights.
!   thpr is a predictor vector.
!   covmatrix is the estimated covariance matrix of the underlying latent vector corresponding to beta.
!   meanvec is the estimated mean of the underlying latent vector.
!   Note that this routine is for the case of an underlying multinomial latent vector.
subroutine meancovstatm(beta,quadpoint,quadweight,thpr,covmatrix,meanvec)
implicit none
interface

!	Get multinomial probabilities for latent vector.  Place in density.
!	Input involves linear component linth, quadrature points quadpoint, quadratic components quadth, and quadrature weights quadweight.
	function densitym(linth,quadpoint,quadth,quadweight)
		real(kind=8),intent(in)::linth(:),quadpoint(:,:),quadth(:,:),quadweight(:)
		real(kind=8)::densitym(size(quadweight))
	end function densitym

!	convert an n by n symmetric matrix from compact to regular form.
!
	function expandmat(n,compmat)
		implicit none
		integer,intent(in)::n
		real(kind=8),intent(in)::compmat(:)
		real(kind=8)::expandmat(n,n)
	end function expandmat

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


!	Multinomial covariance matrix.
!	density is vector of probabilities.
!	meanth is vector of means.
!	quadpoint is array of weights.
	function mcovth(density,meanth,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),meanth(:),quadpoint(:,:)
		real(kind=8)::mcovth(size(meanth),size(meanth))

	end function mcovth
	
!	Multinomial mean.
!	density is vector of probabilities
!	quadpoint is array of weights.
	function mmeanth(density,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),quadpoint(:,:)
		real(kind=8)::mmeanth(size(quadpoint,1))
	end function mmeanth

end interface
real(kind=8),intent(in)::beta(:),quadpoint(:,:),quadweight(:),thpr(:)
real(kind=8),intent(out)::covmatrix(:,:),meanvec(:)


!   linth is the linear commponent of beta for thpr.
!   lintheta is a matrix version of the linear component of beta.
!   quadth is the quadratic component of beta for thpr.
!   quadtheta is an array version of the quadratic component of beta.

real(kind=8)::density(size(quadweight)),linth(size(covmatrix,1)), lintheta(size(covmatrix,1),size(thpr)),&
	quadth(size(covmatrix,1),size(covmatrix,1)),quadtheta(size(covmatrix,1),size(covmatrix,1),size(thpr))
!   Find sizes.



!   Linear and quadratic components of beta.
call getlinquad(beta,lintheta,quadtheta)

!   Linear and quadratic components of beta for thpr.
call getlinquadth(lintheta,quadtheta,thpr,linth,quadth)
density=densitym(linth,quadpoint,quadth,quadweight)





!   Mean of theta.
meanvec=mmeanth(density,quadpoint)
!	Covariance matrix of theta
covmatrix=mcovth(density,meanvec,quadpoint)
return
end subroutine meancovstatm
