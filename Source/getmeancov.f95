!	Find means and covariance matrices corresponding to a predictor list for the normal case.

!	beta is the parameter vector.
!	predlist is the predictor list.
!	tolsing is the tolerance for modified Cholesky decomposition.
!	covs are the covariance matrices.
!	means are the means.

subroutine getmeancov(beta,predlist,tolsing,covs,means)
implicit none
interface
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
		real(kind=8),intent(in)::beta(:),thpr(:),tolsing
		real(kind=8),intent(out)::covmatrix(:,:),meanvec(:)
	end subroutine meancovstat
end interface


real(kind=8),intent(in)::beta(:),predlist(:,:),tolsing
real(kind=8),intent(out)::covs(:,:,:),means(:,:)
!	row counts entries in predlist.
integer::row
do row=1,size(predlist,2)
	call meancovstat(beta,predlist(:,row),tolsing,covs(:,:,row),means(:,row))
end do
return
end subroutine getmeancov