!	Find means and covariance matrices corresponding to a predictor list for the multinomial case.

!	beta is the parameter vector.
!	predlist is the predictor list.
!	quadpoint is the array of quadrature points.
!	quadweight is the array of quadrature weights.

!	covs are the covariance matrices.
!	means are the means.

subroutine getmeancovm(beta,predlist,quadpoint,quadweight,covs,means)
implicit none
interface
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
		real(kind=8),intent(in)::beta(:),quadpoint(:,:),quadweight(:),thpr(:)
		real(kind=8),intent(out)::covmatrix(:,:),meanvec(:)
	end subroutine meancovstatm

end interface


real(kind=8),intent(in)::beta(:),predlist(:,:),quadpoint(:,:),quadweight(:)
real(kind=8),intent(out)::covs(:,:,:),means(:,:)
!	row counts entries in predlist.
integer::row
do row=1,size(predlist,2)
	call meancovstatm(beta,quadpoint,quadweight,predlist(:,row),covs(:,:,row),means(:,row))
end do
return
end subroutine getmeancovm
