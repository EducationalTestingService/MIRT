!	Extract the linear component linth and the quadratic component quadth from lintheta and quadtheta for a given predictor thpr.
subroutine getlinquadth(lintheta,quadtheta,thpr,linth,quadth)
implicit none


real(kind=8),intent(in)::lintheta(:,:),quadtheta(:,:,:),thpr(:)
real(kind=8),intent(out)::linth(:),quadth(:,:)
!	row counts elements of the underlying latent vector.
integer::row
linth=matmul(lintheta,thpr)
do row=1,size(linth)
	quadth(row,:)=matmul(quadtheta(row,:,:),thpr)
end do
return
end subroutine getlinquadth