!	Find the gradient of the conditional log posterior.
!	mask indicates which items are present.
!	condmeans is the vector of conditional means of scale factors.
!	linth is the linear component of the log density.
!	lintran is the linear transformation of the latent vector.
!	means is the vector of means of scale factors.
!	quadth is the quadratic component of the log density.
!	theta is the value of the latent vector.
function gradtheta(mask,condmeans,linth,lintran,means,quadth,theta)
implicit none

logical,intent(in)::mask(:)
real(kind=8),intent(in)::condmeans(:,:),linth(:),lintran(:,:),means(:,:),quadth(:,:),theta(:)
real(kind=8)::gradtheta(size(theta))
real(kind=8)::itempart(size(lintran,1))
!	item is an item
integer::item
itempart=0.0_8
do item=1,size(mask)
	if(mask(item))then
		itempart=itempart+(condmeans(:,item)-means(:,item))
	end if
end do
gradtheta=linth+2.0_8*matmul(quadth,theta)+matmul(transpose(lintran),itempart)
return
end function gradtheta
	
