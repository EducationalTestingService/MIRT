!	Get multinomial probabilities for latent vector.  Place in densitym.
!	Input involves linear component linth, quadrature points quadpoint, quadratic components quadth, and quadrature weights quadweight.
function densitym(linth,quadpoint,quadth,quadweight)
implicit none
real(kind=8),intent(in)::linth(:),quadpoint(:,:),quadth(:,:),quadweight(:)
real(kind=8)::densitym(size(quadweight))
!	quad counts quadrature points.
integer::quad
!	rescale density with rescale.
real(kind=8)::rescale

do quad=1,size(quadweight)
	densitym(quad)=quadweight(quad)*exp(dot_product(quadpoint(:,quad),linth+matmul(quadth,quadpoint(:,quad))))
end do
rescale=1.0_8/sum(densitym)
densitym=rescale*densitym
return
end function densitym