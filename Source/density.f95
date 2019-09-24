!	Get normal density for scaled latent vector.  Place in density.
!	Input involves linear component linth, quadrature points quadpoint, quadratic components quadth, quadrature weights quadweight, and scale
!	factor scaletheta.
function density(linth,quadpoint,quadth,quadweight,scaletheta)
implicit none
real(kind=8),intent(in)::linth(:),quadpoint(:,:),quadth(:,:),quadweight(:),scaletheta
real(kind=8)::density(size(quadweight))
!	quad counts quadrature points.
integer::quad
!	rescale density with rescale.


do quad=1,size(quadweight)
	density(quad)=scaletheta*quadweight(quad)*exp(dot_product(quadpoint(:,quad),linth+matmul(quadth,quadpoint(:,quad))))
end do

return
end function density
