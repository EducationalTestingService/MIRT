!	probvecitem is used to compute underlying item probabilities for one probability.

!	locations provides location parameters and scales provides scale factors.
!	theta is the value of the latent vector.
function probvecitem(locations,scales,theta)
implicit none


real(kind=8),intent(in)::locations(:),scales(:,:),theta(:)
real(kind=8)::probvecitem(size(locations))

!	
!	cat is the category number.



integer::cat
!	probdiv is used to normalize probabilities	
real(kind=8)::probdiv
!	Initialize locations.



do cat=1,size(locations)
	probvecitem(cat)=exp(sum(theta*scales(:,cat))+locations(cat))
	
end do
probdiv=sum(probvecitem)
probvecitem=probvecitem/probdiv

return
end function probvecitem
	
