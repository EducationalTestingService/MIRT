!	probvec is used to compute underlying item probabilities.
!   	numcat indicates the number of underlying categories per item.
!	mask indicates which items were presented.
!	locations provides location parameters and scales provides scale factors.
!	theta is the value of the latent vector.
function probvec(numcat,mask,locations,scales,theta)
implicit none
interface
!	probvecitem is used to compute underlying item probabilities for one probability.
!	locations provides location parameters and scales provides scale factors.
!	theta is the value of the latent vector.
	function probvecitem(locations,scales,theta)
		implicit none


		real(kind=8),intent(in)::locations(:),scales(:,:),theta(:)
		real(kind=8)::probvecitem(size(locations))
	end function probvecitem
end interface
integer,intent(in)::numcat(:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::locations(:),scales(:,:),theta(:)
real(kind=8)::probvec(size(locations))

!	

!	item is the item number.

!	nloc1 is the first location of the item category and location parameter.
!	nloc2 is the end location of the item category and location parameter.

integer::item,nloc1,nloc2

!	Initialize locations.

nloc1=1
nloc2=0
probvec=0.0_8
do item=1,size(mask)
	nloc2=nloc2+numcat(item)
	if(mask(item))probvec(nloc1:nloc2)=probvecitem(locations(nloc1:nloc2),scales(:,nloc1:nloc2),theta)
	nloc1=nloc2+1
end do
return
end function probvec

