!	Obtain mean item scale parameters for an observation.
!	numcat is a vector with the number of
!		underlying categories per item.
!	The items presented are indicated by mask.
!	The underlying item probabilities are in probcat.
!	scales contains the scale parameters.
function mean(numcat,mask,probcat,scales)
implicit none
integer,intent(in)::numcat(:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::probcat(:),scales(:,:)
real(kind=8)::mean(size(scales,1),size(mask))
!	item is the item number. nloc points to the item probability.
integer::item,nloc

nloc=1
do item=1,size(mask)
	if(mask(item))then
		mean(:,item)=matmul(scales(:,nloc:nloc+numcat(item)-1),probcat(nloc:nloc+numcat(item)-1))
	end if
	nloc=nloc+numcat(item)
end do
return
end function mean
