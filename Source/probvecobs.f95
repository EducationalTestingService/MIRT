!	probvecobs is used to compute observed item probabilities given the latent vector.
!	catobsrange is the table of ranges of underlying categories per observed category.
!	numcatobs provides the number of observed categories per item.
!	resp is the response vector.
!	mask indicates which items were presented.
!	probcat is the vector of underlying item probabilities.
function probvecobs(catobsrange,numcatobs,resp,mask,probcat)
implicit none

integer,intent(in)::catobsrange(:,:),numcatobs(:),resp(:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::probcat(:)
real(kind=8)::probvecobs(size(mask))
!
!	item is the item number.	
!	nlocobs is the location of the observed item category.
!	position locates the item response.
integer::item,nlocobs,position
probvecobs=0.0_8
nlocobs=1
do item=1,size(mask)
	if(mask(item))then
		position=nlocobs+resp(item)
		probvecobs(item)=sum(probcat(catobsrange(1,position):catobsrange(2,position)))
		
	end if
	nlocobs=nlocobs+numcatobs(item)
end do

return
end function probvecobs


