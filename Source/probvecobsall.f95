!	probvecobsall is used to compute observed item probabilities given the latent vector.
!	catobsrange is the table of ranges of underlying categories per observed category.
!	numcat provides the number of underlying categories per item, and numcatobs provides the number
!		of observed categories per item.
!	mask indicates which items were presented.
!	probcat is the vector of underlying item probabilities.
function probvecobsall(catobsrange,numcatobs,mask,probcat)
implicit none
integer,intent(in)::catobsrange(:,:),numcatobs(:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::probcat(:)
real(kind=8)::probvecobsall(sum(numcatobs))
!	
!	nlocobs is the location of the observed item category.


!	item is the item number.
!	position locates the item response.
integer::nlocobs,item,position
probvecobsall=0.0_8
nlocobs=1
do item=1,size(mask)
	
	if(mask(item))then
		do position=nlocobs,nlocobs+numcatobs(item)-1
			probvecobsall(position)=sum(probcat(catobsrange(1,position):catobsrange(2,position)))
		end do
	end if
	nlocobs=nlocobs+numcatobs(item)
end do

return
end function probvecobsall
