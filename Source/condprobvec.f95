!	condprobvec is used to compute observed item conditional probabilities of underlying item categories given
!	observed item categories.
!	catobsrange provides ranges of underlying item categories that correspond to observed item categories.
!	numcatobs provides the number of observed categories per item.
!	resp is the item response.
!	mask indicates which items were presented.
!	probcat is the vector of underlying category probabilities.
!	probcatobs is the vector of observed category probabilities.
function condprobvec(catobsrange,numcatobs,resp,mask,probcat,probcatobs)
implicit none
integer,intent(in)::catobsrange(:,:),numcatobs(:),resp(:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::probcat(:),probcatobs(:)
real(kind=8)::condprobvec(size(probcat))
!

!	item is the item number.
!	nlocobs is the location of the observed item category.
!	position locates underlying categories for a response.
!
integer::nlocobs,item,position
nlocobs=1
condprobvec=0.0_8
do item=1,size(mask)
	if(mask(item))then
		position=nlocobs+resp(item)
		condprobvec(catobsrange(1,position):catobsrange(2,position))&
			=probcat(catobsrange(1,position):catobsrange(2,position))/probcatobs(item)
	end if
	nlocobs=nlocobs+numcatobs(item)
end do
return
end function condprobvec
