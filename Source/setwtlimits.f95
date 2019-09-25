!	Find weight limits for each item.
!	numcatobs gives number of observed categories per item.
!	weights are category weights.
!	maxw gives the maximum weights per item.
!	minw gives the minimum weights per item.
subroutine setwtlimits(numcatobs,weights,maxw,minw)
implicit none
integer,intent(in)::numcatobs(:),weights(:)
integer,intent(out)::maxw(:),minw(:)
!	counter counts categories.
!	item counts items.
integer::counter,item
counter=1
do item=1,size(numcatobs)
	maxw(item)=maxval(weights(counter:counter+numcatobs(item)-1))
	minw(item)=minval(weights(counter:counter+numcatobs(item)-1))
	counter=counter+numcatobs(item)
end do
return
end subroutine setwtlimits