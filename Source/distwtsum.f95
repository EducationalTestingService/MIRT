!	Obtain the frequency distribution of a weighted sum of responses.
!	The weights are integers.
!	maxnw is the array of minimum scores for items.
!	minw is the array of maximum scores for items.
!	numcatobs is the number of observed item categories.
!	weight is the category weights.
!	prob is the vector of observed category probabilities.
!	dist is the probability distribution of the weighted sum.

subroutine distwtsum(maxw,minw,numcatobs,weight,prob,dist)
implicit none
integer,intent(in)::maxw(:),minw(:),numcatobs(:),weight(:)
real(kind=8),intent(in)::prob(:)
real(kind=8),intent(out)::dist(sum(minw):sum(maxw))

!	cat is a category.
!	cat1 is another category.
!	counter provides the probability position.
!	item is an item.
!	nitems is the number of items.
!	top is the maximum value up to item i
integer::cat,cat1,counter,item,nitems,top
!	distnew is the probability distribution for the first i-1 items.
!	distold is the probability distribution for the first i items.

real(kind=8)::distnew(0:sum(maxw)-sum(minw)),distold(0:sum(maxw)-sum(minw))	



distnew=0.0_8
dist=0.0_8
distnew(0)=1.0_8

counter=0
top=0
nitems=size(numcatobs)
do item=1,nitems
	
   if(maxw(item)>minw(item)) then
		distold(0:top)=distnew(0:top)
		distnew(0:top+maxw(item)-minw(item))=0.0_8
		do cat=0,top
			do cat1=1,numcatobs(item)
				distnew(cat+weight(counter+cat1)-minw(item))=distnew(cat+weight(counter+cat1)-minw(item))+distold(cat)*prob(counter+cat1)
			end do
		end do
		top=top+maxw(item)-minw(item)
	end if
	counter=counter+numcatobs(item)
end do
dist(sum(minw):sum(maxw))=distnew

return

end subroutine distwtsum


