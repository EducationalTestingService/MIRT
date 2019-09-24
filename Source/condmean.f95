!	Obtain conditional mean item scale parameters.
!	catobsrange contains the range of underlying categories per observed category.
!	numcatobs is a vector with the number of observed categories per item.
!	The observation is resp.
!	The items presented are indicated by mask.
!	The underlying conditional item probabilities are in condprobcat.
!	scales is the matrix of scale parameters. 
function condmean(catobsrange,numcatobs,resp,mask,condprobcat,scales)
implicit none
integer,intent(in)::catobsrange(:,:),numcatobs(:),resp(:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::condprobcat(:),scales(:,:)
real(kind=8)::condmean(size(scales,1),size(resp))
!	end completes a range of categories.
!	item is the item number.
!	nlocobs points to the observed  category.
!	position points to the underlying categories for the item.
!	start begins a range of categories.
integer::end,item,nlocobs,position,start

nlocobs=1
do item=1,size(resp)
	if(mask(item))then
		position=nlocobs+resp(item)
		start=catobsrange(1,position)
		end=catobsrange(2,position)
		condmean(:,item)=matmul(scales(:,start:end),condprobcat(start:end))
	end if
	nlocobs=nlocobs+numcatobs(item)

end do
return
end function condmean
