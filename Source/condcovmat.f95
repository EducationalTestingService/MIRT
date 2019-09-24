!	Obtain conditional covariance matrix of item scale parameters for an observation.
!	catobsrange contains the range of underlying categories per observed category.
!	numcatobs is a vector with the number of observed  categories per item.
!	resp is the observation vector.
!	The items presented are indicated by mask.
!	condmeans contains the conditional means.
!	The underlying item conditional probabilities are in condprobcat.
!	The scale parameters are in scales.

function condcovmat(catobsrange,numcatobs,resp,mask,condmeans,condprobcat,scales)
implicit none
integer,intent(in)::numcatobs(:),catobsrange(:,:),resp(:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::condmeans(:,:),condprobcat(:),scales(:,:)
real(kind=8)::condcovmat(size(scales,1),size(scales,1),size(resp))
!	cat corresponds to an underlying category.
!	col is the column number.
!	item is the item number.
!	nlocobs points to the observed item category.
!	position points to the underlying categories for the item response.
!	row is the row number.
!
integer::cat,col,item,nlocobs,position,row
!	diff is the difference between the scale parameters and their conditional means.
real(kind=8)::diff(size(scales,1))
condcovmat=0.0_8
nlocobs=1
do item=1,size(scales,1)
	if(mask(item))then
		position=nlocobs+resp(item)
	
		do cat=catobsrange(1,position),catobsrange(2,position)
			diff=scales(:,cat)-condmeans(:,item)
			do row=1,size(scales,1)
				do col=1,row
					condcovmat(row,col,item)=condcovmat(row,col,item)&
						+diff(row)*diff(col)*condprobcat(cat)
					
					
				end do
				
			end do
			
		end do
	end if

	nlocobs=nlocobs+numcatobs(item)
end do
return
end function condcovmat
