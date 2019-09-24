!	Obtain covariance matrix of item scale parameters for an observation.
!	numcat is a vector with the number of underlying categories per category.
!	The items presented are indicated by mask.
!	means contains the corresponding means.
!	The underlying item probabilities are in probcat.
!	scales contains the scale parameters.   

function covmat(numcat,mask,means,probcat,scales)
implicit none
integer,intent(in)::numcat(:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::means(:,:),probcat(:),scales(:,:)
real(kind=8)::covmat(size(scales,1),size(scales,1),size(mask))
!	cat is the category number.
!	col is the column number.
!	item is the item number. 	
!	nloc points to the item probability.
!	row is the row number. 
integer::cat,col,item,nloc,row
!	diff contains the deviation of the scale parameters from their means.
real(kind=8)::diff(size(scales,1))
covmat=0.0_8
nloc=1
do item=1,size(mask)
	if(mask(item))then
		do cat=nloc,nloc+numcat(item)-1
			diff=scales(:,cat)-means(:,item)
			do row=1,size(diff)
				do col=1,row
					covmat(row,col,item)=covmat(row,col,item)&
						+diff(row)*diff(col)*probcat(cat)
				end do
			end do
		end do
	endif
	nloc=nloc+numcat(item)
end do	
return
end function covmat
