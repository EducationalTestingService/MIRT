!	Find a substitute negative hessian nhessth of the conditional log posterior of an observation
!	relative to the latent vector theta.
!	mask indicates which items are present.
!	condmeans is the vector of conditional means of scale factors.
!	gradth is the gradient.
!	lintran is the linear transformation of the latent vector.
!	means is the vector of means of scale factors.
!	quadth is the matrix of parameters used in the log density.

function nhessthsub(mask,condmeans,gradth,lintran,means,quadth)
implicit none
interface
	subroutine makesym(matrix)
!	Symmetrize n by n matrix based on lower triangle.
		implicit none
		real(kind=8),intent(inout)::matrix(:,:)
	end subroutine makesym
end interface

logical,intent(in)::mask(:)
real(kind=8),intent(in)::condmeans(:,:),gradth(:),lintran(:,:),means(:,:),quadth(:,:)
real(kind=8)::nhessthsub(size(means,1),size(means,1))
!	col is a column number.
!	item is item number.
!	row is a row number.
integer::col,item,row
!	average gradient in avegradth.
!	difference in diff.
!	Accumulate in itempart.
real(kind=8)::avegradth(size(gradth)),diff(size(gradth)),&
	itempart(size(lintran,1),size(lintran,1))

itempart=0.0_8	

	


do item=1,size(mask)
	if(mask(item))then
		diff=condmeans(:,item)-means(:,item)
		do row=1,size(means,1)
			do col=1,row
				itempart(row,col)=itempart(row,col)+diff(row)*diff(col)
			end do
		end do
	end if
end do
call makesym(itempart)
nhessthsub=-2.0_8*quadth+matmul(transpose(lintran),matmul(itempart,lintran))	
end function nhessthsub
