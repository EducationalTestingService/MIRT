!	Find the negative hessian of the conditional log posterior of an observation
!	relative to the latent vector theta.
!	mask indicates which items are present.
!	condcovs is the arrayr of conditional covariance matrics of scale factors.
!	covs is the array of covariance matrices of scale factors.
!	lintran is the linear transformation of the latent vector.
!	quadth is the matrix parameter for the log density. 

function nhesstheta(mask,condcovs,covs,lintran,quadth)
implicit none
interface
	subroutine makesym(matrix)
!	Symmetrize n by n matrix based on lower triangle.
		implicit none
		real(kind=8),intent(inout)::matrix(:,:)
	end subroutine makesym
end interface


logical,intent(in)::mask(:)
real(kind=8),intent(in)::condcovs(:,:,:),covs(:,:,:),lintran(:,:),quadth(:,:)
real(kind=8)::nhesstheta(size(quadth,1),size(quadth,1))

!	col is a column number
!	item is item number.
!	row is a row number
integer::col,item,row
!	Accumulate in itempart.
real(kind=8)::itempart(size(lintran,1),size(lintran,1))

itempart=0.0_8

do item=1,size(mask)
	if(mask(item))then
		do row=1,size(lintran,1)
			do col=1,row
				itempart(row,col)=itempart(row,col)+(covs(row,col,item)-condcovs(row,col,item))
			end do
		end do
	end if
end do
call makesym(itempart)
nhesstheta=-2.0_8*quadth+matmul(transpose(lintran),matmul(itempart,lintran))
return
end function nhesstheta
	
