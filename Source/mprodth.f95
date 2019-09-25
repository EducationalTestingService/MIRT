!	Means of products.  meanth contains means.
!	covth contains covariances.
!	quadpoint contains points.  density contains weights.
function mprodth(covth,meanth)
implicit none
real(kind=8),intent(in)::covth(:,:),meanth(:)
real(kind=8)::mprodth(size(meanth),size(meanth))
!	col for column number.
!	row for row number.
integer::col,row
do row=1,size(meanth)
	do col=1,row
		mprodth(row,col)=covth(row,col)+meanth(row)*meanth(col)
	end do
end do
return
end function mprodth
