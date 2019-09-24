subroutine makesym(matrix)
!	Symmetrize n by n matrix based on lower triangle.
implicit none
real(kind=8),intent(inout)::matrix(:,:)
!	col is column.
!	n is matrix size.
!	row is row.
integer::row,col,n
n=size(matrix,1)
if(n>1) then
	do row=2,n
		do col=1,row-1
			matrix(col,row)=matrix(row,col)
		end do
	end do
end if
return
end subroutine makesym			
