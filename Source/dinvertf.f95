!	Compute inverse for n by n matrix with modified Cholesky
!	decomposition tri.
function invert(tri)
implicit none
real(kind=8),intent(in)::tri(:,:)
real(kind=8)::invert(size(tri,1),size(tri,1))
integer::row,col,n
n=size(tri,1)
!	Invert diagonals.
do  row=1,n
	if(tri(row,row)>0.0_8) invert(row,row) = 1.0_8/tri(row,row)
	if(tri(row,row)<=0.0) invert(row,row) = 0.0_8
end do
if(n==1) return
!	Off-diagonals.
do  row=2,n
	do  col=1,row-1
		invert(row,col) = -tri(row,col)
	 	if(col<row-1)invert(row,col)=invert(row,col)-dot_product(tri(row,col+1:row-1),invert(col+1:row-1,col))
		invert(col,row) = invert(row,col)*invert(row,row)
	end do
end do
!	Full inverse
do  row=1,n
	do  col=1,row
		invert(row,col) = invert(col,row)
		if(row<n) invert(row,col)=invert(row,col)+dot_product(invert(row+1:n,row),invert(col,row+1:n))
	end do
end do
do row=2,n
	do col=1,row-1
		invert(col,row) = invert(row,col)
	end do
end do
return
end function invert
