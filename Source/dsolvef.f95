
!		This program uses the modified Cholesky decomposition
!		in lu to solve the equation lusolve = b.
function solve(lu,b)
implicit none

real(kind=8),intent(in)::lu(:,:)
real(kind=8),intent(in)::b(:)
real(kind=8)::solve(size(b))
integer:: row,col,n
n=size(b)
!	 Lower triangle.
solve(1) = b(1)
if(lu(1,1)<=0.0_8) solve(1) = 0.0_8
if(n>1)then
	do row=2,n
		if(lu(row,row)<=0.0_8) then
			solve(row)=0.0_8
		else
			solve(row)=b(row)-dot_product(lu(row,1:row-1),solve(1:row-1))
	 	end if
	 end do
!	 Upper triangle.
end if
if(lu(n,n)>0.0_8)then
 	solve(n) = solve(n)/lu(n,n)
else
 	solve(n) = 0.0_8
end if
if(n==1) return
do col=n,2,-1
	row=col-1
 	if(lu(row,row)<=0.0_8)then 
 		solve(row) = 0.0_8
 	else
  		solve(row) = (solve(row)-dot_product(lu(row,col:n),solve(col:n)))/lu(row,row)
 	end if
end do 
return
end function solve
