!	This function computes the modified Cholesky decomposition
!	of the n by n nonnegative definite symmecholc macholx sym.
!	tolsing is the singularity tolerance.
function chol(sym,tolsing)
implicit none
real(kind=8),intent(in)::sym(:,:),tolsing
real(kind=8)::chol(size(sym,1),size(sym,1))
real(kind=8)::minsize
integer::row,col
do row=1,size(sym,1)
	minsize = tolsing*sym(row,row)
	chol(1,row) = sym(row,1)
	if(row>1) then
		if(chol(1,1)>0.0_8) chol(row,1) = chol(1,row)/chol(1,1)
	 	do col=2,row
	 		if(col<row.and.chol(col,col)<=0.0)cycle
	 		chol(col,row) = sym(row,col)-dot_product(chol(row,1:col-1),chol(1:col-1,col))
	 		if(col<row) chol(row,col) = chol(col,row)/chol(col,col)
	 	end do
	end if
	if(chol(row,row)>minsize)cycle
	
	do col=1,size(sym,1)
		if(row/=col) then
			chol(row,col) = 0.0_8
			chol(col,row) = 0.0_8
		end if
	end do
end do
return
end function chol
