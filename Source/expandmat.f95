!	convert an n by n symmetric matrix from compact to regular form.
function expandmat(n,compmat)
implicit none
integer,intent(in)::n
real(kind=8),intent(in)::compmat(:)
real(kind=8)::expandmat(n,n)
integer::row,position
position=1
do row=1,n
	expandmat(row,1:row)=compmat(position:position+row-1)
	if(row>1)then
		expandmat(row,1:row-1)=expandmat(row,1:row-1)/2.0_8
		expandmat(1:row-1,row)=expandmat(row,1:row-1)
	end if
	position=position+row
end do
return
end function expandmat