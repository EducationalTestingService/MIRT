!	reduce a sorted integer array a to distinct entries.
! 
subroutine reduce(na,a)
implicit none
integer, intent(in out)::a(:)
integer,intent(out)::na
!	i is a counter
integer::i
na=0
do i=1,size(a)
	if(i==1)then
		na=1
		cycle
	end if
	if(a(i)>a(na))then
		na=na+1
		a(na)=a(i)
	end if
end do
return
end subroutine reduce