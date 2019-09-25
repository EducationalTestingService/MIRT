!	reduce a sorted array a of integer pairs to distinct entries.
! 
subroutine reduce8(na,a)
implicit none
integer, intent(in out)::a(:,:)
integer,intent(out)::na
!	i is a counter
integer::i
na=0
do i=1,size(a,2)
	if(i==1)then
		na=1
		cycle
	end if
	if(a(1,i)>a(1,na).or.(a(1,i)==a(1,na).and.a(2,i)>a(2,na)))then
		na=na+1
		a(:,na)=a(:,i)
	end if
end do
return
end subroutine reduce8

