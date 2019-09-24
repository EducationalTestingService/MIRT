
!	Look up array a of integer pairs.
integer function lookup8(x,a)
implicit none
integer, intent(in) :: x(2),a(:,:)

!	pivot is pivot point.
integer :: lower,pivot,upper

!	The pivot cannot be less than 1 or greater than size(a).
lower= 1
upper=size(a,2)
pivot=(lower+upper)/2
!	The pivot cannot exceed size(a,2).
do

    if  (a(1,pivot)==x(1).and.a(2,pivot)==x(2)) then
        lookup8=pivot
        return
    end if
    if(a(1,pivot)<x(1).or.(a(1,pivot)==x(1).and.a(2,pivot)<x(2)))then
        lower=pivot
        pivot=max(lower+1,(lower+upper)/2)
    else
        upper=pivot
        pivot=min((lower+upper)/2,upper-1)
    end if
	
	

	
end do
end function lookup8





