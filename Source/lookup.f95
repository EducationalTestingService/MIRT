!	Look up sorted integer array a.
integer function lookup(x,a)
implicit none
integer, intent(in) :: x,a(:)
!	pivot is pivot point.
integer :: lower,pivot,upper

!	The pivot cannot be less than 1 or greater than size(a).
lower= 1
upper=size(a)
!	The pivot cannot exceed size(a).
pivot=(lower+upper)/2
do

	
		
	
	if(a(pivot)==x) then
		lookup=pivot
		return
    end if
	if(a(pivot)<x) then
        lower=pivot
		pivot=max(lower+1,(lower+upper)/2)
	else
	    upper=pivot
	    pivot=min((lower+upper)/2,upper-1)
    end if

end do
end function lookup




