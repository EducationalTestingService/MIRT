
!	Sort array a of integer pairs in ascending order.
recursive subroutine sort8(a)
implicit none
interface
	subroutine partition8(a,pivot)
!	Array to sort is a.
	  integer, intent(in out)::a(:,:)
!	Pivot is the pivot point.
	  integer, intent(out) :: pivot
	end subroutine partition8
end interface
  integer, intent(in out) :: a(:,:)
!	pivot is pivot point.
  integer :: pivot
integer::junk
  if(size(a,2) > 1) then
	
!	partition is used to select pivot point
	call partition8(a, pivot)
	
	
!	Sort part of a below pivot.
	call sort8(a(:,:pivot-1))
	
!	Sort part of a above pivot.
	call sort8(a(:,pivot:))
	
	
  endif
end subroutine sort8
!	Find pivot point.
subroutine partition8(a, pivot)
!	Array to sort is a.
  integer, intent(in out)::a(:,:)
!	Pivot is the pivot point.
  integer, intent(out) :: pivot
!	lower is lower bound and upper is upper bound for pivot.
!	a values below pivot must not exceed x.  a values above
!	pivot must not be less than x.
  integer :: lower, upper
!	The pivot value is x.  Hold is used to store value to be exchanged.
  integer :: hold(2),x(2)
  x=a(:,size(a,2)/2)
!	The pivot cannot be less than 0.
  lower= 0
!	The pivot cannot exceed size(a).
  upper= size(a,2) + 1

  do
!	Drop upper by 1 for new upper bound on pivot.
     upper = upper-1
     do
!	If a(upper)<=x, then pivot does not exceed upper.
        if (a(1,upper) < x(1).or.(a(1,upper)==x(1).and.a(2,upper)<=x(2))) exit
!	Drop upper by 1 because pivot must be less than upper.
        upper = upper-1
     end do
!	Now check on lower bound.
     lower = lower+1
     do
        if (a(1,lower) > x(1).or.(a(1,lower)==x(1).and.a(2,lower)>=x(2))) exit
!	Increase lower by 1.
        lower = lower+1
     end do
     if (lower < upper) then
!	Exchange a(lower) and a(upper).
        hold = a(:,lower)
        a(:,lower) = a(:,upper)
        a(:,upper) = hold
     else
		if (lower == upper) then
			pivot = upper+1
			return
		else
			pivot = lower
			return
		endif
     endif
  end do

end subroutine partition8




