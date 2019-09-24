!	Gauss-Hermite quadrature points.  The ith point out of n.
real(kind=8) function gausspoint(i,n)
use gausstable
implicit none
integer,intent(in)::i,n
!	i is the point number and position is the table position.
!	n is the number of elements of th.
gausspoint=xgh(1,i+n*(n-1)/2)
return
end
