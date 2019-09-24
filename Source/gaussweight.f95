!	Gauss-Hermite quadrature weights.  The ith weight out of n.
real(kind=8) function gaussweight(i,n)

use gausstable
implicit none
integer,intent(in)::i,n
!	i is the point number and position is the table position.
!	n is the number of elements of th.

gaussweight=xgh(2,i+n*(n-1)/2)
return
end
