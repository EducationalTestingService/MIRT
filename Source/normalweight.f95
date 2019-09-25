!	Normal weights for even spacing with centering at 0.
!	pspread is specified as the range of the quadrature points.
!		It is computed automatically  if the input is not positive.
!	Quadrature points are returned in points.
!	Quadrature weights are returned in weights.
subroutine normalweight(pspread,points,weights)
implicit none

real(kind=8),intent(in)::pspread
real(kind=8),intent(out)::points(:),weights(:)

!	i counts points.

integer::i
!	ave is the average of the integers 1 to size(points)
!	b is the multiplier.
!	base is the raw grid.




real(kind=8)::ave,b,base(size(points)),div

ave=(1.0_8+size(points))/2.0_8
b=0.0_8
do i=1,size(points)
	base(i)=i-ave
end do
!	If pspread is positive, then points are determined from spread.
if(pspread>0.0_8) then
	b=pspread/(size(points)-1)
	
else
    if(size(points)>0) b=2.0_8/(size(points)-1.0_8)**0.333333333333333_8
end if
points=b*base
weights=exp(-points*points/2.0_8)
div=sum(weights)
weights=weights/div
return
end subroutine normalweight


