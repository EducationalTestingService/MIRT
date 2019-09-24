!	Check on number of points for item response functions and whether quadrature points are used.
!	quadpoint provides regular quadrature points.
!	irfcount is the desired number of points for item response functions.  The default is the number of quadrature points for
!	quadrature.
!	custompoints is only true if special points are to be read.
subroutine getthetasize(quadpoint,irfcount,custompoints)
implicit none
integer,intent(out)::irfcount
logical,intent(out)::custompoints
real(kind=8),intent(in)::quadpoint(:,:)
!	Error flag for input.
integer::io

namelist/irfspecs/irfcount,custompoints
irfcount=size(quadpoint,2)
custompoints=.FALSE.
read(*,nml=irfspecs,iostat=io)
if(io/=0) stop 'Failure to read specifications for points for item response functions.'
if(irfcount/=size(quadpoint,2)) custompoints=.TRUE.
end subroutine getthetasize