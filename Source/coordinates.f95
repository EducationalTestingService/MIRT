!	Get coordinates for quadrature points and weights.
!	dimpoints contains number of points per dimension.
!	coord is the coordinate vectors.
!	cweight is the correspond weight multipliers.
subroutine coordinates(dimpoints,coord,cweight)
implicit none
interface
    subroutine readcoordinates(coords,cw)
        implicit none
        integer,intent(inout)::coords(:,:),cw(:)
    end subroutine readcoordinates
end interface
integer,intent(in)::dimpoints(:)
integer,intent(out)::coord(:,:)
real(kind=8),intent(out)::cweight(:)
!	dimno counts dimensions.
!	gridsize is the size of the complete grid.
!	io is the input error flag.
!	quad counts points.
integer::counter,dimno,gridsize,io,quad
!	coordinates to read.
!	weights to read.
integer,allocatable::coords(:,:),cw(:)

!	Inital coordinate.
coord(:,1)=1
gridsize=product(dimpoints)
cweight=1.0_8
do quad=2,size(cweight)
	coord(:,quad)=coord(:,quad-1)
	do dimno=1,size(dimpoints)
		coord(dimno,quad)=coord(dimno,quad)+1
		if(coord(dimno,quad)<=dimpoints(dimno))exit
		coord(dimno,quad)=1
	end do
end do
!	Complete case.
if(gridsize==size(cweight))return

!	Incomplete case.
allocate(coords(size(coord,1),size(coord,2)),&
	cw(size(cweight)),stat=io)
if(io/=0) stop "Space for coordinates not allocated."
cw=1.0_8
coords=coord
call readcoordinates(coords,cw)

!	Check for plausible results
do quad=1,size(cweight)
	if(minval(coord(:,quad))<1.or.minval(dimpoints-coord(:,quad))<0&
		.or.cweight(quad)<=0.0_8) stop "Coordinate input invalid."
end do
coord=coords
cweight=cw
return
end subroutine coordinates
subroutine readcoordinates(coords,cw)
    implicit none
    integer,intent(inout)::coords(:,:),cw(:)
    integer::io
    namelist/griddata/coords,cw
    read(*,nml=griddata,iostat=io)
    if(io/=0)stop "Coordinate input invalid."

    return
end subroutine readcoordinates
